"""Review tests of the soft exit (PR #617): interleaved replica exchange over
a socket, the Fortran driver, the MPI transport, and the watchdog with a real
simulation run by bin/i-pi.

Every run is a subprocess in a scratch folder. The criterion of all tests
(`assert_continuous`) is that, after restarting from RESTART and finishing,
every output has each step exactly once, in order, with no partial row.
"""

import importlib.util
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import uuid

import pytest

ROOT = Path(__file__).resolve().parents[3]
IPI = ROOT / "bin" / "i-pi"
PY_DRIVER = ROOT / "drivers" / "py" / "driver.py"
F90_DRIVER = ROOT / "bin" / "i-pi-driver"

# watchdog settings used in the subprocesses (the module constant is 10 s)
GRACE = 1.0

# runs bin/i-pi with other watchdog constants, without editing the sources
WRAPPER = """
import os, runpy, sys
import ipi.utils.softexit as softexit_module
softexit_module.GRACE = float(os.environ.get("IPI_TEST_GRACE", softexit_module.GRACE))
softexit_module.SOFTEXITLATENCY = float(
    os.environ.get("IPI_TEST_LATENCY", softexit_module.SOFTEXITLATENCY)
)
sys.argv = [{ipi!r}] + sys.argv[1:]
runpy.run_path({ipi!r}, run_name="__main__")
"""

# a harmonic potential that can be slow (driver), or hang / be slow once (direct)
PES = """
import os, signal, time
from pathlib import Path
from ipi.pes.harmonic import Harmonic_driver

__DRIVER_NAME__ = "slowharmonic"
__DRIVER_CLASS__ = "SlowHarmonic_driver"


class SlowHarmonic_driver(Harmonic_driver):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("k1", 0.1)
        super().__init__(*args, **kwargs)
        self.ncalls = 0
        self.sleep = float(os.environ.get("IPI_TEST_SLEEP", 0))
        self.special = int(os.environ.get("IPI_TEST_SPECIAL_CALL", -1))
        self.special_sleep = float(os.environ.get("IPI_TEST_SPECIAL_SLEEP", 0))
        self.block = os.environ.get("IPI_TEST_BLOCK_SIGNALS", "") == "1"

    def compute_structure(self, cell, pos):
        if self.sleep > 0:
            time.sleep(self.sleep)
        if self.ncalls == self.special:
            if self.block:
                # as when the main thread is stuck in native code
                signal.pthread_sigmask(
                    signal.SIG_BLOCK, [signal.SIGTERM, signal.SIGINT]
                )
            Path("sleeping").touch()
            time.sleep(self.special_sleep)
            Path("awake").touch()
        self.ncalls += 1
        return super().compute_structure(cell, pos)

    def compute(self, cell, pos):
        if isinstance(cell, list):
            return [self.compute_structure(c, p) for c, p in zip(cell, pos)]
        return self.compute_structure(cell, pos)
"""


def system_xml(prefix="", temperature=300, natoms=1, nbeads=1, scale=1.0):
    q = []
    for b in range(nbeads):
        for a in range(natoms):
            q += [(0.1 + 0.05 * a + 0.01 * b) * scale, 0.02 * a, -0.03 * b]
    p = [1.0, 0.2, -0.1] * (natoms * nbeads)
    prefix = f' prefix="{prefix}"' if prefix else ""
    return f"""<system{prefix}>
    <forces><force forcefield="ff"/></forces>
    <ensemble><temperature units="kelvin">{temperature}</temperature></ensemble>
    <motion mode="dynamics">
      <fixcom>False</fixcom>
      <dynamics mode="nvt">
        <timestep units="femtosecond">0.5</timestep>
        <thermostat mode="langevin"><tau units="femtosecond">20</tau></thermostat>
      </dynamics>
    </motion>
    <beads natoms="{natoms}" nbeads="{nbeads}">
      <q shape="({nbeads},{3 * natoms})">{q}</q>
      <p shape="({nbeads},{3 * natoms})">{p}</p>
      <m shape="({natoms})">{[1837.36223469] * natoms}</m>
      <names shape="({natoms})">[{", ".join(["H"] * natoms)}]</names>
    </beads>
    <cell shape="(3,3)">[18.897261,0,0,0,18.897261,0,0,0,18.897261]</cell>
  </system>"""


def simulation_xml(forcefield, systems, nsteps, extra="", interleave="true"):
    return f"""<simulation verbosity="low" interleave="{interleave}" safe_stride="1">
  <output prefix="sim">
    <properties stride="1" filename="out">[ step, potential, conserved ]</properties>
    <trajectory stride="2" filename="pos">positions</trajectory>
  </output>
  <total_steps>{nsteps}</total_steps>
  <prng><seed>31415</seed></prng>
  {forcefield}
  {systems}
  {extra}
</simulation>"""


def ffsocket(address):
    return f"""<ffsocket name="ff" mode="unix">
    <address>{address}</address><latency>0.001</latency><timeout>60</timeout>
  </ffsocket>"""


def ffdirect(pes_path):
    return f"""<ffdirect name="ff">
    <pes>custom</pes><pes_path>{pes_path}</pes_path>
  </ffdirect>"""


class Processes:
    """Starts subprocesses in a scratch folder, and kills what is left."""

    def __init__(self, folder):
        self.folder = Path(folder)
        self.started = []
        self.addresses = []
        self.nlogs = 0
        (self.folder / "pes.py").write_text(PES)
        (self.folder / "run_ipi.py").write_text(WRAPPER.format(ipi=str(IPI)))

    def env(self, **variables):
        env = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONUNBUFFERED="1")
        env["PATH"] = os.pathsep.join(
            [os.path.dirname(sys.executable), str(ROOT / "bin"), env.get("PATH", "")]
        )
        for name in list(env):
            if name.startswith("IPI_TEST_"):
                del env[name]
        env.update({name: str(value) for name, value in variables.items()})
        return env

    def address(self):
        address = "rv" + uuid.uuid4().hex[:12]
        self.addresses.append(address)
        return address

    def start(self, command, name, **variables):
        self.nlogs += 1
        logfile = self.folder / f"{self.nlogs:02d}_{name}.log"
        log = open(logfile, "w")
        process = subprocess.Popen(
            [str(c) for c in command],
            cwd=self.folder,
            env=self.env(**variables),
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        log.close()
        process.logfile = logfile
        self.started.append(process)
        return process

    def ipi(self, input_file, name="ipi", **variables):
        return self.start(
            [sys.executable, self.folder / "run_ipi.py", input_file], name, **variables
        )

    def py_driver(self, address, name="driver", **variables):
        command = [sys.executable, PY_DRIVER, "-u", "-a", address]
        command += ["-m", "custom", "-P", self.folder / "pes.py"]
        return self.start(command, name, **variables)

    def descendants(self):
        """The processes that run in the scratch folder."""

        pids = []
        for entry in os.listdir("/proc"):
            if not entry.isdigit() or int(entry) == os.getpid():
                continue
            try:
                cwd = os.readlink(f"/proc/{entry}/cwd")
                with open(f"/proc/{entry}/cmdline", "rb") as stream:
                    cmdline = stream.read().replace(b"\0", b" ").decode()
            except OSError:
                continue
            if cwd == str(self.folder):
                pids.append((int(entry), cmdline))
        return pids

    def cleanup(self):
        for process in self.started:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except OSError:
                    pass
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
        for pid, _ in self.descendants():
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass
        for address in self.addresses:
            Path("/tmp/ipi_" + address).unlink(missing_ok=True)


@pytest.fixture
def run(tmp_path):
    processes = Processes(tmp_path)
    yield processes
    processes.cleanup()


def wait_for(condition, timeout=30, message="condition"):
    deadline = time.monotonic() + timeout
    while not condition():
        if time.monotonic() > deadline:
            raise TimeoutError(f"Timed out waiting for {message}")
        time.sleep(0.005)


def wait_exit(process, timeout):
    """Waits for a process, and returns the time it took."""

    start = time.monotonic()
    try:
        process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        pytest.fail(
            f"{process.args[:3]} still runs after {timeout} s\n"
            + process.logfile.read_text()[-3000:]
        )
    return time.monotonic() - start


def property_steps(filename):
    """The steps of the rows of a properties file, checking they are complete."""

    text = Path(filename).read_text()
    assert text == "" or text.endswith("\n"), f"{filename}: partial last row"
    steps = []
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        assert len(fields) == 3, f"{filename}: partial row {line!r}"
        [float(field) for field in fields]
        assert float(fields[0]) == int(float(fields[0])), f"{filename}: {line!r}"
        steps.append(int(float(fields[0])))
    return steps


def trajectory_steps(filename, natoms):
    """The steps of the frames of an xyz file, checking they are complete."""

    if not Path(filename).exists():
        return []
    text = Path(filename).read_text()
    assert text == "" or text.endswith("\n"), f"{filename}: partial last line"
    lines = text.splitlines()
    size = natoms + 2
    assert len(lines) % size == 0, f"{filename}: partial frame ({len(lines)} lines)"
    steps = []
    for i in range(0, len(lines), size):
        assert lines[i].split() == [str(natoms)], f"{filename}: line {i}: {lines[i]!r}"
        assert "Step:" in lines[i + 1], f"{filename}: line {i + 1}: {lines[i + 1]!r}"
        steps.append(int(lines[i + 1].split("Step:")[1].split()[0]))
        for line in lines[i + 2 : i + size]:
            fields = line.split()
            assert len(fields) == 4, f"{filename}: partial atom line {line!r}"
            [float(field) for field in fields[1:]]
    return steps


def restart_step(folder):
    """The step RESTART continues from (read without importing i-PI state)."""

    import xml.etree.ElementTree as ET

    restart = Path(folder) / "RESTART"
    assert restart.exists(), "RESTART was not written"
    return int(ET.parse(restart).getroot().find("step").text)


def assert_outputs(folder, last, natoms=1, nbeads=1, prefixes=("",)):
    """Every output has each step up to `last` exactly once, in order."""

    for prefix in prefixes:
        name = str(Path(folder) / (prefix + "_sim" if prefix else "sim"))
        assert property_steps(name + ".out") == list(range(last + 1)), name + ".out"
        for bead in range(nbeads):
            trajectory = f"{name}.pos_{bead}.xyz"
            assert trajectory_steps(trajectory, natoms) == list(
                range(0, last + 1, 2)
            ), trajectory


def assert_continuous(folder, nsteps, **kwargs):
    """Criterion C, on the outputs of a run that was restarted and finished."""

    assert_outputs(folder, nsteps, **kwargs)


def nrows(filename):
    try:
        with open(filename) as stream:
            return sum(1 for line in stream if not line.startswith("#"))
    except OSError:
        return 0


def assert_clean_log(process):
    log = process.logfile.read_text()
    for error in ("Traceback", "GreenletExit", "RuntimeError", "Broken pipe"):
        assert error not in log, f"{error} in {process.logfile}:\n{log[-3000:]}"
    return log


# 16


def test_interleaved_remd_socket_exit(run):
    """SIGTERM while the tasks of four replicas wait for two drivers."""

    nsteps, replicas = 30, ["r0", "r1", "r2", "r3"]
    folder, address = run.folder, run.address()
    systems = "".join(
        system_xml(prefix=prefix, temperature=300 + 40 * i, natoms=2, scale=1 + 0.1 * i)
        for i, prefix in enumerate(replicas)
    )
    remd = '<smotion mode="remd"><remd><stride>1</stride></remd></smotion>'
    (folder / "input.xml").write_text(
        simulation_xml(ffsocket(address), systems, nsteps, extra=remd)
    )

    server = run.ipi("input.xml")
    wait_for(lambda: os.path.exists("/tmp/ipi_" + address), message="the socket")
    drivers = [
        run.py_driver(address, f"driver{i}", IPI_TEST_SLEEP=0.05) for i in range(2)
    ]
    wait_for(
        lambda: nrows(folder / "r3_sim.out") >= 8 or server.poll() is not None,
        timeout=60,
        message="8 steps",
    )
    assert server.poll() is None, server.logfile.read_text()[-3000:]
    # two requests are running and two are queued most of the time
    time.sleep(0.02)
    server.send_signal(signal.SIGTERM)
    elapsed = wait_exit(server, 15)
    for driver in drivers:
        wait_exit(driver, 10)

    log = assert_clean_log(server)
    assert server.returncode == 0, log[-3000:]
    assert "Restartable as is: YES" in log
    assert elapsed < 5
    for driver in drivers:
        assert driver.returncode == 0, driver.logfile.read_text()[-3000:]
        assert "Traceback" not in driver.logfile.read_text()
    assert not os.path.exists("/tmp/ipi_" + address)

    step = restart_step(folder)
    assert 7 <= step < nsteps
    # all the replicas are at the step of RESTART
    assert_outputs(folder, step, natoms=2, prefixes=replicas)
    swapfiles = [p for p in folder.iterdir() if "PARATEMP" in p.name]

    server = run.ipi("RESTART", "restart")
    wait_for(lambda: os.path.exists("/tmp/ipi_" + address), message="the socket")
    drivers = [run.py_driver(address, f"driver{i}") for i in range(2)]
    wait_exit(server, 60)
    for driver in drivers:
        wait_exit(driver, 10)
    log = assert_clean_log(server)
    assert server.returncode == 0, log[-3000:]
    for driver in drivers:
        assert driver.returncode == 0, driver.logfile.read_text()[-3000:]

    assert_continuous(folder, nsteps, natoms=2, prefixes=replicas)
    # the swaps are recorded once per step
    for swapfile in swapfiles:
        steps = [
            int(float(line.split()[0]))
            for line in swapfile.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert steps == sorted(set(steps)), f"{swapfile}: {steps}"


# 17


@pytest.mark.skipif(not F90_DRIVER.exists(), reason="needs the compiled i-pi-driver")
def test_fortran_driver_survives_close(run):
    """i-PI exits while the Fortran driver computes: the driver finds the
    connection closed, and ends as at the end of a run."""

    nsteps = 12
    folder, address = run.folder, run.address()
    (folder / "input.xml").write_text(
        simulation_xml(ffsocket(address), system_xml(), nsteps)
    )

    def driver(sleep, name):
        command = [F90_DRIVER, "-u", "-a", address, "-m", "gas", "-o", sleep]
        return run.start(command, name)

    server = run.ipi("input.xml")
    wait_for(lambda: os.path.exists("/tmp/ipi_" + address), message="the socket")
    client = driver("0.3", "f90")
    wait_for(
        lambda: nrows(folder / "sim.out") >= 3 or server.poll() is not None,
        timeout=60,
        message="3 steps",
    )
    assert server.poll() is None, server.logfile.read_text()[-3000:]
    # in the middle of the next evaluation, which takes 0.3 s
    time.sleep(0.1)
    server.send_signal(signal.SIGTERM)
    elapsed = wait_exit(server, 15)
    wait_exit(client, 10)

    log = assert_clean_log(server)
    assert server.returncode == 0, log[-3000:]
    assert elapsed < 5
    driver_log = client.logfile.read_text()
    assert client.returncode != -signal.SIGPIPE, driver_log[-2000:]
    assert client.returncode == 0, (client.returncode, driver_log[-2000:])

    step = restart_step(folder)
    assert 2 <= step < nsteps
    assert_outputs(folder, step)

    server = run.ipi("RESTART", "restart")
    wait_for(lambda: os.path.exists("/tmp/ipi_" + address), message="the socket")
    client = driver("0", "f90")
    wait_exit(server, 60)
    wait_exit(client, 10)
    assert server.returncode == 0, server.logfile.read_text()[-3000:]
    assert client.returncode == 0, client.logfile.read_text()[-2000:]
    assert_continuous(folder, nsteps)


# 18


def mpirun_command(folder, input_file):
    mpirun = shutil.which("mpirun", path=os.path.dirname(sys.executable)) or (
        shutil.which("mpirun")
    )
    version = subprocess.run(
        [mpirun, "--version"], capture_output=True, text=True, timeout=30
    )
    command = [mpirun]
    if "HYDRA" not in version.stdout + version.stderr:
        # Open MPI: allows more ranks than cores
        command.append("--oversubscribe")
    command += ["-n", "1", sys.executable, folder / "run_ipi.py", input_file]
    command += [":", "-n", "2", sys.executable, PY_DRIVER, "--mpi"]
    command += ["-m", "custom", "-P", folder / "pes.py"]
    return command


@pytest.mark.skipif(
    importlib.util.find_spec("mpi4py") is None
    or (
        shutil.which("mpirun", path=os.path.dirname(sys.executable)) is None
        and shutil.which("mpirun") is None
    ),
    reason="needs mpi4py and mpirun",
)
def test_mpi_exit(run):
    """SIGTERM to the i-PI rank of an MPI job: the job ends cleanly, and
    restarts."""

    nsteps = 30
    folder = run.folder
    forcefield = '<ffmpi name="ff" mode="mpi"></ffmpi>'
    (folder / "input.xml").write_text(
        simulation_xml(forcefield, system_xml(nbeads=2), nsteps)
    )

    job = run.start(mpirun_command(folder, "input.xml"), "mpi", IPI_TEST_SLEEP=0.05)
    wait_for(
        lambda: nrows(folder / "sim.out") >= 8 or job.poll() is not None,
        timeout=120,
        message="8 steps",
    )
    assert job.poll() is None, job.logfile.read_text()[-3000:]
    servers = [
        pid
        for pid, cmdline in run.descendants()
        if "run_ipi.py" in cmdline and "mpirun" not in cmdline.split()[0]
        and "hydra" not in cmdline.split()[0]
        and "mpiexec" not in cmdline.split()[0]
    ]
    assert len(servers) == 1, run.descendants()
    time.sleep(0.02)
    os.kill(servers[0], signal.SIGTERM)
    elapsed = wait_exit(job, 60)

    log = job.logfile.read_text()
    assert job.returncode == 0, (job.returncode, log[-3000:])
    for error in ("Traceback", "BAD TERMINATION", "Abort", "RuntimeError"):
        assert error not in log, f"{error} in the log:\n{log[-3000:]}"
    assert "Restartable as is: YES" in log
    assert elapsed < 20
    assert run.descendants() == []

    step = restart_step(folder)
    assert 7 <= step < nsteps
    assert_outputs(folder, step, nbeads=2)

    job = run.start(mpirun_command(folder, "RESTART"), "mpi_restart")
    wait_exit(job, 120)
    log = job.logfile.read_text()
    assert job.returncode == 0, (job.returncode, log[-3000:])
    assert "Traceback" not in log, log[-3000:]
    assert_continuous(folder, nsteps, nbeads=2)


# 19


def direct_input(folder, nsteps):
    (folder / "input.xml").write_text(
        simulation_xml(ffdirect(folder / "pes.py"), system_xml(), nsteps)
    )


def test_watchdog_with_simulation(run):
    """The potential hangs at a step, with the signals blocked in the main
    thread: the watchdog writes the stored checkpoint and ends the process."""

    nsteps, hang_call = 10, 5
    folder = run.folder
    direct_input(folder, nsteps)

    server = run.ipi(
        "input.xml",
        IPI_TEST_GRACE=GRACE,
        IPI_TEST_SPECIAL_CALL=hang_call,
        IPI_TEST_SPECIAL_SLEEP=60,
        IPI_TEST_BLOCK_SIGNALS=1,
    )
    wait_for(
        lambda: (folder / "sleeping").exists() or server.poll() is not None,
        message="the hang",
    )
    assert server.poll() is None, server.logfile.read_text()[-3000:]
    time.sleep(0.1)
    # call 0 is made for the initial outputs, call n in the step n - 1
    assert nrows(folder / "sim.out") == hang_call
    server.send_signal(signal.SIGTERM)
    elapsed = wait_exit(server, GRACE + 15)

    log = server.logfile.read_text()
    assert server.returncode == 1, (server.returncode, log[-3000:])
    assert "has not responded" in log, log[-3000:]
    assert "could not be written" not in log, log[-3000:]
    assert not (folder / "awake").exists()
    assert GRACE * 0.9 < elapsed < GRACE + 5

    # the step that hung is discarded
    assert restart_step(folder) == hang_call - 1
    assert_outputs(folder, hang_call - 1)

    server = run.ipi("RESTART", "restart")
    wait_exit(server, 60)
    assert server.returncode == 0, server.logfile.read_text()[-3000:]
    assert_continuous(folder, nsteps)


# 20

REPETITIONS = 10


def slow_not_hung(run, repetition, offset):
    """One run with a potential that is slow once, for GRACE + offset, and an
    EXIT file made meanwhile. Returns the exit status and what is wrong."""

    nsteps, slow_call = 8, 4
    folder = run.folder / f"rep{repetition}"
    folder.mkdir()
    shutil.copy(run.folder / "pes.py", folder / "pes.py")
    shutil.copy(run.folder / "run_ipi.py", folder / "run_ipi.py")
    direct_input(folder, nsteps)
    here = Processes.__new__(Processes)
    here.__dict__.update(run.__dict__, folder=folder)

    server = here.ipi(
        "input.xml",
        f"rep{repetition}",
        IPI_TEST_GRACE=GRACE,
        IPI_TEST_LATENCY=0.01,
        IPI_TEST_SPECIAL_CALL=slow_call,
        IPI_TEST_SPECIAL_SLEEP=GRACE + offset,
    )
    wait_for(
        lambda: (folder / "sleeping").exists() or server.poll() is not None,
        message="the slow evaluation",
    )
    (folder / "EXIT").touch()
    wait_exit(server, 2 * GRACE + 20)
    (folder / "EXIT").unlink()
    status = server.returncode
    log = server.logfile.read_text()
    kind = f"rep {repetition} (offset {offset:+.3f} s, status {status})"

    problems = []
    try:
        assert status in (0, 1), log[-2000:]
        assert "Traceback" not in log, log[-2000:]
        step = restart_step(folder)
        # the slow step is discarded, or was completed
        assert step in (slow_call - 1, slow_call)
        assert_outputs(folder, step)
    except Exception as err:  # also a RESTART that cannot be parsed
        problems.append(f"{kind}: after the exit: {type(err).__name__}: {err}")
        return status, problems

    server = here.ipi("RESTART", f"rep{repetition}_restart")
    wait_exit(server, 60)
    try:
        assert server.returncode == 0, server.logfile.read_text()[-2000:]
        assert_continuous(folder, nsteps)
    except Exception as err:
        problems.append(f"{kind}: after the restart: {type(err).__name__}: {err}")
    return status, problems


def test_watchdog_slow_not_hung(run):
    """A potential that takes about GRACE, but returns, races with the
    watchdog: RESTART and the outputs agree whatever the exit status."""

    # around the time at which the watchdog acts, and clearly later
    offsets = [-0.010 + 0.005 * i for i in range(REPETITIONS - 2)] + [0.3, 1.0]
    statuses, problems = [], []
    for repetition, offset in enumerate(offsets):
        status, found = slow_not_hung(run, repetition, offset)
        statuses.append(status)
        problems += found
    print("exit statuses:", statuses)
    assert not problems, (
        f"{len(problems)} of {REPETITIONS} repetitions inconsistent "
        f"(statuses {statuses}):\n" + "\n".join(problems)
    )

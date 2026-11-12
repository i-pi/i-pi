"""End-to-end checks of the soft exit, through bin/i-pi and socket clients.

i-PI runs as a program (bin/i-pi, so with the hard exit and the watchdog),
with the Python driver (a slowed-down harmonic PES) or a fake client that
stops replying. Each test interrupts a run, and then checks criterion C:
after restarting from RESTART and finishing, every output has each step
exactly once, in order, with no partial row (and, where the forces are those
of the reference, the same values as an uninterrupted run).
"""

import os
from pathlib import Path
import re
import signal
import socket
import subprocess
import sys
import time
import uuid

import numpy as np
import pytest

from ipi_tests.unit_tests.engine.test_softexit_sockets import hang, message

ROOT = Path(__file__).resolve().parents[3]
IPI = str(ROOT / "bin" / "i-pi")
DRIVER = str(ROOT / "drivers" / "py" / "driver.py")
ENV = dict(
    os.environ,
    PYTHONPATH=str(ROOT),
    PYTHONUNBUFFERED="1",
    PATH=os.pathsep.join(
        [str(Path(sys.executable).parent), str(ROOT / "bin"), os.environ["PATH"]]
    ),
)
K = 0.1  # spring constant of the harmonic PES

SLOW_PES = '''"""A harmonic PES that takes a while."""
import time
from ipi.pes.harmonic import Harmonic_driver

__DRIVER_NAME__ = "slowharm"
__DRIVER_CLASS__ = "Slow_driver"


class Slow_driver(Harmonic_driver):
    def __init__(self, k1, delay=0.05, *args, **kwargs):
        self.delay = float(delay)
        super().__init__(k1, *args, **kwargs)

    def compute_structure(self, cell, pos):
        time.sleep(self.delay)
        return super().compute_structure(cell, pos)
'''


def forcefield_xml(name, address=None, port=None, options=""):
    """A socket forcefield (unix or inet), or the in-process reference."""

    if address is None and port is None:
        return f"""<ffdirect name="{name}"><pes>harmonic</pes>
    <parameters>{{k1: {K}}}</parameters></ffdirect>"""
    if port is None:
        where = f"<address>{address}</address>"
        mode = "unix"
    else:
        where = f"<address>localhost</address><port>{port}</port>"
        mode = "inet"
    return f"""<ffsocket name="{name}" mode="{mode}">{where}
    <latency>0.001</latency><timeout>60</timeout>{options}</ffsocket>"""


def input_xml(forcefields, nsteps, nbeads=1, extra="", outputs=""):
    """A small NVE run of one atom, with properties and bead trajectories."""

    names = re.findall(r'<ff\w+ name="(\w+)"', forcefields)
    forces = "".join(f'<force forcefield="{name}"/>' for name in names)
    q = ",".join("%g,0,0" % (0.1 + 0.01 * b) for b in range(nbeads))
    p = ",".join(["1,0.2,-0.1"] * nbeads)
    return f"""<simulation verbosity="medium">
  <output prefix="sim">
    <properties stride="1" filename="out">[ step, potential ]</properties>
    <trajectory stride="1" filename="pos">positions</trajectory>
    {outputs}
  </output>
  <total_steps>{nsteps}</total_steps>
  {extra}
  <prng><seed>31415</seed></prng>
  {forcefields}
  <system>
    <forces>{forces}</forces>
    <motion mode="dynamics">
      <fixcom>False</fixcom>
      <dynamics mode="nve">
        <timestep units="femtosecond">0.5</timestep>
      </dynamics>
    </motion>
    <ensemble><temperature units="kelvin">300</temperature></ensemble>
    <beads natoms="1" nbeads="{nbeads}">
      <q shape="({nbeads},3)">[{q}]</q>
      <p shape="({nbeads},3)">[{p}]</p>
      <m shape="(1)">[1837.36223469]</m>
      <names shape="(1)">[H]</names>
    </beads>
    <cell shape="(3,3)">
      [18.897261,0,0,0,18.897261,0,0,0,18.897261]
    </cell>
  </system>
</simulation>"""


class Run:
    """The processes of a test, in a scratch directory. Everything that is
    started is killed by `cleanup`."""

    def __init__(self, path):
        self.path = Path(path)
        self.processes = []
        self.addresses = []
        self.nlogs = 0
        (self.path / "slowharm.py").write_text(SLOW_PES)

    def address(self):
        address = "rev617e2e_" + uuid.uuid4().hex[:16]
        self.addresses.append(address)
        return address

    def socket_file(self, address):
        return Path("/tmp/ipi_" + address)

    def spawn(self, command, name):
        self.nlogs += 1
        logfile = self.path / ("log.%02d.%s" % (self.nlogs, name))
        with open(logfile, "w") as log:
            process = subprocess.Popen(
                command,
                cwd=self.path,
                env=ENV,
                stdout=log,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
            )
        process.logfile = logfile
        self.processes.append(process)
        return process

    def ipi(self, input_file, wait_for_log=()):
        """Starts bin/i-pi, and returns once its sockets listen."""

        server = self.spawn([sys.executable, IPI, str(self.path / input_file)], "ipi")
        wait_until(
            lambda: all(self.socket_file(a).exists() for a in self.addresses)
            and all(text in server.logfile.read_text() for text in wait_for_log),
            server,
            what="the sockets of i-PI",
        )
        return server

    def driver(self, address=None, port=None, delay=None):
        """Starts the Python driver, slowed down by `delay` seconds per
        evaluation."""

        command = [sys.executable, DRIVER]
        if port is None:
            command += ["-u", "-a", address]
        else:
            command += ["-a", "localhost", "-p", str(port)]
        if delay is None:
            command += ["-m", "harmonic", "-o", str(K)]
        else:
            command += ["-m", "custom", "-P", str(self.path / "slowharm.py")]
            command += ["-o", "%s,%s" % (K, delay)]
        return self.spawn(command, "driver")

    def logs(self):
        return "\n".join(
            "=== %s (rc %s)\n%s" % (p.logfile.name, p.poll(), p.logfile.read_text())
            for p in self.processes
        )

    def cleanup(self):
        for process in self.processes:
            if process.poll() is None:
                process.kill()
        for process in self.processes:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
        for address in self.addresses:
            self.socket_file(address).unlink(missing_ok=True)

    # --- checks

    def rows(self, filename="sim.out"):
        """The complete rows written so far."""

        path = self.path / filename
        if not path.exists():
            return []
        lines = path.read_text().split("\n")[:-1]  # drops an incomplete row
        return [line for line in lines if not line.startswith("#")]

    def restart_step(self):
        text = (self.path / "RESTART").read_text()
        match = re.search(r"<step>\s*(\d+)\s*</step>", text)
        return int(match.group(1)) if match else 0

    def frames(self, bead):
        """The steps of the frames of the trajectory of a bead, checking
        that each frame is complete."""

        path = self.path / ("sim.pos_%d.xyz" % bead)
        if not path.exists():
            return [], np.zeros((0, 3))
        text = path.read_text()
        assert text == "" or text.endswith("\n"), "partial line in %s" % path.name
        lines = text.splitlines()
        assert len(lines) % 3 == 0, "partial frame in %s" % path.name
        steps, positions = [], []
        for count, comment, atom in zip(lines[0::3], lines[1::3], lines[2::3]):
            assert count.split() == ["1"], (path.name, count)
            steps.append(int(comment.split("Step:")[1].split()[0]))
            fields = atom.split()
            assert len(fields) == 4 and fields[0] == "H", (path.name, atom)
            positions.append([float(x) for x in fields[1:]])
        return steps, np.array(positions)

    def table(self, filename="sim.out", ncolumns=2):
        """The rows of a property file, checking that each is complete."""

        path = self.path / filename
        text = path.read_text()
        assert text.endswith("\n"), "partial row in %s" % filename
        rows = [line for line in text.splitlines() if not line.startswith("#")]
        for row in rows:
            assert len(row.split()) == ncolumns, "partial row in %s: %r" % (
                filename,
                row,
            )
        return np.array([[float(x) for x in row.split()] for row in rows])

    def assert_rolled_back(self, nbeads=1):
        """After the exit, the outputs hold the steps up to that of RESTART,
        and nothing else. Returns that step."""

        step = self.restart_step()
        expected = list(range(step + 1))
        assert list(self.table()[:, 0]) == expected, "sim.out vs RESTART at %d" % step
        for bead in range(nbeads):
            assert self.frames(bead)[0] == expected, "bead %d vs RESTART" % bead
        return step

    def assert_complete(self, nsteps, nbeads=1, reference=None):
        """Criterion C."""

        expected = list(range(nsteps + 1))
        table = self.table()
        assert list(table[:, 0]) == expected, "steps of sim.out"
        if reference is not None:
            np.testing.assert_allclose(
                table[:, 1], reference.table()[:, 1], rtol=1e-7, atol=1e-12
            )
        for bead in range(nbeads):
            steps, positions = self.frames(bead)
            assert steps == expected, "frames of bead %d" % bead
            if reference is not None:
                np.testing.assert_allclose(
                    positions, reference.frames(bead)[1], rtol=2e-5, atol=1e-9
                )


def wait_until(condition, *processes, timeout=30, what="i-PI"):
    """Waits for condition(), failing if one of the processes ends first."""

    deadline = time.monotonic() + timeout
    while not condition():
        for process in processes:
            if process.poll() is not None:
                raise AssertionError(
                    "A process ended (rc %s) while waiting for %s:\n%s"
                    % (process.returncode, what, process.logfile.read_text()[-3000:])
                )
        if time.monotonic() > deadline:
            raise TimeoutError("Timed out waiting for " + what)
        time.sleep(0.005)


def wait_exit(process, timeout, run):
    """Waits for a process to end, and returns how long it took."""

    start = time.monotonic()
    try:
        process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        pytest.fail(
            "%s still runs after %g s\n%s" % (process.logfile.name, timeout, run.logs())
        )
    return time.monotonic() - start


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("localhost", 0))
        return probe.getsockname()[1]


@pytest.fixture
def run(tmp_path):
    run = Run(tmp_path)
    yield run
    run.cleanup()


_references = {}


@pytest.fixture
def reference(tmp_path_factory):
    """Returns the outputs of an uninterrupted run with the in-process PES."""

    def make(nsteps, nbeads=1):
        key = (nsteps, nbeads)
        if key not in _references:
            ref = Run(tmp_path_factory.mktemp("reference"))
            try:
                xml = input_xml(forcefield_xml("harmonic"), nsteps, nbeads)
                (ref.path / "input.xml").write_text(xml)
                server = ref.ipi("input.xml")
                wait_exit(server, 60, ref)
                assert server.returncode == 0, ref.logs()
                ref.assert_complete(nsteps, nbeads)
            finally:
                ref.cleanup()
            _references[key] = ref
        return _references[key]

    return make


def restart_and_finish(run, addresses=(), ports=(), ndrivers=1, edit=None):
    """Restarts from RESTART with new (fast) drivers, and checks that i-PI
    and the drivers end cleanly."""

    (run.path / "EXIT").unlink(missing_ok=True)
    if edit is not None:
        restart = run.path / "RESTART"
        restart.write_text(edit(restart.read_text()))
    server = run.ipi("RESTART")
    drivers = [run.driver(address=a) for a in addresses for _ in range(ndrivers)]
    drivers += [run.driver(port=p) for p in ports]
    wait_exit(server, 60, run)
    assert server.returncode == 0, run.logs()
    for driver in drivers:
        wait_exit(driver, 15, run)
        assert driver.returncode == 0, run.logs()
    for address in addresses:
        assert not run.socket_file(address).exists()


# --- 8


@pytest.mark.parametrize("trigger", ["SIGINT", "SIGTERM", "EXIT", "total_time"])
def test_exit_triggers_subprocess(run, reference, trigger):
    """Each way of asking for a soft exit ends bin/i-pi and its driver
    cleanly, with a RESTART that continues the outputs."""

    nsteps = 60
    address = run.address()
    extra = "<total_time>2.5</total_time>" if trigger == "total_time" else ""
    xml = input_xml(forcefield_xml("harmonic", address), nsteps, extra=extra)
    (run.path / "input.xml").write_text(xml)

    server = run.ipi("input.xml")
    driver = run.driver(address=address, delay=0.05)
    if trigger != "total_time":
        wait_until(lambda: len(run.rows()) >= 10, server, driver, what="10 steps")
    if trigger in ("SIGINT", "SIGTERM"):
        server.send_signal(getattr(signal, trigger))
        limit = 2.0
    elif trigger == "EXIT":
        (run.path / "EXIT").touch()
        limit = 3.5  # the file is looked for once per second
    else:
        limit = 15.0
    elapsed = wait_exit(server, 20, run)
    wait_exit(driver, 10, run)
    log = server.logfile.read_text()

    assert server.returncode == 0, run.logs()
    assert driver.returncode == 0, run.logs()
    assert elapsed < limit, "i-PI took %.2f s to exit" % elapsed
    assert (run.path / "RESTART").exists()
    assert not run.socket_file(address).exists(), "socket file left behind"
    step = run.assert_rolled_back()
    assert step < nsteps
    if trigger != "total_time":
        assert step >= 9

    # the wall-clock limit of RESTART would end the run again
    edit = None
    if trigger == "total_time":

        def edit(text):
            assert "<total_time>" in text
            return re.sub(r"<total_time>.*?</total_time>", "", text, flags=re.S)

    restart_and_finish(run, [address], edit=edit)
    run.assert_complete(nsteps, reference=reference(nsteps))
    # checked last, so that a wrong message does not hide the checks above
    assert "Restartable as is: YES" in log, log[-1500:]


# --- 9


SCALEDCOORDS = (
    '<properties stride="4" filename="sc">'
    "[ step, scaledcoords(fd_delta=5e-3) ]</properties>"
)


@pytest.mark.parametrize("where", ["step", "outputs"])
def test_exit_during_force_evaluation(run, where):
    """SIGTERM while the only client hangs: in the middle of a step, and in
    the middle of the outputs (an estimator that needs more forces, after
    other outputs of the same step were written)."""

    nsteps = 12
    address = run.address()
    if where == "step":
        xml = input_xml(forcefield_xml("harmonic", address), nsteps)
        # the initial forces and those of 5 steps: hangs in the sixth step
        nreplies, expected = 6, 5
    else:
        xml = input_xml(
            forcefield_xml("harmonic", address),
            nsteps,
            outputs=SCALEDCOORDS,
        )
        # the initial forces, the NSC extra ones of the initial outputs, and
        # 4 steps: hangs in the outputs of step 4, after sim.out and the
        # trajectory were written. The state was last stored at step 3.
        nreplies, expected = 1 + NSC + 4, 3
    (run.path / "input.xml").write_text(xml)

    server = run.ipi("input.xml")
    peer = None
    try:
        peer = hang(address, nreplies=nreplies)
        time.sleep(0.3)
        assert server.poll() is None, run.logs()
        server.send_signal(signal.SIGTERM)
        elapsed = wait_exit(server, 20, run)
        assert server.returncode == 0, run.logs()
        assert elapsed < 2.0, "i-PI took %.2f s to exit" % elapsed
        assert message("exit") in peer.recv(4096)
    finally:
        if peer is not None:
            peer.close()

    assert "Restartable as is: YES" in server.logfile.read_text()
    assert not run.socket_file(address).exists(), "socket file left behind"
    assert run.restart_step() == expected
    run.assert_rolled_back()
    if where == "outputs":
        # only the initial row: the one of step 4 was not complete
        assert list(run.table("sim.sc", 3)[:, 0]) == [0]

    restart_and_finish(run, [address])
    run.assert_complete(nsteps)
    if where == "outputs":
        assert list(run.table("sim.sc", 3)[:, 0]) == [0, 4, 8, 12]


# number of force evaluations that scaledcoords adds to an output (1 bead)
NSC = 2


# --- 12


def test_exit_with_multiple_clients(run, reference):
    """SIGTERM while four drivers compute the forces of eight beads."""

    nsteps, nbeads, ndrivers = 12, 8, 4
    address = run.address()
    xml = input_xml(forcefield_xml("harmonic", address), nsteps, nbeads)
    (run.path / "input.xml").write_text(xml)

    server = run.ipi("input.xml")
    drivers = [run.driver(address=address, delay=0.15) for _ in range(ndrivers)]
    wait_until(lambda: len(run.rows()) >= 4, server, *drivers, what="3 steps")
    # all the drivers have connected and compute
    wait_until(
        lambda: server.logfile.read_text().count("Handshaking was successful")
        >= ndrivers,
        server,
        what="the handshakes",
    )
    time.sleep(0.2)
    server.send_signal(signal.SIGTERM)
    elapsed = wait_exit(server, 20, run)
    for driver in drivers:
        wait_exit(driver, 10, run)

    assert server.returncode == 0, run.logs()
    assert [driver.returncode for driver in drivers] == [0] * ndrivers, run.logs()
    assert elapsed < 2.0, "i-PI took %.2f s to exit" % elapsed
    assert "Restartable as is: YES" in server.logfile.read_text()
    assert not run.socket_file(address).exists(), "socket file left behind"
    step = run.assert_rolled_back(nbeads)
    assert 3 <= step < nsteps

    # one driver: more could find the run finished before they connect
    restart_and_finish(run, [address])
    run.assert_complete(nsteps, nbeads, reference=reference(nsteps, nbeads))


# --- 13


def test_exit_with_two_sockets_one_hung(run):
    """SIGTERM with a unix and an inet socket, the client of the first one
    hanging: both sockets are closed and both clients are told to exit."""

    nsteps = 12
    address = run.address()
    port = free_port()
    forcefields = forcefield_xml("harmonic", address) + forcefield_xml(
        "second", port=port
    )
    (run.path / "input.xml").write_text(input_xml(forcefields, nsteps))

    server = run.ipi("input.xml", wait_for_log=["Created inet socket"])
    driver = run.driver(port=port)
    peer = None
    try:
        peer = hang(address, nreplies=6)
        time.sleep(0.3)
        assert server.poll() is None, run.logs()
        server.send_signal(signal.SIGTERM)
        elapsed = wait_exit(server, 20, run)
        wait_exit(driver, 10, run)
        assert server.returncode == 0, run.logs()
        assert driver.returncode == 0, run.logs()
        assert elapsed < 2.0, "i-PI took %.2f s to exit" % elapsed
        assert message("exit") in peer.recv(4096)
    finally:
        if peer is not None:
            peer.close()

    assert "Restartable as is: YES" in server.logfile.read_text()
    assert not run.socket_file(address).exists(), "socket file left behind"
    # nothing listens on the port any more, and it can be bound at once
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        assert probe.connect_ex(("localhost", port)) != 0
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        probe.bind(("localhost", port))
    assert run.restart_step() == 5
    run.assert_rolled_back()

    # i-PI itself listens on the same port right away
    restart_and_finish(run, [address], ports=[port])
    run.assert_complete(nsteps)


# --- 14


def test_exit_on_disconnect(run, reference):
    """With exit_on_disconnect, killing the driver in the middle of a step
    makes a soft exit, from which a new driver continues."""

    nsteps = 40
    address = run.address()
    options = "<exit_on_disconnect>True</exit_on_disconnect>"
    xml = input_xml(forcefield_xml("harmonic", address, options=options), nsteps)
    (run.path / "input.xml").write_text(xml)

    server = run.ipi("input.xml")
    driver = run.driver(address=address, delay=0.05)
    wait_until(lambda: len(run.rows()) >= 8, server, driver, what="8 steps")
    driver.kill()
    driver.wait(timeout=5)
    elapsed = wait_exit(server, 20, run)
    log = server.logfile.read_text()

    assert server.returncode == 0, run.logs()
    assert elapsed < 2.0, "i-PI took %.2f s to exit" % elapsed
    assert "Client shutdown" in log, log[-3000:]
    assert "Restartable as is: YES" in log, log[-3000:]
    assert not run.socket_file(address).exists(), "socket file left behind"
    step = run.assert_rolled_back()
    assert 7 <= step < nsteps

    restart_and_finish(run, [address])
    run.assert_complete(nsteps, reference=reference(nsteps))


# --- 15


def test_client_dies_and_is_replaced(run, reference):
    """Without exit_on_disconnect, the request of a driver that is killed is
    computed by the other one, and the run finishes with the right results."""

    nsteps, nbeads = 40, 2
    address = run.address()
    xml = input_xml(forcefield_xml("harmonic", address), nsteps, nbeads)
    (run.path / "input.xml").write_text(xml)

    server = run.ipi("input.xml")
    drivers = [run.driver(address=address, delay=0.05) for _ in range(2)]
    wait_until(
        lambda: server.logfile.read_text().count("Handshaking was successful") >= 2
        and len(run.rows()) >= 8,
        server,
        *drivers,
        what="8 steps with two drivers",
    )
    drivers[0].kill()
    drivers[0].wait(timeout=5)
    killed_at = len(run.rows())
    wait_exit(server, 60, run)
    wait_exit(drivers[1], 10, run)
    log = server.logfile.read_text()

    assert killed_at < nsteps, "the driver was killed too late"
    assert server.returncode == 0, run.logs()
    assert drivers[1].returncode == 0, run.logs()
    assert "Restartable as is: NO" in log and "success" in log, log[-3000:]
    assert not run.socket_file(address).exists(), "socket file left behind"
    run.assert_complete(nsteps, nbeads, reference=reference(nsteps, nbeads))

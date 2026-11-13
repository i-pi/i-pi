#!/usr/bin/env python
"""Benchmark of the time i-PI spends per MD step, excluding force evaluation.

With `--forcefield dummy` (the default) each case uses a direct forcefield
with the `dummy` PES (zero forces), so the whole force request/gather pipeline
runs, but there is no socket communication and no actual force computation.
With `--forcefield none` the forcefield has zero weight and is never
evaluated, so only the integration of the equations of motion is timed.
With `--forcefield metatomic` the forces come from a metatomic model
(`--model`, evaluated on `--model-device`) through a direct forcefield, which
evaluates all the beads in one batch. The structure is then a supercell of
the 32-water box or of the 32-atom nickel box (`--structure`), rather than
hydrogen atoms at random positions, and the sizes are given as supercells.
For direct forcefields, the time taken by the driver alone on the same
structures is also reported.
The timed loop calls `Simulation.run_step` directly, so outputs and
checkpointing are excluded. Every (case, size) combination runs in a
separate subprocess.

Usage:
    python benchmark.py [-c nve,pimd_npt] [-n 100,1000] [-o results.json]
    python benchmark.py --forcefield none -o noff.json
    python benchmark.py --backend torch --device cuda -o torch-cuda.json
    python benchmark.py --forcefield metatomic --model pet-mad-xs-v1.6.0.pt \\
        --model-device cuda --structure water -c nve,pimd4_nvt -s 1,2,3
    python benchmark.py --compare reference.json new.json

The `ipi` package is imported from the source tree containing this script, or
from `--ipi-root` if given, which allows benchmarking different source trees
with the same script.
"""

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from importlib import metadata as importlib_metadata

import numpy as np

LANGEVIN = "<thermostat mode='langevin'><tau units='femtosecond'>100</tau></thermostat>"
SVR = "<thermostat mode='svr'><tau units='femtosecond'>100</tau></thermostat>"
PILE_L = "<thermostat mode='pile_l'><tau units='femtosecond'>100</tau></thermostat>"
BARO_ISOTROPIC = (
    "<barostat mode='isotropic'><tau units='femtosecond'>1000</tau>"
    f"{LANGEVIN}</barostat>"
)
BARO_FLEXIBLE = (
    "<barostat mode='flexible'><tau units='femtosecond'>1000</tau>"
    f"{LANGEVIN}</barostat>"
)

CASES = {
    "nve": dict(nbeads=1, mode="nve"),
    "nvt_langevin": dict(nbeads=1, mode="nvt", thermostat=LANGEVIN),
    "nvt_svr": dict(nbeads=1, mode="nvt", thermostat=SVR),
    "npt_isotropic": dict(
        nbeads=1, mode="npt", thermostat=LANGEVIN, barostat=BARO_ISOTROPIC
    ),
    "npt_flexible": dict(
        nbeads=1, mode="npt", thermostat=LANGEVIN, barostat=BARO_FLEXIBLE
    ),
    "pimd4_nvt": dict(nbeads=4, mode="nvt", thermostat=PILE_L),
    "pimd_nvt": dict(nbeads=8, mode="nvt", thermostat=PILE_L),
    "pimd_npt": dict(nbeads=8, mode="npt", thermostat=PILE_L, barostat=BARO_ISOTROPIC),
    "pimd32_npt": dict(
        nbeads=32, mode="npt", thermostat=PILE_L, barostat=BARO_ISOTROPIC
    ),
}

NATOMS = [100, 316, 1000, 3162, 10000]
SUPERCELLS = ["1", "2x1x1", "2x2x1", "2", "3", "4", "5"]

# forcefield definition and force component for each forcefield mode, with
# placeholders for the model, its device and the batch size
FORCEFIELDS = {
    "dummy": (
        "<ffdirect name='ff'><pes>dummy</pes></ffdirect>",
        "<force forcefield='ff'/>",
    ),
    "none": (
        "<fflj name='ff' pbc='False'><parameters>{{eps: 0.0, sigma: 1.0}}"
        "</parameters></fflj>",
        "<force forcefield='ff' weight='0.0'/>",
    ),
    "metatomic": (
        "<ffdirect name='ff' pbc='false'><pes>metatomic</pes>"
        "<parameters>{{template:init.xyz,model:{model},device:{model_device}}}"
        "</parameters><batch_size>{nbeads}</batch_size></ffdirect>",
        "<force forcefield='ff'/>",
    ),
}

# periodic structures that are replicated into supercells
EXAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../examples")
STRUCTURES = {
    "water": os.path.join(EXAMPLES, "init_files/water_32.pdb"),
    "nickel": os.path.join(EXAMPLES, "clients/metatomic/lj/nickel.xyz"),
}


def parse_supercell(spec):
    """Converts `n` or `axbxc` into a tuple of three integers."""

    reps = [int(n) for n in spec.split("x")]
    return tuple(reps * 3) if len(reps) == 1 else tuple(reps)


def write_input(args, seed=12345):
    """Writes `init.xyz` in the current folder, and returns the XML input and
    the number of atoms."""

    if args.structure == "random":
        # hydrogen atoms at random positions, with a density of 0.05 atoms/A^3
        natoms = int(args.size)
        box = (natoms / 0.05) ** (1.0 / 3.0)
        pos = np.random.default_rng(seed).random((natoms, 3)) * box
        with open("init.xyz", "w") as f:
            f.write(f"{natoms}\n")
            f.write(
                f"# CELL(abcABC): {box} {box} {box} 90 90 90 "
                "cell{angstrom} positions{angstrom}\n"
            )
            for x, y, z in pos:
                f.write(f"H {x:.8f} {y:.8f} {z:.8f}\n")
        initfile = "<file mode='xyz'>init.xyz</file>"
    else:
        import ase.io

        atoms = ase.io.read(STRUCTURES[args.structure])
        atoms.pbc = True
        atoms = atoms.repeat(parse_supercell(args.size))
        ase.io.write("init.xyz", atoms)
        natoms = len(atoms)
        initfile = "<file mode='ase' units='angstrom'>init.xyz</file>"

    setup = CASES[args.case]
    pressure = "<pressure units='bar'>1</pressure>" if "barostat" in setup else ""
    ffxml, forcexml = FORCEFIELDS[args.forcefield]
    ffxml = ffxml.format(
        model=args.model, model_device=args.model_device, nbeads=setup["nbeads"]
    )
    xml = f"""
<simulation verbosity='quiet' threading='false' safe_stride='100000000'>
  {ffxml}
  <total_steps>100000000</total_steps>
  <output prefix='benchmark'/>
  <prng><seed>{seed}</seed></prng>
  <system>
    <initialize nbeads='{setup["nbeads"]}'>
      {initfile}
      <velocities mode='thermal' units='kelvin'>300</velocities>
    </initialize>
    <forces>{forcexml}</forces>
    <ensemble>
      <temperature units='kelvin'>300</temperature>{pressure}
    </ensemble>
    <motion mode='dynamics'>
      <dynamics mode='{setup["mode"]}'>
        <timestep units='femtosecond'>0.5</timestep>
        {setup.get("thermostat", "")}
        {setup.get("barostat", "")}
      </dynamics>
    </motion>
  </system>
</simulation>
"""
    return xml, natoms


def time_driver(sim, nsteps, repeats, synchronize):
    """Returns the time (ms) taken by the driver of the direct forcefield to
    evaluate the structures of one step, or None for other forcefields."""

    ff = sim.fflist["ff"]
    if not hasattr(ff, "driver"):
        return None
    try:
        from ipi.utils.array_backend import to_numpy
    except ImportError:
        to_numpy = np.asarray

    # the same structures that the forcefield receives, as a batch if it batches,
    # as arrays of the backend if the driver works on them
    h, q = sim.syslist[0].cell.h, sim.syslist[0].beads.q
    if getattr(ff, "backend_arrays", False):
        h = h.clone() if hasattr(h, "clone") else h.copy()
        pos = [
            (p.clone() if hasattr(p, "clone") else p.copy()).reshape(-1, 3) for p in q
        ]
    else:
        h = np.array(to_numpy(h))
        pos = [np.array(to_numpy(p)).reshape(-1, 3) for p in q]
    if ff.batch_size > 1:
        calls = [([h] * len(pos), pos)]
    else:
        calls = [(h, p) for p in pos]

    timings = []
    for _ in range(repeats):
        synchronize()
        t = time.perf_counter()
        for _ in range(nsteps):
            for cell, positions in calls:
                ff.driver(cell, positions)
        synchronize()
        timings.append((time.perf_counter() - t) / nsteps * 1e3)
    return timings


def run_worker(args):
    """Times one (case, size) combination and writes a JSON result."""

    import ipi

    synchronize = lambda: None
    if args.backend != "numpy" or args.device != "cpu":
        from ipi.utils.array_backend import set_backend, synchronize

        set_backend(args.backend, device=args.device, jit=args.jit)
    if args.threads is not None and args.backend == "torch":
        import torch

        torch.set_num_threads(args.threads)

    from ipi.scripting import InteractiveSimulation

    os.chdir(tempfile.mkdtemp(prefix="ipi-benchmark-"))
    xml, natoms = write_input(args)

    tinit = time.perf_counter()
    sim = InteractiveSimulation(xml)
    tinit = time.perf_counter() - tinit

    step = 0

    def advance(nsteps):
        """Runs `nsteps` MD steps and returns the elapsed wall-clock time."""
        nonlocal step
        synchronize()
        t = time.perf_counter()
        for _ in range(nsteps):
            sim.run_step(step)
            step += 1
        synchronize()
        return time.perf_counter() - t

    advance(args.warmup)

    # estimates the time per step, doubling the number of steps until the
    # probe lasts at least a tenth of the block time
    nprobe = 1
    while True:
        elapsed = advance(nprobe)
        if elapsed > 0.1 * args.min_time:
            break
        nprobe *= 2
    nblock = max(1, int(np.ceil(args.min_time * nprobe / elapsed)))

    timings = [advance(nblock) / nblock * 1e3 for _ in range(args.repeats)]
    driver_timings = time_driver(sim, nblock, args.repeats, synchronize)

    result = dict(
        case=args.case,
        size=args.size,
        natoms=natoms,
        nbeads=CASES[args.case]["nbeads"],
        steps_per_block=nblock,
        ms_per_step=timings,
        median=statistics.median(timings),
        min=min(timings),
        driver_ms_per_step=driver_timings,
        driver_median=driver_timings and statistics.median(driver_timings),
        init_time=tinit,
        ipi_path=os.path.dirname(os.path.dirname(os.path.abspath(ipi.__file__))),
    )
    with open(args.result, "w") as f:
        json.dump(result, f)


def source_version(root):
    """Returns the git commit of an i-PI source tree, if available."""

    try:
        commit = subprocess.run(
            ["git", "-C", root, "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", root, "status", "--porcelain", "-uno"],
            capture_output=True,
            text=True,
        ).stdout.strip()
        return commit + ("-dirty" if dirty else "")
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def metadata(args):
    """Collects information on the machine and the software versions."""

    meta = dict(
        label=args.label,
        forcefield=args.forcefield,
        structure=args.structure,
        model=args.model,
        model_device=args.model_device,
        backend=args.backend,
        device=args.device,
        jit=args.jit,
        threads=args.threads,
        hostname=platform.node(),
        python=platform.python_version(),
        numpy=np.__version__,
        date=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    try:
        with open("/proc/cpuinfo") as f:
            meta["cpu"] = next(
                line.split(":", 1)[1].strip() for line in f if "model name" in line
            )
    except (OSError, StopIteration):
        meta["cpu"] = platform.processor()
    for package in ["torch", "metatomic-torch", "upet"]:
        try:
            meta[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            pass
    if args.backend == "torch" or args.model_device.startswith("cuda"):
        import torch

        if torch.cuda.is_available():
            meta["gpu"] = torch.cuda.get_device_name()
    return meta


def run_all(args):
    """Runs all the requested combinations, each in its own subprocess."""

    env = dict(os.environ)
    if args.ipi_root is None:
        args.ipi_root = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../.."
        )
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.abspath(args.ipi_root), env.get("PYTHONPATH", "")]
    )
    if args.threads is not None:
        for var in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
            env[var] = str(args.threads)

    output = dict(metadata=metadata(args), results=[])
    for case in args.cases:
        for size in args.sizes:
            with tempfile.NamedTemporaryFile(suffix=".json") as tmp:
                cmd = [sys.executable, os.path.abspath(__file__), "--worker"]
                cmd += ["--case", case, "--size", size, "--result", tmp.name]
                cmd += ["--backend", args.backend, "--device", args.device]
                cmd += ["--forcefield", args.forcefield, "--structure", args.structure]
                cmd += ["--model", str(args.model), "--model-device", args.model_device]
                cmd += ["--warmup", str(args.warmup), "--repeats", str(args.repeats)]
                cmd += ["--min-time", str(args.min_time)]
                if args.jit:
                    cmd += ["--jit"]
                if args.threads is not None:
                    cmd += ["--threads", str(args.threads)]
                proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
                if proc.returncode != 0 or os.path.getsize(tmp.name) == 0:
                    print(f"{case:>14s} {size:>6s}  FAILED", file=sys.stderr)
                    print(proc.stdout[-2000:], proc.stderr[-2000:], file=sys.stderr)
                    continue
                with open(tmp.name) as f:
                    result = json.load(f)
            output["results"].append(result)
            driver = result.get("driver_median")
            print(
                f"{case:>14s} {result['natoms']:6d}  {result['median']:10.4f} ms/step"
                f"  (min {result['min']:.4f}, init {result['init_time']:.1f} s"
                + (f", driver {driver:.4f} ms/step)" if driver else ")"),
                flush=True,
            )

    if output["results"]:
        root = output["results"][0]["ipi_path"]
        output["metadata"]["ipi_path"] = root
        output["metadata"]["ipi_version"] = source_version(root)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=1)


def load_timings(filenames, label, stat):
    """Reads comma-separated result files and returns the timings of each
    (case, natoms), taking the minimum or the median over the files."""

    timings = {}
    for filename in filenames.split(","):
        with open(filename) as f:
            data = json.load(f)
        m = data["metadata"]
        print(
            f"# {label:9s}: {m.get('label') or ''} {m.get('ipi_version', '')} "
            f"{m['backend']}/{m['device']} forcefield={m.get('forcefield', 'dummy')} "
            f"({m['date']})"
        )
        for r in data["results"]:
            timings.setdefault((r["case"], r["natoms"]), []).append(r[stat])
    aggregate = min if stat == "min" else statistics.median
    return {key: aggregate(values) for key, values in timings.items()}


def compare(reference, new, threshold, stat):
    """Prints the ratio of the timings of two (sets of) benchmark runs."""

    ref = load_timings(reference, "reference", stat)
    cur = load_timings(new, "new", stat)

    print(f"{'case':>14s} {'natoms':>6s} {'ref':>10s} {'new':>10s} {'new/ref':>8s}")
    nslow = 0
    for key, value in cur.items():
        if key not in ref:
            continue
        ratio = value / ref[key]
        flag = " *" if ratio > 1 + threshold else ""
        nslow += bool(flag)
        print(
            f"{key[0]:>14s} {key[1]:6d} {ref[key]:10.4f} {value:10.4f}"
            f" {ratio:8.3f}{flag}"
        )
    print(f"# {nslow} combinations slower than reference by more than {threshold:.0%}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--cases",
        default=",".join(CASES),
        help="comma-separated list of cases, among: " + ", ".join(CASES),
    )
    parser.add_argument(
        "-n",
        "--natoms",
        default=",".join(str(n) for n in NATOMS),
        help="comma-separated list of system sizes, for the random structure",
    )
    parser.add_argument(
        "-s",
        "--supercells",
        default=",".join(SUPERCELLS),
        help="comma-separated list of supercells (n or axbxc) of the periodic "
        "structures",
    )
    parser.add_argument("-o", "--output", default="benchmark.json")
    parser.add_argument("--label", default="", help="free-text label for the run")
    parser.add_argument(
        "--forcefield",
        choices=list(FORCEFIELDS),
        default="dummy",
        help="'dummy' evaluates zero forces, 'none' skips force evaluation, "
        "'metatomic' evaluates --model",
    )
    parser.add_argument(
        "--structure",
        choices=["random"] + list(STRUCTURES),
        default="random",
        help="hydrogen atoms at random positions, or a periodic structure",
    )
    parser.add_argument("--model", default=None, help="metatomic model file")
    parser.add_argument(
        "--model-device", default="cpu", help="device of the metatomic model"
    )
    parser.add_argument("--backend", default="numpy")
    parser.add_argument(
        "--jit",
        action="store_true",
        help="compiles the steps (CUDA graphs with torch on CUDA, jax.jit with jax)",
    )
    parser.add_argument("--device", default="cpu")
    parser.add_argument(
        "--threads", type=int, default=None, help="number of threads for numpy/torch"
    )
    parser.add_argument("--ipi-root", default=None, help="i-PI source tree to test")
    parser.add_argument("--warmup", type=int, default=5, help="untimed steps")
    parser.add_argument("--repeats", type=int, default=5, help="number of timed blocks")
    parser.add_argument(
        "--min-time", type=float, default=1.0, help="minimum duration of a block (s)"
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("REFERENCE", "NEW"),
        help="compare two result files (or comma-separated lists of files from "
        "repeated runs) instead of running",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.03,
        help="relative slowdown flagged by --compare",
    )
    parser.add_argument(
        "--stat",
        choices=["min", "median"],
        default="min",
        help="timing statistic used by --compare",
    )
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--case", help=argparse.SUPPRESS)
    parser.add_argument("--size", help=argparse.SUPPRESS)
    parser.add_argument("--result", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.compare:
        compare(*args.compare, args.threshold, args.stat)
    elif args.worker:
        run_worker(args)
    else:
        args.cases = args.cases.split(",")
        unknown = set(args.cases) - set(CASES)
        if unknown:
            parser.error(f"unknown cases: {', '.join(sorted(unknown))}")
        if args.forcefield == "metatomic":
            if args.model is None or args.structure == "random":
                parser.error("metatomic needs --model and a periodic --structure")
            # the workers run in temporary folders
            args.model = os.path.abspath(args.model)
        sizes = args.natoms if args.structure == "random" else args.supercells
        args.sizes = sizes.split(",")
        run_all(args)


if __name__ == "__main__":
    main()

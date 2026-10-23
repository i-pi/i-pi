"""Run independent mathematical checks and write reports outside the checkout.

Example:
  python ipi_tests/run_math_verification.py --baseline /outside/inventory.json \
      --report-dir /outside/verification --mutations

The baseline is the external mathematical-check inventory. Test functions
declare the entries they exercise with @verifies. Reports keep reference
presence, actual execution and expected failures separate. Mutation checks
modify only an isolated temporary copy of the Python package and unit tests.
"""

import argparse
import ast
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
TEST_DIR = Path("ipi_tests/unit_tests/engine")
TEST_FILES = [
    TEST_DIR / ("test_" + name + "_mathematical_references.py")
    for name in ["property", "force", "initializer", "virial"]
]
GROUPS = ["properties", "trajectories", "properties_helpers", "forces", "initializer"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_index():
    index = {}
    tests = {}
    for relative in TEST_FILES:
        for node in ast.walk(ast.parse((ROOT / relative).read_text())):
            if not isinstance(node, ast.FunctionDef):
                continue
            for decorator in node.decorator_list:
                if (
                    isinstance(decorator, ast.Call)
                    and isinstance(decorator.func, ast.Name)
                    and decorator.func.id == "verifies"
                ):
                    entries = [ast.literal_eval(arg) for arg in decorator.args]
                    tests[node.name] = dict(
                        path=str(ROOT / relative), line=node.lineno, entries=entries
                    )
                    for entry in entries:
                        index.setdefault(entry, []).append(node.name)
    return index, tests


def independence_check():
    reference = ROOT / TEST_DIR / "mathematical_references.py"
    tree = ast.parse(reference.read_text())
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports += [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
        elif (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in ["eval", "exec", "__import__"]
        ):
            raise ValueError(
                "Reference calculations must not dynamically load production code"
            )
    allowed = {"itertools", "math", "numpy", "numpy.polynomial"}
    if set(imports) - allowed:
        raise ValueError(
            "Unexpected reference dependency: " + repr(set(imports) - allowed)
        )
    return dict(
        path=str(reference), sha256=sha(reference), imports=imports, status="passed"
    )


def case_checks_entry(case, entry, tests):
    """Restrict parameterized output tests to the output actually selected."""
    name = case["function"]
    if name not in tests:
        return False
    identifiers = case["name"].partition("[")[2].removesuffix("]")
    if name == "test_sc_corrected_virial_differentiates_complete_energy":
        part = "virssc_part_2" if "virssc_part_2" in identifiers else "virssc"
        return entry == "forces:Forces.get_" + part
    mentioned = {
        candidate
        for candidate in tests[name]["entries"]
        if re.search(
            r"(^|-)" + re.escape(candidate.split(":", 1)[1]) + r"(-|$)",
            identifiers,
        )
    }
    if mentioned:
        return entry in mentioned
    if name == "test_empty_mts_weights_place_unpartitioned_component_on_outer_level":
        suffix = {"f": "get_forces_mts", "virs": "get_virs_mts", "vir": "get_vir_mts"}[
            identifiers
        ]
        return entry.endswith("." + suffix)
    return entry in tests[name]["entries"]


def run_pytest(root, targets, directory, label):
    xml = directory / (label + ".xml")
    log = directory / (label + ".log")
    command = [
        sys.executable,
        "-m",
        "pytest",
        *map(str, targets),
        "-q",
        "--tb=short",
        "-o",
        "cache_dir=" + str(directory / "pytest-cache"),
        "--junitxml=" + str(xml),
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(root)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    with log.open("w") as stream:
        result = subprocess.run(
            command,
            cwd=root,
            env=environment,
            stdout=stream,
            stderr=subprocess.STDOUT,
            timeout=300,
        )
    cases = []
    if xml.exists():
        for case in ET.parse(xml).iter("testcase"):
            skipped, failure, error = (
                case.find("skipped"),
                case.find("failure"),
                case.find("error"),
            )
            status = (
                "error"
                if error is not None
                else (
                    "failed"
                    if failure is not None
                    else (
                        "expected failure"
                        if skipped is not None and skipped.get("type") == "pytest.xfail"
                        else "skipped" if skipped is not None else "passed"
                    )
                )
            )
            detail = (
                error
                if error is not None
                else failure if failure is not None else skipped
            )
            cases.append(
                dict(
                    name=case.get("name"),
                    function=case.get("name").split("[")[0],
                    status=status,
                    message="" if detail is None else detail.get("message", ""),
                    assertion_failure=failure is not None
                    and (
                        "AssertionError"
                        in (failure.get("message", "") + (failure.text or ""))
                        or "assert "
                        in (failure.get("message", "") + (failure.text or ""))
                    ),
                )
            )
    return dict(
        returncode=result.returncode,
        command=command,
        log=str(log),
        xml=str(xml),
        counts=dict(Counter(c["status"] for c in cases)),
        cases=cases,
    )


MUTATIONS = [
    (
        "density-division",
        "properties",
        None,
        "self.beads.m.sum() / self.cell.V",
        "self.beads.m.sum() * self.cell.V",
        "property",
        "test_cell_and_density_outputs_match_lattice_geometry",
    ),
    (
        "spring-bead-factor",
        "properties",
        "__init__",
        "self.nm.vspring / self.beads.nbeads",
        "self.nm.vspring * self.beads.nbeads",
        "property",
        "test_spring_primitive_conserved_and_bias_outputs_match_explicit_actions",
    ),
    (
        "mode-kinetic-half",
        "properties",
        "get_kinmd",
        "2.0 * dm3[inm, k]",
        "4.0 * dm3[inm, k]",
        "property",
        "test_mass_scaled_kinetic_energy_and_stress_use_prescribed_mode_masses",
    ),
    (
        "sc-linear-sign",
        "forces",
        "get_coeffsc_part_1",
        "rc[0::2] = -1.0 / 3.0",
        "rc[0::2] = 1.0 / 3.0",
        "force",
        "test_individual_sc_coefficients_and_corrections_follow_polynomial_derivatives",
    ),
    (
        "sc-force-coefficient",
        "forces",
        "get_fsc_part_1",
        "self.coeffsc_part_1 *",
        "self.coeffsc_part_2 *",
        "force",
        "test_individual_sc_coefficients_and_corrections_follow_polynomial_derivatives",
    ),
    (
        "mts-level-weight",
        "forces",
        "get_forces_mts",
        "* mts_weights[level]",
        "* 1.0",
        "force",
        "test_noncentroid_components_and_mts_levels_match_explicit_fourier_projection",
    ),
    (
        "thermal-mass-square-root",
        "initializer",
        "init_stage1",
        "np.sqrt(rnm.dynm3)",
        "rnm.dynm3",
        "initializer",
        "test_thermal_initialization_with_nonphysical_mode_masses_uses_fixed_draws",
    ),
    (
        "momentum-bead-square-root",
        "initializer",
        "init_stage1",
        "rp *= np.sqrt(self.nbeads / nbeads)",
        "rp *= self.nbeads / nbeads",
        "initializer",
        "test_nonconstant_path_initialization_matches_independent_fourier_projection",
    ),
    (
        "velocity-mass-multiply",
        "initializer",
        "init_stage1",
        "ev *= simul.beads.m3[0]",
        "ev /= simul.beads.m3[0]",
        "initializer",
        "test_initial_mass_cell_and_velocity_units_use_independent_si_anchors",
    ),
    (
        "fourth-order-isotope-ti",
        "properties",
        "get_isotope_zetatd_4th",
        "(1.0 / 24.0)",
        "(1.0 / 12.0)",
        "property",
        "test_fourth_order_isotope_weights_include_nonzero_force_and_potential_changes",
    ),
]


def mutate_source(source, function, old, new):
    if function is None:
        start, stop = 0, len(source)
    else:
        node = next(
            n
            for n in ast.walk(ast.parse(source))
            if isinstance(n, ast.FunctionDef) and n.name == function
        )
        lines = source.splitlines(keepends=True)
        start = sum(map(len, lines[: node.lineno - 1]))
        stop = sum(map(len, lines[: node.end_lineno]))
    part = source[start:stop]
    if part.count(old) != 1:
        raise ValueError("Mutation target is absent or ambiguous: " + old)
    result = source[:start] + part.replace(old, new, 1) + source[stop:]
    ast.parse(result)
    return result


def run_mutations(directory):
    results = []
    with tempfile.TemporaryDirectory(prefix="ipi-math-mutants-") as temporary:
        copy = Path(temporary)
        ignore = shutil.ignore_patterns("__pycache__", ".pytest_cache")
        shutil.copytree(ROOT / "ipi", copy / "ipi", ignore=ignore)
        (copy / "ipi_tests").mkdir()
        shutil.copy2(ROOT / "ipi_tests/__init__.py", copy / "ipi_tests/__init__.py")
        shutil.copytree(
            ROOT / "ipi_tests/unit_tests", copy / "ipi_tests/unit_tests", ignore=ignore
        )
        for label, module, function, old, new, test_module, test_name in MUTATIONS:
            target = copy / "ipi/engine" / (module + ".py")
            original = target.read_text()
            try:
                target.write_text(mutate_source(original, function, old, new))
                node = (
                    str(
                        TEST_DIR
                        / ("test_" + test_module + "_mathematical_references.py")
                    )
                    + "::"
                    + test_name
                )
                run = run_pytest(copy, [node], directory, "mutation-" + label)
                killed = (
                    run["returncode"] == 1
                    and run["counts"].get("failed", 0) > 0
                    and not run["counts"].get("error", 0)
                    and all(
                        case["assertion_failure"]
                        for case in run["cases"]
                        if case["status"] == "failed"
                    )
                    and all("XPASS" not in case["message"] for case in run["cases"])
                )
                results.append(
                    dict(
                        label=label,
                        module=module,
                        function=function,
                        old=old,
                        new=new,
                        detected=killed,
                        execution=run,
                    )
                )
                print(
                    "Mutation "
                    + label
                    + ": "
                    + ("detected" if killed else "NOT DETECTED"),
                    flush=True,
                )
            finally:
                target.write_text(original)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--mutations", action="store_true")
    args = parser.parse_args()
    directory = args.report_dir.resolve()
    if directory.is_relative_to(ROOT):
        parser.error("The report directory must be outside the repository")
    directory.mkdir(parents=True, exist_ok=True)
    baseline = json.loads(args.baseline.read_text())
    pending = baseline.get("deferred_definitions", {})
    marked = {
        group + ":" + key: row
        for group in GROUPS
        for key, row in baseline[group].items()
        if row["status"].startswith("MARK:")
    }
    index, tests = reference_index()
    missing = sorted(set(marked) - set(index) - set(pending))
    if missing:
        raise ValueError("Unmapped marked entries: " + repr(missing))
    independence = independence_check()
    production = [
        ROOT / "ipi/engine" / (name + ".py")
        for name in ["properties", "forces", "initializer", "barostats"]
    ]
    before = {str(path): sha(path) for path in production}
    execution = run_pytest(ROOT, TEST_FILES, directory, "references")
    print("References: " + repr(execution["counts"]), flush=True)
    mutations = (
        run_mutations(directory)
        if args.mutations and execution["returncode"] == 0
        else []
    )
    assert before == {
        str(path): sha(path) for path in production
    }, "Production source changed during verification"
    rows = {}
    for entry, previous in marked.items():
        names = index.get(entry, [])
        cases = [
            case for case in execution["cases"] if case_checks_entry(case, entry, tests)
        ]
        counts = dict(Counter(case["status"] for case in cases))
        status = (
            "DEFERRED SCIENTIFIC DEFINITION"
            if entry in pending
            else (
                "FAILED CHECK"
                if counts.get("failed") or counts.get("error")
                else (
                    "REFERENCE EXISTS; LINKED CASES INCLUDE KNOWN FAILURES"
                    if counts.get("expected failure")
                    else (
                        "INDEPENDENT REFERENCE CASES PASS"
                        if counts.get("passed")
                        else "NOT EXECUTED"
                    )
                )
            )
        )
        rows[entry] = dict(
            status=status,
            previous=previous,
            tests=[dict(name=name, **tests[name]) for name in names],
            linked_case_counts=counts,
            pending=pending.get(entry),
            scope="These are finite reference cases, not a proof for all inputs. Shared tests may exercise several entries; counts are not additive across rows.",
        )
    payload = dict(
        baseline=str(args.baseline.resolve()),
        production_sha256=before,
        independence=independence,
        execution=execution,
        mutations=mutations,
        entries=rows,
    )
    (directory / "verification.json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = [
        "# Automated independent mathematical verification",
        "",
        "Reports are outside the repository. Expected failures remain unfixed bugs; no production implementation was changed.",
        "",
        "## Results",
        "",
        "- Marked entries: " + str(len(marked)),
        "- Entries with reference tests: " + str(len(set(marked) & set(index))),
        "- Deferred definitions: " + str(len(set(marked) & set(pending))),
        "- Executed cases: " + repr(execution["counts"]),
        "- Temporary mutations detected: "
        + str(sum(m["detected"] for m in mutations))
        + "/"
        + str(len(mutations)),
        "",
        "Reference dependency guard: "
        + independence["status"]
        + ". The independent numerical helper imports only NumPy and the standard library. This structural guard does not prove scientific correctness; equations and cases still require review.",
        "",
        "The reference cases use nonzero convex quartic forces, unequal masses, explicit Fourier projection, nonunit component/MTS weights, fixed Gaussian draws, independent SI anchors, and explicit two-particle permutation sums. SC coefficients come from arXiv:1606.00920 Eqs.7-8 and Appendix C; beta derivatives from physics/0505109 Eqs.18,22,33,34; SC operator double virial from JCTC 15,3237 (2019), Eqs.10,14,15.",
        "Ordinary and SC force virials additionally use minus the cell-strain derivative of the complete energy at fixed fractional coordinates. The volume-coupled convex quartic model changes both energy and forces with cell volume. The SC expression is cross-checked against JCTC 15,3237 (2019), Eq.9 / SI S1.26--S1.31, accounting for the opposite paper sign. Separate checks cover force components, weighting, pressure outputs, MTS stresses with and without bias, BZP/SCBZP/MTK/RGB momentum contributions and the classical SCR pressure force. This does not establish full SC-NPT conserved-energy behavior or all contracted-path combinations.",
        "",
        "## Every marked entry",
        "",
        "Rows state the results of linked test groups. A group may cover several functions; an expected failure in that group is not automatically a distinct bug in every linked function. Read the failure list below. Passing reference cases establish the listed examples, not all possible parameters or feature combinations.",
        "",
        "| Entry | Status | Linked cases | Previous gap | Tests / source |",
        "| --- | --- | --- | --- | --- |",
    ]
    for entry, row in rows.items():
        links = (
            ", ".join(
                "[" + t["name"] + "](" + t["path"] + ":" + str(t["line"]) + ")"
                for t in row["tests"]
            )
            or row["pending"]
        )
        lines.append(
            "| `"
            + entry
            + "` | "
            + row["status"]
            + " | "
            + str(row["linked_case_counts"])
            + " | "
            + row["previous"]["gap"].replace("|", "\\|")
            + " | "
            + links
            + " |"
        )
    lines += ["", "## Known failures observed", ""]
    failures = sorted(
        {
            case["message"]
            for case in execution["cases"]
            if case["status"] == "expected failure"
        }
    )
    lines += ["- " + reason for reason in failures]
    lines += ["", "## Deferred definitions", ""] + [
        "- `" + entry + "`: " + reason for entry, reason in pending.items()
    ]
    lines += [
        "",
        "## Mutation sensitivity",
        "",
        "Each mutation was made only in an isolated package copy and required an assertion failure, not a setup/import error or an unexpected pass. This demonstrates sensitivity to those mistakes only.",
        "",
    ]
    lines += [
        "- "
        + m["label"]
        + ": "
        + ("detected" if m["detected"] else "NOT DETECTED")
        + "; `"
        + m["old"]
        + "` -> `"
        + m["new"]
        + "`."
        for m in mutations
    ]
    lines += [
        "",
        "## Reproduce",
        "",
        "```sh",
        shlex.join(
            [
                sys.executable,
                "ipi_tests/run_math_verification.py",
                "--baseline",
                str(args.baseline.resolve()),
                "--report-dir",
                str(directory),
            ]
            + (["--mutations"] if args.mutations else [])
        ),
        "```",
        "",
        "Machine-readable mapping, commands, source hashes and results: [verification.json]("
        + str(directory / "verification.json")
        + ").",
        "",
    ]
    (directory / "verification.md").write_text("\n".join(lines))
    return (
        1
        if execution["returncode"] != 0 or any(not m["detected"] for m in mutations)
        else 0
    )


if __name__ == "__main__":
    raise SystemExit(main())

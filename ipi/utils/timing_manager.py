"""Lightweight hierarchical timing utilities.

The timing manager is intentionally opt-in. Timers are active only when the
environment variable ``ENABLE_TIMING_MANAGER=1`` is set.

Basic usage:

.. code-block:: python

    from timing_manager import timers

    timers.start("Integrate")
    # ... code to benchmark ...
    timers.stop("Integrate")

    # Write ``run_timing.log`` and ``run_timing.dat`` and print the table.
    timers.summary("run_timing")

Hierarchy convention:

- A timer named ``Parent(+)`` declares ``+`` as the marker for its children.
- A timer named ``[+]Child`` is attached under the most recent parent that
  declares marker ``+``.

Supported markers are listed in ``TimingManager.group_symbols``.
"""

import time
from collections import defaultdict
import numpy as np
import os
import pandas as pd
from io import StringIO
import sys

ENABLED = os.getenv("ENABLE_TIMING_MANAGER", "0") == "1"
try:
    DEFAULT_SKIP_FIRST = int(os.getenv("TIMING_MANAGER_SKIP_FIRST", "5"))
except ValueError:
    DEFAULT_SKIP_FIRST = 5

class TimingManager:
    """Collects elapsed times and writes a human-readable + tabular summary.

    Notes:
    - Elapsed values are stored in milliseconds.
    - Individual measurements are rounded and stored as ``np.float32``.
    - When disabled, ``start`` and ``stop`` become no-op callables to keep
      instrumentation overhead minimal.
    """

    def __init__(self):
        if ENABLED:
            self._start_times = {}
            self._all_times = defaultdict(
                lambda: {
                    "sum": 0.0,
                    "count": 0,
                    "max": 0.0,
                    "sum_after": 0.0,
                    "count_after": 0,
                    "max_after": 0.0,
                }
            )
            self.group_symbols = [
                                    "+", "++", "+++", "++++","+++++","++++++",
                                    "*", "**", "***", "****",
                                    "#", "##", "###", "####"
                                  ]  # supported hierarchy roots
            self.decimals = 4
            self.skip_first_default = max(DEFAULT_SKIP_FIRST, 0)
            self.current_step = 0
            self.steps_seen = 0
        else:
            # empty functions in isabled mode to avoid overhead
            self.start = lambda *a, **k: None
            self.stop = lambda *a, **k: None
            self.set_step = lambda *a, **k: None

    def start(self, name):
        """Start or restart a named timer."""
        self._start_times[name] = time.perf_counter()

    def stop(self, name):
        """Stop a started timer and record elapsed milliseconds.

        Args:
            name: Timer name passed to :meth:`start`.

        Returns:
            Elapsed time in milliseconds as ``np.float32``.

        Raises:
            ValueError: If ``start(name)`` was not called first.
        """
        if name not in self._start_times:
            raise ValueError(f"Timer '{name}' was not started.")
        elapsed_sec = time.perf_counter() - self._start_times.pop(name)
        elapsed_ms = np.float32(round(elapsed_sec * 1000, self.decimals))
        stats = self._all_times[name]
        elapsed = float(elapsed_ms)
        stats["sum"] += elapsed
        stats["count"] += 1
        if stats["count"] == 1 or elapsed > stats["max"]:
            stats["max"] = elapsed

        self.steps_seen = max(self.steps_seen, self.current_step + 1)
        if self.current_step >= self.skip_first_default:
            stats["sum_after"] += elapsed
            stats["count_after"] += 1
            if stats["count_after"] == 1 or elapsed > stats["max_after"]:
                stats["max_after"] = elapsed
        return elapsed_ms

    def set_step(self, step):
        """Set the current simulation step for step-normalized statistics."""
        self.current_step = max(int(step), 0)
        self.steps_seen = max(self.steps_seen, self.current_step + 1)

    def get_total(self, name):
        """Return total measured time (ms) for ``name`` across all samples."""
        return round(self._all_times[name]["sum"], self.decimals)

    def get_count(self, name):
        """Return number of recorded samples for ``name``."""
        return self._all_times[name]["count"]

    def get_average(self, name, skip_first=0):
        """Return average time (ms) for ``name``, optionally skipping warmup.

        Args:
            name: Timer name.
            skip_first: If ``> 0``, use the configured "after N steps" average.
        """
        if skip_first <= 0:
            steps = max(self.steps_seen, 1)
            return round(self._all_times[name]["sum"] / steps, self.decimals)
        steps_after = max(self.steps_seen - skip_first, 0)
        if steps_after == 0:
            return 0.0
        return round(self._all_times[name]["sum_after"] / steps_after, self.decimals)

    def get_max(self, name, skip_first=0):
        """Return maximum recorded value (ms) for ``name``.

        Args:
            name: Timer name.
            skip_first: If ``> 0``, use the configured "after N steps" max.
        """
        stats = self._all_times[name]
        if skip_first <= 0:
            return stats["max"]
        return stats["max_after"]

    def reset(self, name=None):
        """Clear measurements for one timer or all timers.

        Args:
            name: If provided, clear only this timer. If ``None``, clear all
                running and recorded timers.
        """
        if name is None:
            self._start_times.clear()
            self._all_times.clear()
            self.current_step = 0
            self.steps_seen = 0
        else:
            self._start_times.pop(name, None)
            self._all_times.pop(name, None)

    def _parse_marker(self, name):
        """Parse hierarchy markers embedded in a timer name.

        Returns:
            A tuple ``(marker_of_parent, marker_of_children)`` where:
            - ``marker_of_parent`` is extracted from ``[marker]Name``.
            - ``marker_of_children`` is extracted from ``Name(marker)``.
            Missing markers are returned as ``None``.
        """
        marker_of_parent = None
        marker_of_children = None

        if "(" in name and name.endswith(")"):
            _, marker = name.rsplit("(", 1)
            marker_of_children =  marker[:-1]
        if "]" in name and name.startswith("["):
            marker, _ = name.rsplit("]", 1)
            marker_of_parent =  marker[1:]

        return marker_of_parent, marker_of_children

    def _display_name(self, name):
        """Strip hierarchy markers from timer name for user-facing output."""
        marker_of_parent, marker_of_children = self._parse_marker(name)
        display = name

        if marker_of_parent is not None and display.startswith("[") and "]" in display:
            display = display.split("]", 1)[1]

        if marker_of_children is not None and display.endswith(")") and "(" in display:
            display = display.rsplit("(", 1)[0]

        return display.strip()



    def summary(self, filename, skip_first=None, tree_numbers=False):
        """Generate and print a hierarchical timing report.

        This method prints a formatted table, writes ``{filename}.log`` with the
        same console output, and writes ``{filename}.dat`` as ``;``-separated
        tabular data.

        Args:
            filename: Output path prefix for ``.log`` and ``.dat`` files.
            skip_first: Number of simulation steps to ignore in ``*_aftern``
                columns. Defaults to ``TIMING_MANAGER_SKIP_FIRST`` or ``5``.
            tree_numbers: If ``True``, visually indent the numeric "after n"
                columns to align with the tree depth of each timer.

        Returns:
            ``None``. In disabled mode it exits immediately.
        """

        if not ENABLED:
            return None
        if skip_first is None:
            skip_first = self.skip_first_default
        skip_first = max(int(skip_first), 0)
        if skip_first != self.skip_first_default:
            raise ValueError(
                "In accumulation mode, skip_first is fixed from startup. "
                "Set TIMING_MANAGER_SKIP_FIRST before running."
            )

        # Capture console output and also write to log file
        log_buffer = StringIO()
        old_stdout = sys.stdout
        sys.stdout = log_buffer

        # We also prepare DataFrame rows
        rows = []

        # Determine global total time (sum of top-level timers)
        top_level = [n for n in self._all_times.keys()
                     if self._parse_marker(n)[0] is None]
        global_total = sum(self.get_total(n) for n in top_level)

        # Helper to register a row
        def record_row(name, avg_skip, avg_all, max_skip, max_all,
                       count, total, indent_level, parent_total):
            total_pct = (total / global_total * 100) if global_total > 0 else 0.0
            parent_pct = (total / parent_total * 100) if parent_total > 0 else 0.0

            rows.append({
                "Name": name,
                "Indent": indent_level,
                "Avg_aftern(ms/step)": avg_skip,
                "Total%": total_pct,
                "Parent%": parent_pct,
                "Avg(ms/step)": avg_all,
                "Max_aftern(ms/step)": max_skip,
                "Max(ms/step)": max_all,
                "Count": count,
                "Total(ms)": total
            })


        # ---------------- Print identical header ----------------
        print("\n")
        name_w = 40
        avg_w = 20
        max_w = 20
        count_w = 8
        total_w = 12
        table_w = name_w + avg_w + max_w + count_w + total_w + 4
        header_fmt = (f"{{:<{name_w}}} {{:>{avg_w}}} {{:>{max_w}}} "
                      f"{{:>{count_w}}} {{:>{total_w}}}")
        row_fmt = (f"{{:<{name_w}}} {{:>{avg_w}}} {{:>{max_w}}} "
                   f"{{:>{count_w}}} {{:>{total_w}}}")

        print(header_fmt.format("Timer", "Avg_aftern(ms/step)",
                                "Max_aftern(ms/step)", "Count", "Total(ms)"))
        print("-" * table_w)

        all_names = sorted(self._all_times.keys())
        visited = set()

        # Build tree connectors for the name column.
        def build_tree_prefix(ancestor_has_more, is_last):
            if is_last is None:
                return "", 0
            stem = "".join("│   " if has_more else "    "
                           for has_more in ancestor_has_more)
            branch = "└─ " if is_last else "├─ "
            indent_width = 4 * len(ancestor_has_more)
            return stem + branch, indent_width

        def _numeric_cell(value, width, prefix="", is_int=False):
            if is_int:
                value_str = f"{value:d}"
            else:
                value_str = f"{value:.3f}"

            if len(prefix) >= width:
                return prefix[:width]
            return f"{prefix}{value_str:>{width - len(prefix)}}"

        # Local helper for printing line AND recording DF row
        def print_line(name, skip_first, indent_level, parent_name=None,
                       ancestor_has_more=None, is_last=None):
            if ancestor_has_more is None:
                ancestor_has_more = []
            avg_all = self.get_average(name)
            avg_skip = self.get_average(name, skip_first)
            max_all = self.get_max(name)
            max_skip = self.get_max(name, skip_first)
            total = self.get_total(name)
            count = self.get_count(name)
            #for some reason indented number print doesnt work yet
            tree_prefix, indent_width = build_tree_prefix(ancestor_has_more, is_last)
            indent_width = len(tree_prefix)
            rendered_name = f"{tree_prefix}{self._display_name(name)}"
            indent_prefix = (" " * indent_width) if tree_numbers else ""
            avg_cell = _numeric_cell(avg_skip, avg_w, indent_prefix, is_int=False)
            max_cell = _numeric_cell(max_skip, max_w, "", is_int=False)
            count_cell = _numeric_cell(count, count_w, "", is_int=True)
            total_cell = _numeric_cell(total, total_w, "", is_int=False)

            print(row_fmt.format(rendered_name, avg_cell, max_cell,
                                 count_cell, total_cell))

            parent_total = self.get_total(parent_name) if parent_name else total
            record_row(name, avg_skip, avg_all, max_skip, max_all, count,
                       total, indent_level, parent_total)


        # ---------------- Recursive children printer ----------------
        def print_children(parent_name, marker, skip_first, indent_level,
                           ancestor_has_more):
            direct_children = [
                n for n in self._all_times.keys()
                if self._parse_marker(n)[0] == marker
            ]
            direct_children.sort()

            total_child_time = 0.0

            for i, child in enumerate(direct_children):
                visited.add(child)
                is_last_child = (i == len(direct_children) - 1)
                print_line(child, skip_first, indent_level, parent_name,
                           ancestor_has_more, is_last=is_last_child)

                # get_average already returns per-step cost, including repeated
                # calls within a step, so no extra count scaling is needed.
                total_child_time += self.get_average(child, skip_first)

                _, child_marker = self._parse_marker(child)
                if child_marker:
                    child_ancestor_has_more = ancestor_has_more + [not is_last_child]
                    print_children(child, child_marker, skip_first,
                                   indent_level + 1,
                                   ancestor_has_more=child_ancestor_has_more)

            # Remainder line (console only)
            if direct_children:
                parent_avg = self.get_average(parent_name, skip_first)
                remainder = round(parent_avg - total_child_time, 3)
                if parent_avg > 0:
                    percent = f"Remainder: {remainder/parent_avg*100:3.1f}%"
                else:
                    percent = 0
                    print(f"Error in timing manager: Parent time of [{marker}] is 0.")
                remainder_indent = "".join("│   " if has_more else "    "
                                           for has_more in ancestor_has_more)
                remainder_indent_width = 4 * len(ancestor_has_more)
                remainder_num_prefix = (" " * remainder_indent_width) if tree_numbers else ""
                print(row_fmt.format(remainder_indent + str(percent),
                                     _numeric_cell(remainder,
                                                   avg_w,
                                                   remainder_num_prefix,
                                                   is_int=False),
                                     "", "", ""))


        # ---------------- Top-level timers ----------------
        for name in all_names:
            if name in visited:
                continue

            marker_of_parent, marker_of_children = self._parse_marker(name)

            if marker_of_parent is None and marker_of_children:
                visited.add(name)
                print_line(name, skip_first, indent_level=0, parent_name=None)
                print_children(name, marker_of_children, skip_first,
                               indent_level=1, ancestor_has_more=[])
            else:
                visited.add(name)
                print_line(name, skip_first, indent_level=0, parent_name=None)

        # ---------------- Restore stdout ----------------
        sys.stdout = old_stdout

        # Write log file
        with open(f"{filename}.log", "w") as f:
            f.write(log_buffer.getvalue())

        # Write dataframe
        df = pd.DataFrame(rows)
        df.to_csv(f"{filename}.dat", sep=";", index=False, float_format='%.4f')

        # Print to console what user would normally see
        print(log_buffer.getvalue(), end="")


# Global singleton
timers = TimingManager()

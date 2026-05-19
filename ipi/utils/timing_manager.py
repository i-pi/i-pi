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
"""

import time
from collections import defaultdict
import os
import pandas as pd
from io import StringIO
import sys
import threading

ENABLED = os.getenv("ENABLE_TIMING_MANAGER", "0") == "1"
try:
    DEFAULT_SKIP_FIRST = int(os.getenv("TIMING_MANAGER_SKIP_FIRST", "5"))
except ValueError:
    DEFAULT_SKIP_FIRST = 5


class _TimerStats:
    __slots__ = ("sum", "count", "max", "sum_after", "count_after", "max_after")

    def __init__(self):
        self.sum = 0.0
        self.count = 0
        self.max = 0.0
        self.sum_after = 0.0
        self.count_after = 0
        self.max_after = 0.0


class _TimerNodeStats:
    __slots__ = (
        "name",
        "total_ns",
        "self_ns",
        "count",
        "max_ns",
        "total_after_ns",
        "self_after_ns",
        "count_after",
        "max_after_ns",
    )

    def __init__(self, name=""):
        self.name = name
        self.total_ns = 0
        self.self_ns = 0
        self.count = 0
        self.max_ns = 0
        self.total_after_ns = 0
        self.self_after_ns = 0
        self.count_after = 0
        self.max_after_ns = 0


class _DependNodeStats:
    __slots__ = ("name", "owner", "total_ns", "self_ns", "count", "max_ns")

    def __init__(self, name="", owner=""):
        self.name = name
        self.owner = owner
        self.total_ns = 0
        self.self_ns = 0
        self.count = 0
        self.max_ns = 0


class _DependEdgeStats:
    __slots__ = (
        "total_ns",
        "count",
        "max_ns",
        "total_after_ns",
        "count_after",
        "max_after_ns",
    )

    def __init__(self):
        self.total_ns = 0
        self.count = 0
        self.max_ns = 0
        self.total_after_ns = 0
        self.count_after = 0
        self.max_after_ns = 0


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
            self._all_times = defaultdict(_TimerStats)
            self._thread_state = threading.local()
            self._timer_nodes = {}
            self._timer_edges = defaultdict(_DependEdgeStats)
            self._depend_nodes = {}
            self._depend_edges = defaultdict(_DependEdgeStats)
            self.decimals = 4
            self.skip_first_default = max(DEFAULT_SKIP_FIRST, 0)
            self.current_step = 0
            self.steps_seen = 0
        else:
            # empty functions in isabled mode to avoid overhead
            self.start = lambda *a, **k: None
            self.stop = lambda *a, **k: None
            self.set_step = lambda *a, **k: None
            self.enter_depend = lambda *a, **k: None
            self.exit_depend = lambda *a, **k: None

    def _get_timer_stack(self):
        stack = getattr(self._thread_state, "timer_stack", None)
        if stack is None:
            stack = []
            self._thread_state.timer_stack = stack
        return stack

    def _get_depend_stack(self):
        stack = getattr(self._thread_state, "depend_stack", None)
        if stack is None:
            stack = []
            self._thread_state.depend_stack = stack
        return stack

    def start(self, name):
        """Start or restart a named timer."""
        if name not in self._timer_nodes:
            self._timer_nodes[name] = _TimerNodeStats(name=name)
        self._get_timer_stack().append(
            {
                "name": name,
                "start_ns": time.perf_counter_ns(),
                "child_ns": 0,
            }
        )

    def stop(self, name):
        """Stop a started timer and record elapsed milliseconds.

        Args:
            name: Timer name passed to :meth:`start`.

        Returns:
            Elapsed time in milliseconds as ``np.float32``.

        Raises:
            ValueError: If ``start(name)`` was not called first.
        """
        timer_stack = self._get_timer_stack()
        if not timer_stack:
            raise ValueError(f"Timer '{name}' was not started.")
        frame = timer_stack.pop()
        if frame["name"] != name:
            raise ValueError(
                f"Timer stack mismatch: stopping '{name}' but active timer is '{frame['name']}'."
            )
        elapsed_ns = time.perf_counter_ns() - frame["start_ns"]
        self_ns = elapsed_ns - frame["child_ns"]
        elapsed_ms = elapsed_ns * 1.0e-6
        stats = self._all_times[name]
        stats.sum += elapsed_ms
        stats.count += 1
        if stats.count == 1 or elapsed_ms > stats.max:
            stats.max = elapsed_ms

        next_step = self.current_step + 1
        if next_step > self.steps_seen:
            self.steps_seen = next_step
        node = self._timer_nodes[name]
        node.total_ns += elapsed_ns
        node.self_ns += self_ns
        node.count += 1
        if node.count == 1 or elapsed_ns > node.max_ns:
            node.max_ns = elapsed_ns
        if self.current_step >= self.skip_first_default:
            stats.sum_after += elapsed_ms
            stats.count_after += 1
            if stats.count_after == 1 or elapsed_ms > stats.max_after:
                stats.max_after = elapsed_ms
            node.total_after_ns += elapsed_ns
            node.self_after_ns += self_ns
            node.count_after += 1
            if node.count_after == 1 or elapsed_ns > node.max_after_ns:
                node.max_after_ns = elapsed_ns

        if timer_stack:
            parent = timer_stack[-1]
            parent["child_ns"] += elapsed_ns
            edge = self._timer_edges[(parent["name"], name)]
            edge.total_ns += elapsed_ns
            edge.count += 1
            if edge.count == 1 or elapsed_ns > edge.max_ns:
                edge.max_ns = elapsed_ns
            if self.current_step >= self.skip_first_default:
                edge.total_after_ns += elapsed_ns
                edge.count_after += 1
                if edge.count_after == 1 or elapsed_ns > edge.max_after_ns:
                    edge.max_after_ns = elapsed_ns
        return elapsed_ms

    def set_step(self, step):
        """Set the current simulation step for step-normalized statistics."""
        self.current_step = max(int(step), 0)
        self.steps_seen = max(self.steps_seen, self.current_step + 1)

    def _depend_node_id(self, dep):
        return id(dep)

    def _depend_node_meta(self, dep):
        name = getattr(dep, "_timing_label", None)
        if name is None:
            name = getattr(dep, "_name", type(dep).__name__)
        owner = type(dep).__name__
        return name, owner

    def enter_depend(self, dep):
        """Enter a depend refresh region.

        Returns an opaque token or ``None`` when depend timing is inactive for
        the current step.
        """
        if self.current_step < self.skip_first_default:
            return None

        node_id = self._depend_node_id(dep)
        if node_id not in self._depend_nodes:
            name, owner = self._depend_node_meta(dep)
            self._depend_nodes[node_id] = _DependNodeStats(name=name, owner=owner)

        frame = {
            "node_id": node_id,
            "start_ns": time.perf_counter_ns(),
            "child_ns": 0,
        }
        self._get_depend_stack().append(frame)
        return frame

    def exit_depend(self, token):
        """Exit a depend refresh region started with :meth:`enter_depend`."""
        if token is None:
            return 0.0
        depend_stack = self._get_depend_stack()
        frame = depend_stack.pop()
        if frame is not token:
            raise ValueError("Depend timing stack mismatch.")

        elapsed_ns = time.perf_counter_ns() - frame["start_ns"]
        self_ns = elapsed_ns - frame["child_ns"]
        node = self._depend_nodes[frame["node_id"]]
        node.total_ns += elapsed_ns
        node.self_ns += self_ns
        node.count += 1
        if elapsed_ns > node.max_ns:
            node.max_ns = elapsed_ns

        if depend_stack:
            parent = depend_stack[-1]
            parent["child_ns"] += elapsed_ns
            edge = self._depend_edges[(parent["node_id"], frame["node_id"])]
            edge.total_ns += elapsed_ns
            edge.count += 1
            if elapsed_ns > edge.max_ns:
                edge.max_ns = elapsed_ns
        return elapsed_ns * 1.0e-6

    def get_total(self, name):
        """Return total measured time (ms) for ``name`` across all samples."""
        return round(self._all_times[name].sum, self.decimals)

    def get_count(self, name):
        """Return number of recorded samples for ``name``."""
        return self._all_times[name].count

    def get_average(self, name, skip_first=0):
        """Return average time (ms) for ``name``, optionally skipping warmup.

        Args:
            name: Timer name.
            skip_first: If ``> 0``, use the configured "after N steps" average.
        """
        if skip_first <= 0:
            steps = max(self.steps_seen, 1)
            return round(self._all_times[name].sum / steps, self.decimals)
        steps_after = max(self.steps_seen - skip_first, 0)
        if steps_after == 0:
            return 0.0
        return round(self._all_times[name].sum_after / steps_after, self.decimals)

    def get_max(self, name, skip_first=0):
        """Return maximum recorded value (ms) for ``name``.

        Args:
            name: Timer name.
            skip_first: If ``> 0``, use the configured "after N steps" max.
        """
        stats = self._all_times[name]
        if skip_first <= 0:
            return stats.max
        return stats.max_after

    def reset(self, name=None):
        """Clear measurements for one timer or all timers.

        Args:
            name: If provided, clear only this timer. If ``None``, clear all
                running and recorded timers.
        """
        if name is None:
            self._all_times.clear()
            self._get_timer_stack().clear()
            self._timer_nodes.clear()
            self._timer_edges.clear()
            self._get_depend_stack().clear()
            self._depend_nodes.clear()
            self._depend_edges.clear()
            self.current_step = 0
            self.steps_seen = 0
        else:
            self._all_times.pop(name, None)
            self._timer_nodes.pop(name, None)

    def _display_name(self, name):
        """Return a user-facing timer name."""
        return str(name).strip()

    def _summary_layout(self):
        name_w = 40
        avg_w = 20
        max_w = 20
        count_w = 8
        total_w = 12
        table_w = name_w + avg_w + max_w + count_w + total_w + 4
        header_fmt = (
            f"{{:<{name_w}}} {{:>{avg_w}}} {{:>{max_w}}} "
            f"{{:>{count_w}}} {{:>{total_w}}}"
        )
        row_fmt = (
            f"{{:<{name_w}}} {{:>{avg_w}}} {{:>{max_w}}} "
            f"{{:>{count_w}}} {{:>{total_w}}}"
        )
        return {
            "name_w": name_w,
            "avg_w": avg_w,
            "max_w": max_w,
            "count_w": count_w,
            "total_w": total_w,
            "table_w": table_w,
            "header_fmt": header_fmt,
            "row_fmt": row_fmt,
        }

    def _summary_header_lines(
        self, title=None, avg_header="Avg(ms/MDstep)", max_header="Max(ms/MDstep)"
    ):
        layout = self._summary_layout()
        lines = []
        if title is None:
            lines.append("\n")
        else:
            lines.extend(["", title])
        lines.append(
            layout["header_fmt"].format(
                "Timer",
                avg_header,
                max_header,
                "Count",
                "Total(ms)",
            )
        )
        lines.append("-" * layout["table_w"])
        return lines

    def _build_tree_prefix(self, ancestor_has_more, is_last):
        if is_last is None:
            return "", 0
        stem = "".join("│   " if has_more else "    " for has_more in ancestor_has_more)
        branch = "└─ " if is_last else "├─ "
        indent_width = 4 * len(ancestor_has_more)
        return stem + branch, indent_width

    def _numeric_cell(self, value, width, prefix="", is_int=False):
        if is_int:
            value_str = f"{value:d}"
        else:
            value_str = f"{value:.3f}"

        if len(prefix) >= width:
            return prefix[:width]
        return f"{prefix}{value_str:>{width - len(prefix)}}"

    def _format_summary_row(
        self,
        rendered_name,
        avg_value,
        max_value,
        count_value,
        total_value,
        tree_numbers=False,
        indent_width=0,
    ):
        layout = self._summary_layout()
        indent_prefix = (" " * indent_width) if tree_numbers else ""
        return layout["row_fmt"].format(
            rendered_name,
            self._numeric_cell(avg_value, layout["avg_w"], indent_prefix, is_int=False),
            self._numeric_cell(max_value, layout["max_w"], "", is_int=False),
            self._numeric_cell(count_value, layout["count_w"], "", is_int=True),
            self._numeric_cell(total_value, layout["total_w"], "", is_int=False),
        )

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

        incoming = set()
        edge_children = defaultdict(list)
        for (parent_name, child_name), edge in self._timer_edges.items():
            incoming.add(child_name)
            edge_children[parent_name].append((child_name, edge))
        for parent_name in edge_children:
            edge_children[parent_name].sort(
                key=lambda item: item[1].total_ns, reverse=True
            )

        roots = [name for name in self._timer_nodes.keys() if name not in incoming]
        roots.sort(key=lambda name: self._timer_nodes[name].total_ns, reverse=True)
        global_total = sum(self.get_total(name) for name in roots)
        if global_total <= 0.0:
            global_total = sum(
                self.get_total(name) for name in self._timer_nodes.keys()
            )

        # Helper to register a row
        def record_row(
            name,
            avg_skip,
            avg_all,
            max_skip,
            max_all,
            count,
            total,
            indent_level,
            parent_total,
        ):
            total_pct = (total / global_total * 100) if global_total > 0 else 0.0
            parent_pct = (total / parent_total * 100) if parent_total > 0 else 0.0

            rows.append(
                {
                    "Name": name,
                    "Indent": indent_level,
                    "Avg(ms/MDstep)": avg_skip,
                    "Total%": total_pct,
                    "Parent%": parent_pct,
                    "Avg_all(ms/MDstep)": avg_all,
                    "Max(ms/MDstep)": max_skip,
                    "Max_all(ms/MDstep)": max_all,
                    "Count": count,
                    "Total(ms)": total,
                }
            )

        # ---------------- Print identical header ----------------
        for line in self._summary_header_lines():
            print(line)
        visited = set()

        # Local helper for printing line AND recording DF row
        def print_line(
            name,
            skip_first,
            indent_level,
            parent_name=None,
            ancestor_has_more=None,
            is_last=None,
        ):
            if ancestor_has_more is None:
                ancestor_has_more = []
            avg_all = self.get_average(name)
            avg_skip = self.get_average(name, skip_first)
            max_all = self.get_max(name)
            max_skip = self.get_max(name, skip_first)
            total = self.get_total(name)
            count = self.get_count(name)
            # for some reason indented number print doesnt work yet
            tree_prefix, indent_width = self._build_tree_prefix(
                ancestor_has_more, is_last
            )
            indent_width = len(tree_prefix)
            rendered_name = f"{tree_prefix}{self._display_name(name)}"
            print(
                self._format_summary_row(
                    rendered_name,
                    avg_skip,
                    max_skip,
                    count,
                    total,
                    tree_numbers=tree_numbers,
                    indent_width=indent_width,
                )
            )

            parent_total = self.get_total(parent_name) if parent_name else total
            record_row(
                name,
                avg_skip,
                avg_all,
                max_skip,
                max_all,
                count,
                total,
                indent_level,
                parent_total,
            )

        def print_children(parent_name, indent_level, ancestor_has_more):
            parent_node = self._timer_nodes[parent_name]
            children = edge_children.get(parent_name, [])
            child_ancestor_has_more = ancestor_has_more

            for idx, (child_name, _) in enumerate(children):
                visited.add(child_name)
                is_last_child = idx == len(children) - 1
                print_line(
                    child_name,
                    skip_first,
                    indent_level,
                    parent_name,
                    ancestor_has_more=child_ancestor_has_more,
                    is_last=is_last_child,
                )
                next_ancestor_has_more = child_ancestor_has_more + [not is_last_child]
                print_children(
                    child_name,
                    indent_level + 1,
                    ancestor_has_more=next_ancestor_has_more,
                )

            if children:
                steps_after = max(self.steps_seen - skip_first, 0)
                remainder = (
                    (parent_node.self_after_ns * 1.0e-6) / steps_after
                    if steps_after > 0
                    else 0.0
                )
                parent_avg = self.get_average(parent_name, skip_first)
                if parent_avg > 0:
                    percent = f"Remainder: {remainder/parent_avg*100:3.1f}%"
                else:
                    percent = "Remainder: 0.0%"
                remainder_indent = "".join(
                    "│   " if has_more else "    " for has_more in ancestor_has_more
                )
                print(
                    self._format_summary_row(
                        remainder_indent + percent,
                        remainder,
                        0.0,
                        0,
                        0.0,
                        tree_numbers=tree_numbers,
                        indent_width=4 * len(ancestor_has_more),
                    )
                )

        for name in roots:
            if name in visited:
                continue
            visited.add(name)
            print_line(name, skip_first, indent_level=0, parent_name=None)
            print_children(name, indent_level=1, ancestor_has_more=[])

        for name in sorted(self._timer_nodes.keys()):
            if name in visited:
                continue
            visited.add(name)
            print_line(name, skip_first, indent_level=0, parent_name=None)
            print_children(name, indent_level=1, ancestor_has_more=[])

        # ---------------- Restore stdout ----------------
        sys.stdout = old_stdout

        # Write log file
        with open(f"{filename}.log", "w") as f:
            f.write(log_buffer.getvalue())

        # Write dataframe
        df = pd.DataFrame(rows)
        df.to_csv(f"{filename}.dat", sep=";", index=False, float_format="%.4f")

        # Print to console what user would normally see
        print(log_buffer.getvalue(), end="")

        self.summary_depend(filename)

    def summary_depend(self, filename):
        """Write and print a dynamic depend-refresh timing summary."""
        if not ENABLED or not self._depend_nodes:
            return None

        skip_first = self.skip_first_default
        roots = set(self._depend_nodes.keys())
        for _, child_id in self._depend_edges.keys():
            roots.discard(child_id)
        roots = sorted(
            roots, key=lambda nid: self._depend_nodes[nid].total_ns, reverse=True
        )

        edge_children = defaultdict(list)
        for (parent_id, child_id), edge in self._depend_edges.items():
            edge_children[parent_id].append((child_id, edge))
        for parent_id in edge_children:
            edge_children[parent_id].sort(
                key=lambda item: item[1].total_ns, reverse=True
            )

        rows = []
        lines = self._summary_header_lines(
            "Depend Refresh Timings",
            avg_header="Avg(ms/call)",
            max_header="Max(ms/call)",
        )

        global_total_ms = (
            sum(node.total_ns for node in self._depend_nodes.values()) * 1.0e-6
        )

        visited = set()

        def visit(
            node_id, parent_id=None, indent=0, ancestor_has_more=None, is_last=None
        ):
            if ancestor_has_more is None:
                ancestor_has_more = []
            node = self._depend_nodes[node_id]
            total_ms = node.total_ns * 1.0e-6
            self_ms = node.self_ns * 1.0e-6
            steps_after = max(self.steps_seen - skip_first, 0)
            avg_call_ms = total_ms / node.count if node.count > 0 else 0.0
            max_after_ms = node.max_ns * 1.0e-6

            tree_prefix, indent_width = self._build_tree_prefix(
                ancestor_has_more, is_last
            )
            label = f"{tree_prefix}{node.name}"
            lines.append(
                self._format_summary_row(
                    label,
                    avg_call_ms,
                    max_after_ms,
                    node.count,
                    total_ms,
                    indent_width=len(tree_prefix),
                )
            )

            parent_total_ms = (
                self._depend_nodes[parent_id].total_ns * 1.0e-6
                if parent_id is not None
                else total_ms
            )
            total_pct = (
                (total_ms / global_total_ms * 100.0) if global_total_ms > 0 else 0.0
            )
            parent_pct = (
                (total_ms / parent_total_ms * 100.0) if parent_total_ms > 0 else 0.0
            )
            rows.append(
                {
                    "Name": node.name,
                    "Indent": indent,
                    "Avg(ms/call)": round(avg_call_ms, 4),
                    "Total%": round(total_pct, 4),
                    "Parent%": round(parent_pct, 4),
                    "Avg_all(ms/call)": round(avg_call_ms, 4),
                    "Max(ms/call)": round(max_after_ms, 4),
                    "Max_all(ms/call)": round(max_after_ms, 4),
                    "Count": node.count,
                    "Total(ms)": round(total_ms, 4),
                    "Self(ms)": round(self_ms, 4),
                }
            )
            visited.add(node_id)
            children = edge_children.get(node_id, [])
            child_ancestor_has_more = (
                ancestor_has_more + [not is_last]
                if is_last is not None
                else ancestor_has_more
            )
            child_total_avg = 0.0
            for idx, (child_id, edge) in enumerate(children):
                is_last_child = idx == len(children) - 1
                child_total_avg += (
                    (edge.total_ns * 1.0e-6) / steps_after if steps_after > 0 else 0.0
                )
                visit(
                    child_id,
                    parent_id=node_id,
                    indent=indent + 1,
                    ancestor_has_more=child_ancestor_has_more,
                    is_last=is_last_child,
                )

            if children:
                remainder = self_ms / node.count if node.count > 0 else 0.0
                if avg_call_ms > 0:
                    percent = f"Remainder: {remainder/avg_call_ms*100:3.1f}%"
                else:
                    percent = "Remainder: 0.0%"
                remainder_indent = "".join(
                    "│   " if has_more else "    "
                    for has_more in child_ancestor_has_more
                )
                lines.append(
                    self._format_summary_row(
                        remainder_indent + percent,
                        remainder,
                        0.0,
                        0,
                        0.0,
                        indent_width=4 * len(child_ancestor_has_more),
                    )
                )

        for idx, node_id in enumerate(roots):
            visit(node_id, is_last=None)
        for node_id in sorted(self._depend_nodes.keys()):
            if node_id not in visited:
                visit(node_id, is_last=None)

        text = "\n".join(lines) + "\n"
        with open(f"{filename}_depend.log", "w") as f:
            f.write(text)
        pd.DataFrame(rows).to_csv(
            f"{filename}_depend.dat", sep=";", index=False, float_format="%.4f"
        )
        print(text, end="")


# Global singleton
timers = TimingManager()

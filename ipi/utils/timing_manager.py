"""
#use:
from timing_manager import timers

timers.start("Name") #start a timer
timers.stop("Name") #stop a timer

timers.summary() # print timings at the end, take averages skipping the first 5 measure by default
timers.summary(skip_first=1) # take averages skipping the first 1

#timer keys with the format [+]Name2 are automatically the children of Name1(+)
"""

import time
from collections import defaultdict
import numpy as np
import os
import pandas as pd
from io import StringIO
import sys

ENABLED = os.getenv("ENABLE_TIMING_MANAGER", "0") == "1"

class TimingManager:
    def __init__(self):
        if ENABLED:
            self._start_times = {}
            self._all_times = defaultdict(list)  # ms values, low precision
            self.group_symbols = [
                                    "+", "++", "+++", "++++","+++++","++++++",
                                    "*", "**", "***", "****",
                                    "#", "##", "###", "####"
                                  ]  # supported hierarchy roots
            self.decimals = 4
            self.majorstep = None
        else:
            # empty functions in isabled mode to avoid overhead
            self.start = lambda *a, **k: None
            self.stop = lambda *a, **k: None

    def start(self, name):
        self._start_times[name] = time.perf_counter()

    def stop(self, name):
        if name not in self._start_times:
            raise ValueError(f"Timer '{name}' was not started.")
        elapsed_sec = time.perf_counter() - self._start_times.pop(name)
        elapsed_ms = np.float32(round(elapsed_sec * 1000, self.decimals))
        self._all_times[name].append(elapsed_ms)
        return elapsed_ms

    def get_total(self, name):
        return round(sum(self._all_times[name]), self.decimals)

    def get_count(self, name):
        return len(self._all_times[name])

    def set_majorstep(self):
        # to handle multiple calls in a step, we get the major step number by getting the minimum of step numbers
        for i,name in enumerate(self._all_times.keys()):
            size = self.get_count(name)
            if i == 0:
                self.majorstep = size
            if size < self.majorstep:
                self.majorstep = size
        print(f"@timingmanager set majorstep for {self.majorstep}") 

    def get_average(self, name, skip_first=0):
        skip = self.get_count(name) // self.majorstep * skip_first # to handle multiple calls

        times = self._all_times[name][skip:]
        return round(sum(times) / len(times), self.decimals) if times else 0.0

    def get_max(self, name, skip_first=0):
        times = self._all_times[name][skip_first:]
        return max(times) if times else 0.0

    def reset(self, name=None):
        if name is None:
            self._start_times.clear()
            self._all_times.clear()
        else:
            self._start_times.pop(name, None)
            self._all_times.pop(name, None)

    def _parse_marker(self, name):
        """Return (base_name, marker_of_parent, marker_of_children) if '(marker)' present, else (name, None)."""
        marker_of_parent = None
        marker_of_children = None

        if "(" in name and name.endswith(")"):
            _, marker = name.rsplit("(", 1)
            marker_of_children =  marker[:-1]
        if "]" in name and name.startswith("["):
            marker, _ = name.rsplit("]", 1)
            marker_of_parent =  marker[1:]

        return marker_of_parent, marker_of_children



    def summary(self, skip_first=5):
        """Print hierarchical timing summary and write .log and .dat files."""

        if not ENABLED:
            return None

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
        self.set_majorstep()
        print("\n")
        print(f"{'Timer':<40} {'Avg_aftern(ms/step)':>15} "
              f"{'Avg(ms/step)':>12} "
              f"{'Max_aftern(ms/step)':>15} "
              f"{'Max(ms/step)':>12} "
              f"{'Count':>6} "
              f"{'Total(ms)':>12}")
        print("-" * 130)

        all_names = sorted(self._all_times.keys())
        visited = set()

        # Local helper for printing line AND recording DF row
        def print_line(name, skip_first, indent_level, parent_name=None):
            avg_all = self.get_average(name)
            avg_skip = self.get_average(name, skip_first)
            max_all = self.get_max(name)
            max_skip = self.get_max(name, skip_first)
            total = self.get_total(name)
            count = self.get_count(name)

            print(f"{'  ' * indent_level}{name:<40} {avg_skip:15.3f} "
                  f"{avg_all:12.3f}"
                  f"{max_skip:15.3f}"
                  f"{max_all:12.3f} "
                  f"{count:6d} "
                  f"{total:12.3f} ")

            parent_total = self.get_total(parent_name) if parent_name else total
            record_row(name, avg_skip, avg_all, max_skip, max_all, count,
                       total, indent_level, parent_total)


        # ---------------- Recursive children printer ----------------
        def print_children(parent_name, marker, skip_first, indent_level):
            direct_children = [
                n for n in self._all_times.keys()
                if self._parse_marker(n)[0] == marker
            ]
            direct_children.sort()

            total_child_time = 0.0

            for child in direct_children:
                visited.add(child)
                print_line(child, skip_first, indent_level, parent_name)

                count_ratio = self.get_count(child)/self.get_count(parent_name)
                total_child_time += self.get_average(child, skip_first)*count_ratio

                _, child_marker = self._parse_marker(child)
                if child_marker:
                    print_children(child, child_marker, skip_first, indent_level + 1)

            # Remainder line (console only)
            if direct_children:
                parent_avg = self.get_average(parent_name, skip_first)
                remainder = round(parent_avg - total_child_time, 3)
                percent = f"[{marker}]Remainder: {remainder/parent_avg*100:3.1f}%"
                print(f"{'  ' * indent_level}{percent:<40} {remainder:15.3f}")


        # ---------------- Top-level timers ----------------
        for name in all_names:
            if name in visited:
                continue

            marker_of_parent, marker_of_children = self._parse_marker(name)

            if marker_of_parent is None and marker_of_children:
                visited.add(name)
                print_line(name, skip_first, indent_level=0, parent_name=None)
                print_children(name, marker_of_children, skip_first,
                               indent_level=1)
            else:
                visited.add(name)
                print_line(name, skip_first, indent_level=0, parent_name=None)

        # ---------------- Restore stdout ----------------
        sys.stdout = old_stdout

        # Write log file
        with open("timing_breakdown.log", "w") as f:
            f.write(log_buffer.getvalue())

        # Write dataframe
        df = pd.DataFrame(rows)
        df.to_csv("timing_breakdown.dat", sep=";", index=False, float_format='%.4f')

        # Print to console what user would normally see
        print(log_buffer.getvalue(), end="")


# Global singleton
timers = TimingManager()


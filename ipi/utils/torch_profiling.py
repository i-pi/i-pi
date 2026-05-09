"""Optional torch.profiler integration.

Enabled by the ``i-pi -p torch`` command-line flag (see ``bin/i-pi``).
The CLI calls :func:`configure` before ``Simulation.run`` starts; everything
else keys off the module-level state set there. When the feature is not
configured all helpers are zero-cost no-ops, so the single call site in
``Simulation.run`` can be unconditional.

The profiler records CPU + CUDA activities and captures Python call
stacks for every op (``with_stack=True``). In the chrome trace every
``aten::*`` node shows the Python source line that launched it, so no
in-engine annotations are needed to find the host/device transfers.

Schedule is ``wait,warmup,active`` (see ``torch.profiler.schedule``):
``wait`` steps of the MD loop are skipped, ``warmup`` steps prime kernels,
and the next ``active`` steps are recorded into the chrome trace. The
window fires exactly once from the start of the run.
"""

__all__ = [
    "configure",
    "is_enabled",
    "build_profiler",
    "summarize",
]

_STATE = {
    "enabled": False,
    "schedule": (5, 5, 40),
    "prefix": "profile",
}


def _parse_schedule(spec):
    if spec is None:
        return (5, 5, 40)
    if isinstance(spec, (tuple, list)) and len(spec) == 3:
        return tuple(int(x) for x in spec)
    try:
        parts = [int(p.strip()) for p in str(spec).split(",")]
    except ValueError:
        raise ValueError(
            "torch profiler schedule must be 'wait,warmup,active' integers"
        )
    if len(parts) != 3:
        raise ValueError("torch profiler schedule must have exactly three integers")
    return tuple(parts)


def configure(enabled=True, schedule=None, prefix="profile"):
    """Arm the torch profiler. Called by the CLI before ``Simulation.run``."""
    _STATE["enabled"] = bool(enabled)
    _STATE["schedule"] = _parse_schedule(schedule)
    _STATE["prefix"] = prefix


def is_enabled():
    return _STATE["enabled"]


class _NoopProfiler:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def step(self):
        pass


def build_profiler():
    """Return ``(profiler, trace_path)``.

    When disabled returns a no-op profiler and ``None`` so callers can
    unconditionally use ``with prof:`` / ``prof.step()`` and branch on
    ``trace_path is not None`` for post-processing.
    """
    if not _STATE["enabled"]:
        return _NoopProfiler(), None
    try:
        import torch
        from torch.profiler import profile, ProfilerActivity, schedule
    except ImportError:
        return _NoopProfiler(), None
    wait, warmup, active = _STATE["schedule"]
    activities = [ProfilerActivity.CPU]
    if torch.cuda.is_available():
        activities.append(ProfilerActivity.CUDA)
    prof = profile(
        activities=activities,
        schedule=schedule(wait=wait, warmup=warmup, active=active, repeat=1),
        record_shapes=False,
        with_stack=True,
    )
    trace_path = _STATE["prefix"] + ".torch_trace.json"
    return prof, trace_path


def summarize(prof, trace_path):
    """Export the chrome trace and return top-N tables as a list of strings."""
    lines = []
    try:
        prof.export_chrome_trace(trace_path)
        lines.append(f" @torch_profiler: chrome trace written to {trace_path}")
    except Exception as exc:
        lines.append(f" @torch_profiler: failed to export trace: {exc}")
    try:
        import torch

        if torch.cuda.is_available():
            lines.append(" @torch_profiler: top ops by self CUDA time")
            lines.append(
                prof.key_averages().table(sort_by="self_cuda_time_total", row_limit=25)
            )
        lines.append(" @torch_profiler: top ops by self CPU time")
        lines.append(
            prof.key_averages().table(sort_by="self_cpu_time_total", row_limit=25)
        )
    except Exception as exc:
        lines.append(f" @torch_profiler: failed to format tables: {exc}")
    return lines

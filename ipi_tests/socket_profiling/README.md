A technical socket profiler testsuite for high-performance PIMD runs, with 32 beads.

Timers for timing_manager are currently set up for an NVT run, to measure NPT correctly, you need to add timers into the NPT integrator. 

To measure communication overhead, the Waiting Driver timings, and the corresponding Remainder (not measured explicitly by timers, but significant due to python looping) should be added together.

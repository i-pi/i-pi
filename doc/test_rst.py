import sys 

src_dir = "../"
sys.path.insert(0, src_dir)

from ipi.inputs import simulation

sim = simulation.InputSimulation()

print(sim.help_rst(name="simulation", stop_level=4))

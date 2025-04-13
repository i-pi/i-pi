import numpy
from ipi import read_output

suffixes = ["_species", "_indices"]

for suffix in suffixes:
    data = read_output("simulation.temp" + suffix)
    for key in data[0].keys():
        print(
            "%.30s :  %10.5f +/- %10.5f"
            % (
                key,
                numpy.mean(data[0][key]),
                2 * numpy.std(data[0][key]) / len(data[0][key]) ** 0.5,
            )
        )

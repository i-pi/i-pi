from pylab import *

idtdata = []
# output is a standard output from i-pi
for line in open("output"):
    llist = line.split()
    if llist:
        if llist[0] == "ThermoCL":
            idtdata.append(float(llist[5]))

nsteps = len(idtdata)
maxt = nsteps * 0.1
print(
    (
        "Mean IDTAU for "
        + str(maxt * 0.5)
        + "<t<"
        + str(maxt)
        + " fs: "
        + str(mean(idtdata[nsteps // 2 : nsteps]))
    )
)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(linspace(0, maxt, num=nsteps), idtdata)
ax.set(
    title="Automatic IDTAU Parameter Adjustment",
    xlabel="Time [fs]",
    ylabel="IDTAU [a.u.]",
    xlim=(0, maxt),
)
show()

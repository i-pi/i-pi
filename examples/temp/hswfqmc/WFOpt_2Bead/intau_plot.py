from pylab import *

invdata = []
for line in open("output"):
    llist = line.split()
    if llist:
        if llist[0] == "ThermoCL":
            invdata.append(float(llist[5]))

nsteps = len(invdata)
maxt = nsteps * 0.1
print(
    (
        "Mean INTAU for "
        + str(maxt * 0.5)
        + "<t<"
        + str(maxt)
        + " fs: "
        + str(mean(invdata[nsteps // 2 : nsteps]))
    )
)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(linspace(0, maxt, num=nsteps), invdata)
ax.set(
    title="Automatic INTAU Parameter Adjustment",
    xlabel="Time [fs]",
    ylabel="INTAU [a.u.]",
    xlim=(0, maxt),
)
show()

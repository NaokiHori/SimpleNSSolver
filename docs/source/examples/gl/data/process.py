import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(dname, ras, xmin):
    xs = list()
    ys = list()
    for ra in ras:
        fname = f"{dname}/nusselt_{ra}.dat"
        data = np.loadtxt(fname, usecols=[0, 1])
        x = data[:, 0]
        y = data[:, 1]
        y = y[x > xmin]
        xs.append(float(ra))
        ys.append(np.average(y))
    return xs, ys

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    ras = ["1.0e+4", "3.1e+4", "1.0e+5", "3.1e+5", "1.0e+6", "3.1e+6", "1.0e+7", "3.1e+7", "1.0e+8"]
    xmin = 100.
    idname = argv[1]
    ofname = argv[2]
    x, y = load(idname, ras, xmin)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, color="#FF0000", marker=".")
    # Kooloth et al., PRF, 2021
    ax.plot(x, 0.1 * np.power(x, 0.28), color="#000000", linestyle="--")
    kwrds = {
            "title": "",
            "xlim": [0.5 * np.min(x), 2. * np.max(x)],
            "xlabel": "Ra",
            "ylabel": "Nu",
            "xticks": x,
            "yticks": [1, 2, 4, 8, 16, 32],
            "xscale": "log",
            "yscale": "log",
    }
    ax.set(**kwrds)
    pyplot.savefig(ofname)
    pyplot.close()


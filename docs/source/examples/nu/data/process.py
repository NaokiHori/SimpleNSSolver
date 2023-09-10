import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname):
    data = np.loadtxt(fname)
    x  = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    y2 = data[:, 3]
    y3 = data[:, 4]
    return x, y0, y1, y2, y3

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    ifname = argv[1]
    ofname = argv[2]
    x, y0, y1, y2, y3 = load(ifname)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, np.abs(y1 - y0), color="#FF0000")
    ax.plot(x, np.abs(y2 - y0), color="#0000FF")
    ax.plot(x, np.abs(y3 - y0), color="#33AA00")
    kwrds = {
            "title": "",
            "xlabel": "time",
            "ylabel": "deviations",
            "yticks": [1.e-16, 1.e-12, 1.e-8, 1.e-4, 1.e+0],
            "yscale": "log",
    }
    ax.set(**kwrds)
    pyplot.savefig(ofname)
    pyplot.close()


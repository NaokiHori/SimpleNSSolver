import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname, xmin):
    data = np.loadtxt(fname)
    x  = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    y2 = data[:, 3]
    y3 = data[:, 4]
    y0 = y0[x > xmin]
    y1 = y1[x > xmin]
    y2 = y2[x > xmin]
    y3 = y3[x > xmin]
    x  =  x[x > xmin]
    return x, y0, y1, y2, y3

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    ifname = argv[1]
    ofname = argv[2]
    is2d = True if "2" in ofname else False
    xmin = 200.
    x, y0, y1, y2, y3 = load(ifname, xmin)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y0, color="#FF0000", label="heat flux")
    ax.plot(x, y1, color="#0000FF", label="injection")
    ax.plot(x, y2, color="#33AA00", label="kinetic dissipation")
    ax.plot(x, y3, color="#FF00FF", label="thermal dissipation")
    if is2d:
        # van der Poel et al., JFM, 2013
        y4 = np.full(x.shape, 27.25)
        ax.plot(x, y4, color="#000000", linestyle="--")
    kwrds = {
            "title": "",
            "xlim": [xmin, 300.],
            "xlabel": "time",
            "ylabel": "Nusselt numbers",
            "xticks": [200., 250., 300.],
    }
    ax.set(**kwrds)
    ax.legend()
    pyplot.savefig(ofname)
    pyplot.close()


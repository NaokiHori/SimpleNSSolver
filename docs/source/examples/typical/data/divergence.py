import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname, xmin):
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
    y = y[x > xmin]
    x = x[x > xmin]
    return x, y

if __name__ == "__main__":
    argv = sys.argv
    assert 3 == len(argv)
    ifname = argv[1]
    ofname = argv[2]
    xmin = 200.
    x, y = load(ifname, xmin)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, color="#FF0000")
    kwrds = {
            "title": "",
            "xlim": [200., 300.],
            "ylim": [1.e-16, 1.e-13],
            "xlabel": "time",
            "ylabel": "maximum divergence",
            "xticks": [200., 250., 300.],
            "yticks": [1.e-16, 1.e-15, 1.e-14, 1.e-13],
            "yscale": "log",
    }
    ax.set(**kwrds)
    pyplot.savefig(ofname)
    pyplot.close()


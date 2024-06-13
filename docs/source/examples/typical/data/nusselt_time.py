import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname, xmin):
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
    y = y[xmin < x]
    x = x[xmin < x]
    return x, y

def get_ref_heat_transfer():
    ly = os.environ.get("ly")
    lz = os.environ.get("lz")
    Ra = os.environ.get("Ra")
    Pr = os.environ.get("Pr")
    ly = float(ly)
    lz = float(lz) if lz else 1.
    Ra = float(Ra)
    Pr = float(Pr)
    return ly * lz / Ra**0.5 / Pr**0.5

if __name__ == "__main__":
    argv = sys.argv
    assert 6 == len(argv)
    ifnames = argv[1:-1]
    ofname = argv[-1]
    is2d = True if None == os.environ.get("lz") else False
    xmin = 200.
    q_ref = get_ref_heat_transfer()
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    colors = [
            "#FF0000",
            "#0000FF",
            "#33AA00",
            "#FF00FF",
    ]
    labels = [
            "heat transfer",
            "kinetic energy injection",
            "kinetic energy dissipation",
            "thermal energy dissipation",
    ]
    converters = [
            lambda x: x / q_ref,
            lambda x: (x + q_ref) / q_ref,
            lambda x: (x + q_ref) / q_ref,
            lambda x: x / q_ref,
    ]
    for color, label, converter, ifname in zip(colors, labels, converters, ifnames):
        x, y = load(ifname, xmin)
        y = converter(y)
        ax.plot(x, y, color=color, label=label)
    if is2d:
        # van der Poel et al., JFM, 2013
        y = np.full(x.shape, 27.25)
        ax.plot(x, y, color="#000000", linestyle="dashed")
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


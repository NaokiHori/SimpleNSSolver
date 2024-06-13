import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname):
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
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
    xs = list()
    ys = list()
    for ifname in ifnames:
        x, y = load(ifname)
        xs.append(x)
        ys.append(y)
    q_ref = get_ref_heat_transfer()
    converters = [
            lambda x: x / q_ref,
            lambda x: (x + q_ref) / q_ref,
            lambda x: (x + q_ref) / q_ref,
            lambda x: x / q_ref,
    ]
    for cnt, converter in enumerate(converters):
        ys[cnt] = converter(ys[cnt])
    fig = pyplot.figure()
    ax = fig.add_subplot()
    ax.plot(xs[0], np.abs(ys[1] - ys[0]), color="#FF0000")
    ax.plot(xs[0], np.abs(ys[2] - ys[0]), color="#0000FF")
    ax.plot(xs[0], np.abs(ys[3] - ys[0]), color="#33AA00")
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


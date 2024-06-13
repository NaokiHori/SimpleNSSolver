import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(dname):
    fnames = [f"{dname}/{fname}" for fname in os.listdir(dname) if fname.startswith("energy") and fname.endswith(".dat")]
    fnames = sorted(fnames)
    ls  = list()
    xs  = list()
    y1s = list()
    y2s = list()
    for fname in fnames:
        data = np.loadtxt(fname)
        is3d = 5 == len(data[0])
        t  = data[:, 0]
        if is3d:
            kx = data[:, 1]
            ky = data[:, 2]
            kz = data[:, 3]
            k  = kx + ky + kz
            h  = data[:, 4]
            kref = k[0]
            href = h[0]
        else:
            kx = data[:, 1]
            ky = data[:, 2]
            k  = kx + ky
            h  = data[:, 3]
            kref = k[0]
            href = h[0]
        ls.append(fname.strip().split("energy-")[1].split(".dat")[0])
        xs.append(t[1:])
        y1s.append(k[1:] - kref)
        y2s.append(h[1:] - href)
    return ls, xs, y1s, y2s

if __name__ == "__main__":
    argv = sys.argv
    assert 4 == len(argv)
    idname  = argv[1]
    ofname1 = argv[2]
    ofname2 = argv[3]
    ls, xs, y1s, y2s = load(idname)
    colors = pyplot.rcParams["axes.prop_cycle"].by_key()["color"]
    #
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    for l, x, y1, y2, color in zip(ls, xs, y1s, y2s, colors):
        ax.plot(x, -y1, color=color, linestyle="-",  label=l)
        ax.plot(x, -y2, color=color, linestyle="--")
    kwrds = {
            "title": "",
            "xlabel": "time",
            "ylabel": "$-\Delta K$, $-\Delta H$",
            "yscale": "log",
    }
    ax.legend()
    ax.set(**kwrds)
    pyplot.savefig(ofname1)
    pyplot.close()
    #
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    xs  = list()
    z1s = list()
    z2s = list()
    for l, y1, y2 in zip(ls, y1s, y2s):
        xs.append(float(l))
        z1s.append(-y1[-1])
        z2s.append(-y2[-1])
    xs  = np.array(xs)
    z1s = np.array(z1s)
    z2s = np.array(z2s)
    ax.plot(xs, z1s, marker="o", label="$-\Delta K$", color="#FF0000")
    ax.plot(xs, z2s, marker="o", label="$-\Delta H$", color="#0000FF")
    ax.plot(xs, 2. * z1s[0] / xs[0]**3. * xs**3., label="3rd order", color="#000000")
    kwrds = {
            "title": "",
            "xlabel": "$\Delta t$",
            "ylabel": "$-\Delta K$, $-\Delta H$",
            "xscale": "log",
            "yscale": "log",
    }
    ax.legend()
    ax.set(**kwrds)
    pyplot.savefig(ofname2)
    pyplot.close()


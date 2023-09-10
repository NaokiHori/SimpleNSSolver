import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(dname):
    lengths = np.load(f"{dname}/lengths.npy")
    glsizes = np.load(f"{dname}/glsizes.npy")
    xc = np.load(f"{dname}/xc.npy")
    lx = lengths[0]
    ly = lengths[1]
    ny = glsizes[1]
    dy = ly / ny
    yc = np.linspace(0.5 * dy, ly - 0.5 * dy, ny)
    t = np.load(f"{dname}/t.npy")
    return lx, ly, xc, yc, t

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    idname = argv[1]
    ofname = argv[2]
    lx, ly, x, y, t = load(idname)
    # visualise transposed array since my screen is wider
    fig = pyplot.figure(figsize=(8, 4))
    fig.set_size_inches(8, 4)
    ax = pyplot.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.contourf(y, x, t.T, vmin=-0.5, vmax=+0.5, cmap="bwr", levels=51)
    kwrds = {
            "title": "",
            "aspect": "equal",
            "xlim": [0.0, ly],
            "ylim": [0.0, lx],
            "xlabel": "",
            "ylabel": "",
            "xticks": [],
            "yticks": [],
    }
    ax.set(**kwrds)
    pyplot.savefig(ofname, dpi=150)
    pyplot.close()


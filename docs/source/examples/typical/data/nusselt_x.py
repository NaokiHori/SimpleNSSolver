import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load_and_process(dname):
    coef    = np.load(f"{dname}/t_dif.npy")
    xf      = np.load(f"{dname}/xf.npy")
    xc      = np.load(f"{dname}/xc.npy")
    num     = np.load(f"{dname}/num.npy")
    t1      = np.load(f"{dname}/t1.npy")
    uxt     = np.load(f"{dname}/uxt.npy")
    glsizes = np.load(f"{dname}/glsizes.npy")
    if 3 == len(glsizes):
        ngrids = 1. * glsizes[1] * glsizes[2]
    else:
        ngrids = 1. * glsizes[1]
    # divide by the number of grids in the homogeneous direction(s)
    t1  /= ngrids
    uxt /= ngrids
    # divide by the number of samples
    t1  /= num
    uxt /= num
    # adv
    adv = 1. / coef * uxt
    # dif
    dif = -1. * np.diff(t1) / np.diff(xc)
    # average [0:1/2] and [1/2:1]
    nx = len(xf) // 2
    adv = 0.5 * adv + 0.5 * adv[::-1]
    dif = 0.5 * dif + 0.5 * dif[::-1]
    return xf[:nx+1], adv[:nx+1], dif[:nx+1]

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    idname = argv[1]
    ofname = argv[2]
    x, adv, dif = load_and_process(idname)
    tot = adv + dif
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, adv, color="#FF0000", linestyle="-", label="advection")
    ax.plot(x, dif, color="#0000FF", linestyle="-", label="diffusion")
    ax.plot(x, tot, color="#000000", linestyle=":", label="total")
    kwrds = {
            "title": "",
            "xlim": [0., 0.5],
            "xlabel": "wall-normal position",
            "ylabel": "Nusselt contributions",
            "xticks": [0., 0.25, 0.5],
    }
    ax.set(**kwrds)
    ax.legend()
    pyplot.savefig(ofname, dpi=150)
    pyplot.close()


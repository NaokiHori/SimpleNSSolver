import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

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

def load_and_process(dname):
    xf  = np.load(f"{dname}/xf.npy")
    num = np.load(f"{dname}/num.npy")
    adv = np.load(f"{dname}/adv.npy")
    dif = np.load(f"{dname}/dif.npy")
    is3d = 3 == len(adv.shape)
    # normalise
    adv /= get_ref_heat_transfer()
    dif /= get_ref_heat_transfer()
    # divide by the number of samples
    adv /= num
    dif /= num
    # average in the homogeneous direction(s)
    adv = np.sum(adv, axis=0)
    dif = np.sum(dif, axis=0)
    if is3d:
        adv = np.sum(adv, axis=0)
        dif = np.sum(dif, axis=0)
    # average bottom-half and top-half
    nx = len(xf) // 2
    adv = 0.5 * adv + 0.5 * adv[::-1]
    dif = 0.5 * dif + 0.5 * dif[::-1]
    return xf[:nx+1], adv[:nx+1], dif[:nx+1]

if __name__ == "__main__":
    argv = sys.argv
    assert 3 == len(argv)
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


import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(dname):
    glsizes = np.load(f"{dname}/glsizes.npy")
    is3d = True if 3 == len(glsizes) else False
    if is3d:
        ngrids = 1. * glsizes[1] * glsizes[2]
    else:
        ngrids = 1. * glsizes[1]
    num = np.load(f"{dname}/num.npy")
    xf  = np.load(f"{dname}/xf.npy")
    xc  = np.load(f"{dname}/xc.npy")
    ux1 = np.load(f"{dname}/ux1.npy")
    ux2 = np.load(f"{dname}/ux2.npy")
    uy1 = np.load(f"{dname}/uy1.npy")
    uy2 = np.load(f"{dname}/uy2.npy")
    if is3d:
        uz1 = np.load(f"{dname}/uz1.npy")
        uz2 = np.load(f"{dname}/uz2.npy")
    else:
        uz1 = None
        uz2 = None
    t1 = np.load(f"{dname}/t1.npy")
    t2 = np.load(f"{dname}/t2.npy")
    return is3d, num * ngrids, xf, xc, ux1, ux2, uy1, uy2, uz1, uz2, t1, t2

def compute_rms(data1, data2):
    rms2 = data2 - np.power(data1, 2.)
    rms2[rms2 < 0.] = 0.
    return np.sqrt(rms2)

def trunc(arr):
    n = len(arr)
    arr = arr[:n//2+1]
    return arr

if __name__ == "__main__":
    idname = sys.argv[1]
    ofname = sys.argv[2]
    is3d, num, xf, xc, ux1, ux2, uy1, uy2, uz1, uz2, t1, t2 = load(idname)
    #
    t1  /= num
    ux1 /= num
    uy1 /= num
    t2  /= num
    ux2 /= num
    uy2 /= num
    if is3d:
        uz1 /= num
        uz2 /= num
    #
    t2  = compute_rms(t1, t2)
    ux2 = compute_rms(ux1, ux2)
    uy2 = compute_rms(uy1, uy2)
    if is3d:
        uz2   = compute_rms(uz1, uz2)
    t2 = 0.5*(t2[:]+t2[::-1])
    ux2   = 0.5*(ux2[:]+ux2[::-1])
    uy2   = 0.5*(uy2[:]+uy2[::-1])
    if is3d:
        uz2   = 0.5*(uz2[:]+uz2[::-1])
    #
    xf  = trunc(xf)
    xc  = trunc(xc)
    t2  = trunc(t2)
    ux2 = trunc(ux2)
    uy2 = trunc(uy2)
    if is3d:
        uz2 = trunc(uz2)
    #
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(xf, ux2, color="#FF0000", label="ux")
    ax.plot(xc, uy2, color="#0000FF", label="uy")
    if is3d:
        ax.plot(xc, uz2, color="#FF00FF", label="uz")
    ax.plot(xc, t2, color="#33AA00", label="T")
    kwrds = {
            "title": "",
            "xlim": [0., 0.5],
            "xlabel": "wall-normal position",
            "ylabel": "standard deviations",
            "xticks": [0., 0.25, 0.5],
    }
    ax.set(**kwrds)
    ax.legend()
    pyplot.savefig(ofname, dpi=150)
    pyplot.close()


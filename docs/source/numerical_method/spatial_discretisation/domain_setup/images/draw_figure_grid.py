import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import patches

def gen_xf(glsize):
    xf = np.zeros(glsize+1)
    for i in range(glsize):
        xf[i + 1] = xf[i] \
            + 1. - 0.15 * np.cos(2. * np.pi * (2 * i + 1) / (2. * glsize))
    return xf

def gen_xc(glsize):
    xf = gen_xf(glsize)
    xc = 0.5 * xf[1:] + 0.5 * xf[:-1]
    xc = np.append(xf[0], xc)
    xc = np.append(xc, xf[-1])
    xc = np.delete(xc, [4, 8])
    return xc

def general(ax, glsize, ys):
    xf = gen_xf(glsize)
    # vertical lines
    ax.vlines(xf, ys[0], ys[1], colors="#000000")
    # horizontal lines
    for cnt, (x0, x1) in enumerate(zip(xf[:-1], xf[1:])):
        ls = "--" if 3 == cnt or 7 == cnt else "-"
        ax.hlines(ys, x0, x1, colors="#000000", linestyles=ls)
    keywords = {
            "aspect": "equal",
            "xlabel": "",
            "ylabel": "",
            "xlim": [np.min(xf)-1.0, np.max(xf)+1.0],
            "ylim": [np.min(ys)-1.5, np.max(ys)+0.5],
    }
    ax.set(**keywords)

def grid1(ax, glsize, ys):
    # xf
    xf = gen_xf(glsize)
    ax.vlines(xf, ys[0], ys[1], colors="#FF0000", linewidth=4)
    labels = (
            "XF(1)",
            "XF(2)",
            "XF(3)",
            "XF(4)",
            "XF(i-1)",
            "XF(i)",
            "XF(i+1)",
            "XF(i+2)",
            "XF(size-2)",
            "XF(size-1)",
            "XF(size)",
            "XF(size+1)",
    )
    for cnt, (x, label) in enumerate(zip(xf, labels)):
        ax.text(x, ys[0] - 1.0 + 0.5 * (cnt % 2), label, ha="center", va="center")

def grid2(ax, glsize, ys):
    # xc
    xc = gen_xc(glsize)
    ax.scatter(xc, np.full(xc.shape, 0.5 * ys[0] + 0.5 * ys[1]), color="#FF0000")
    labels = (
            "XC(0)",
            "XC(1)",
            "XC(2)",
            "XC(3)",
            "XC(i-1)",
            "XC(i)",
            "XC(i+1)",
            "XC(size-2)",
            "XC(size-1)",
            "XC(size)",
            "XC(size+1)",
    )
    for cnt, (x, label) in enumerate(zip(xc, labels)):
        va = "top"    if len(xc) - 3 == cnt else "center"
        va = "bottom" if len(xc) - 1 == cnt else "center"
        y = ys[0] - 1. + 0.5 * (cnt % 2)
        if len(xc) - 3 == cnt:
            y -= 0.1
        if len(xc) - 1 == cnt:
            y += 0.1
        ax.text(x, y, label, ha="center", va=va)

def grid3(ax, glsize, ys):
    # dxf
    xf = gen_xf(glsize)
    labels = (
            "DXF(1)",
            "DXF(2)",
            "DXF(3)",
            "dummy",
            "DXF(i-1)",
            "DXF(i)",
            "DXF(i+1)",
            "dummy",
            "DXF(size-2)",
            "DXF(size-1)",
            "DXF(size)",
    )
    for cnt, (x0, x1, label) in enumerate(zip(xf[:-1], xf[1:], labels)):
        if "dummy" == label:
            continue
        y = 0.5 * ys[0] + 0.5 * ys[1]
        ax.arrow(x0, y, x1-x0, 0, width=0.04, length_includes_head=True, color="#FF0000")
        ax.arrow(x1, y, x0-x1, 0, width=0.04, length_includes_head=True, color="#FF0000")
        x = 0.5 * x0 + 0.5 * x1
        y = ys[0] - 1.0 + 0.5 * (cnt % 2)
        ax.text(x, y, label, ha="center", va="center")

def grid4(ax, glsize, ys):
    # dxc
    xc = gen_xc(glsize)
    ax.scatter(xc, np.full(xc.shape, 0.5 * ys[0] + 0.5 * ys[1]), color="#000000")
    labels = (
            "DXC(1)",
            "DXC(2)",
            "DXC(3)",
            "dummy",
            "DXC(i)",
            "DXC(i+1)",
            "dummy",
            "DXC(size-1)",
            "DXC(size)",
            "DXC(size+1)",
    )
    for cnt, (x0, x1, label) in enumerate(zip(xc[:-1], xc[1:], labels)):
        if "dummy" == label:
            continue
        y = 0.5 * ys[0] + 0.5 * ys[1]
        ax.arrow(x0, y, x1-x0, 0, width=0.04, length_includes_head=True, color="#FF0000")
        ax.arrow(x1, y, x0-x1, 0, width=0.04, length_includes_head=True, color="#FF0000")
        x = 0.5 * x0 + 0.5 * x1
        y = ys[0] - 1. + 0.5 * (cnt % 2)
        if len(labels) - 3 == cnt:
            x -= 0.1
        if len(labels) - 1 == cnt:
            x += 0.1
        ax.text(x, y, label, ha="center", va="center")

def main(fname):
    glsize = 11
    ys = [0., 1.]
    fig = pyplot.figure(figsize=(glsize+2, np.max(ys)-np.min(ys)+2.))
    ax = pyplot.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    general(ax, glsize, ys)
    if "grid1" == fname:
        grid1(ax, glsize, ys)
    elif "grid2" == fname:
        grid2(ax, glsize, ys)
    elif "grid3" == fname:
        grid3(ax, glsize, ys)
    elif "grid4" == fname:
        grid4(ax, glsize, ys)
    else:
        msg = f"unknown file name: {fname}"
        raise RuntimeError(msg)
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/{fname}.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

if __name__ == "__main__":
    matplotlib.rcParams["lines.linewidth"] = 3
    matplotlib.rcParams["font.size"] = 20
    main("grid1")
    main("grid2")
    main("grid3")
    main("grid4")


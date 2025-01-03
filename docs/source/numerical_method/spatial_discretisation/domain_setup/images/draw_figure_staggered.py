import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import patches

is_periodic_x = False

isize = "isize"
jsize = "jsize"

def staggered1():
    matplotlib.rcParams["font.size"] = 20
    fig = pyplot.figure(figsize=(8, 4))
    ax = pyplot.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.text(1., 1.85, "Math", ha="center", va="center")
    ax.text(3., 1.85, "Code", ha="center", va="center")
    offsets = (
            {"x": 0.5, "y": 0.4},
            {"x": 2.5, "y": 0.4},
    )
    for cnt, offset in enumerate(offsets):
        xs = [offset["x"]-0.125, offset["x"]+1.125]
        ys = [offset["y"]-0.125, offset["y"]+1.125]
        ax.vlines(xs, ys[0], ys[1], colors="#000000", linestyles="-")
        ax.hlines(ys, xs[0], xs[1], colors="#000000", linestyles="-")
        # p and T
        x = 0.5 * xs[0] + 0.5 * xs[1]
        y = 0.5 * ys[0] + 0.5 * ys[1]
        ax.scatter([x], [y], c="#000000", s=10)
        if 0 == cnt:
            ax.text(x, y+0.125, "$T_{i,j}$", ha="center", va="center")
            ax.text(x, y-0.125, "$p_{i,j}$", ha="center", va="center")
        else:
            ax.text(x, y+0.125, "T(i,j)", ha="center", va="center")
            ax.text(x, y-0.125, "p(i,j)", ha="center", va="center")
        # arrows
        arrow_length = 0.5
        # x
        y = 0.5 * ys[0] + 0.5 * ys[1]
        ax.arrow(xs[0]-0.5*arrow_length, y, arrow_length, 0., width=0.07, length_includes_head=True, edgecolor="#FF0000", facecolor="none", linewidth=2., zorder=20.)
        ax.arrow(xs[1]-0.5*arrow_length, y, arrow_length, 0., width=0.07, length_includes_head=True, edgecolor="#FF0000", facecolor="none", linewidth=2., zorder=20.)
        for bgc, zorder in zip(["#FFFFFF", "none"], [10., 30.]):
            if 0 == cnt:
                ax.text(xs[0], y, "$u_{x,i-\\frac{1}{2},j}$", zorder=zorder, ha="center", va="center", backgroundcolor=bgc)
                ax.text(xs[1], y, "$u_{x,i+\\frac{1}{2},j}$", zorder=zorder, ha="center", va="center", backgroundcolor=bgc)
            else:
                ax.text(xs[0], y, "UX(i,j)",   zorder=zorder, ha="center", va="center", backgroundcolor=bgc)
                ax.text(xs[1], y, "UX(i+1,j)", zorder=zorder, ha="center", va="center", backgroundcolor=bgc)
        # y
        x = 0.5 * xs[0] + 0.5 * xs[1]
        ax.arrow(x, ys[0]-0.5*arrow_length, 0., arrow_length, width=0.07, length_includes_head=True, edgecolor="#0000FF", facecolor="none", linewidth=2., zorder=20.)
        ax.arrow(x, ys[1]-0.5*arrow_length, 0., arrow_length, width=0.07, length_includes_head=True, edgecolor="#0000FF", facecolor="none", linewidth=2., zorder=20.)
        if 0 == cnt:
            ax.text(x, ys[0], "$u_{y,i,j-\\frac{1}{2}}$", zorder=10., ha="center", va="center", backgroundcolor="#FFFFFF")
            ax.text(x, ys[1], "$u_{y,i,j+\\frac{1}{2}}$", zorder=10., ha="center", va="center", backgroundcolor="#FFFFFF")
        else:
            ax.text(x, ys[0], "UY(i,j)",   zorder=10., ha="center", va="center", backgroundcolor="#FFFFFF")
            ax.text(x, ys[1], "UY(i,j+1)", zorder=10., ha="center", va="center", backgroundcolor="#FFFFFF")
    keywords = {
            "aspect": "equal",
            "xlabel": "",
            "ylabel": "",
            "xlim": [0., 4.],
            "ylim": [0., 2.],
            "xticks": [],
            "yticks": [],
    }
    ax.set(**keywords)
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/staggered1.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

def gen_ax(sizes):
    matplotlib.rcParams["font.size"] = 15
    fig = pyplot.figure(figsize=(8, 8))
    ax = pyplot.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    keywords = {
            "aspect": "equal",
            "xlabel": "",
            "ylabel": "",
            "xlim": [-3, 11.],
            "ylim": [-3, 11.],
            "xticks": [],
            "yticks": [],
    }
    # water
    ax.add_patch(patches.Rectangle((0., 0.), sizes[0], sizes[1], color="#00FFFF"))
    if not is_periodic_x:
        # walls
        ax.add_patch(patches.Rectangle((-1., -1.), 1., 2.+sizes[1], color="#888888"))
        ax.add_patch(patches.Rectangle((sizes[0], -1.), 1., 2.+sizes[1], color="#888888"))
    # grid lines
    offsets = [-1., 3., 7.]
    for yoffset in offsets:
        ys = np.array([0., 1., 2., 3.])
        ys += yoffset
        for xoffset in offsets:
            xs = np.array([0., 1., 2., 3.])
            xs += xoffset
            ax.vlines(xs, ys[0], ys[-1], colors="#000000", linestyles="-")
            ax.hlines(ys, xs[0], xs[-1], colors="#000000", linestyles="-")
    ax.set(**keywords)
    return ax

def get_xparams(is_face, size):
    if is_face:
        if is_periodic_x:
            xs = np.arange(-1.0, size+1.0, 1)
            xlabels = (
                    "0",
                    "1",
                    "2",
                    "3",
                    "i-1",
                    "i",
                    "i+1",
                    "i+2",
                    f"{isize}-1",
                    f"{isize}",
                    f"{isize}+1",
            )
        else:
            xs = np.arange(0.0, size+1.0, 1)
            xlabels = (
                    "1",
                    "2",
                    "3",
                    "i-1",
                    "i",
                    "i+1",
                    "i+2",
                    f"{isize}-1",
                    f"{isize}",
                    f"{isize}+1",
            )
    else:
        if is_periodic_x:
            xs = np.arange(-1, size+1, 1) + 0.5
        else:
            xs = np.arange(0, size, 1) + 0.5
            xs = np.append(0., xs)
            xs = np.append(xs, size)
        xlabels = (
                "0",
                "1",
                "2",
                "dummy",
                "i-1",
                "i",
                "i+1",
                "dummy",
                f"{isize}-1",
                f"{isize}",
                f"{isize}+1",
        )
    return xs, xlabels

def get_yparams(is_face, size):
    if is_face:
        ys = np.arange(0, size+3, 1) - 1.
        ylabels = (
                "0",
                "1",
                "2",
                "dummy",
                "j-1",
                "j",
                "j+1",
                "dummy",
                f"{jsize}-1",
                f"{jsize}",
                f"{jsize}+1",
        )
    else:
        ys = np.arange(0, size+3, 1) - 0.5
        ylabels = (
                "0",
                "1",
                "2",
                "dummy",
                "j-1",
                "j",
                "j+1",
                "dummy",
                f"{jsize}-1",
                f"{jsize}",
                f"{jsize}+1",
        )
    return ys, ylabels

def attach_labels(ax, xs, ys, xlabels, ylabels):
    for cnt, (x, xlabel) in enumerate(zip(xs, xlabels)):
        if "dummy" == xlabel:
            continue
        y = -2. + 0.5 * (cnt % 2)
        ax.text(x, y, xlabel, ha="center", va="center")
    for cnt, (y, ylabel) in enumerate(zip(ys, ylabels)):
        if "dummy" == ylabel:
            continue
        x = -2.
        ax.text(x, y, ylabel, ha="center", va="center")

def get_grid(xs, ys, xlabels, ylabels):
    arrow_length = 0.5
    scale = 0.7
    qxs = list()
    qys = list()
    for cnt, (y, ylabel) in enumerate(zip(ys, ylabels)):
        for cnt, (x, xlabel) in enumerate(zip(xs, xlabels)):
            if "dummy" == ylabel or "dummy" == xlabel:
                continue
            qxs.append(x)
            qys.append(y)
    return arrow_length, scale, np.array(qxs), np.array(qys)

def staggered2(sizes):
    ax = gen_ax(sizes)
    xs, xlabels = get_xparams(is_face=True,  size=sizes[0])
    ys, ylabels = get_yparams(is_face=False, size=sizes[1])
    attach_labels(ax, xs, ys, xlabels, ylabels)
    length, scale, qxs, qys = get_grid(xs, ys, xlabels, ylabels)
    qxs -= 0.5 * length
    qus = np.full(qxs.shape, length)
    qvs = np.full(qxs.shape, 0.0)
    ax.quiver(qxs, qys, qus, qvs, scale=0.7, scale_units="x", color="#FF0000", zorder=10.)
    if attach_cr:
        ax.text(0.5 * sizes[0], 0.5 * sizes[1], **cr_config)
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/staggered2.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

def staggered3(sizes):
    ax = gen_ax(sizes)
    xs, xlabels = get_xparams(is_face=False, size=sizes[0])
    ys, ylabels = get_yparams(is_face=True,  size=sizes[1])
    attach_labels(ax, xs, ys, xlabels, ylabels)
    length, scale, qxs, qys = get_grid(xs, ys, xlabels, ylabels)
    qys -= 0.5 * length
    qus = np.full(qxs.shape, 0.0)
    qvs = np.full(qxs.shape, length)
    ax.quiver(qxs, qys, qus, qvs, scale=scale, scale_units="x", color="#0000FF", zorder=10.)
    if attach_cr:
        ax.text(0.5 * sizes[0], 0.5 * sizes[1], **cr_config)
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/staggered3.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

def staggered4(sizes):
    ax = gen_ax(sizes)
    xs, xlabels = get_xparams(is_face=False, size=sizes[0])
    ys, ylabels = get_yparams(is_face=False, size=sizes[1])
    attach_labels(ax, xs, ys, xlabels, ylabels)
    _, _, xs, ys = get_grid(xs, ys, xlabels, ylabels)
    ax.scatter(xs, ys, color="#000000", s=20.)
    if attach_cr:
        ax.text(0.5 * sizes[0], 0.5 * sizes[1], **cr_config)
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/staggered4.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

if __name__ == "__main__":
    attach_cr = False
    cr_config = {
        "s": "© 2024, Naoki Hori",
        "ha": "center",
        "va": "center",
        "size": "xx-large",
        "rotation": 44.78236,
        "color": "#aaaaaa",
    }
    matplotlib.rcParams["lines.linewidth"] = 3
    matplotlib.rcParams["axes.linewidth"] = 5
    staggered1()
    sizes = (9, 9)
    staggered2(sizes)
    staggered3(sizes)
    staggered4(sizes)


import os
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def main():
    lx = 3.
    ly = 2.
    fig = pyplot.figure(figsize=(lx, ly))
    ax = pyplot.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.plot([0.5, 2.5], [0.5, 0.5], color="#FF0000", linewidth=5)
    ax.plot([0.5, 2.5], [1.5, 1.5], color="#0000FF", linewidth=5)
    ax.text(1.5, 0.25, "$T_H$", color="#FF0000", ha="center", va="center")
    ax.text(1.5, 1.75, "$T_L$", color="#0000FF", ha="center", va="center")
    # horizontal length
    ax.text(1.0, 0.8, "$l_y$", ha="center", va="center")
    for case in [True, False]:
        x = 0.5 if case else 2.5
        y = 0.6
        dx = +2.0 if case else -2.0
        dy = 0.
        ax.arrow(x, y, dx, dy, width=0.03, length_includes_head=True, head_width=0.15, head_length=0.25, shape="right", edgecolor="#000000", facecolor="#000000")
    # vertical length
    ax.text(2.5, 1.0, "$l_x$", ha="center", va="center")
    for case in [True, False]:
        x = 2.7
        y = 0.5 if case else 1.5
        dx = 0.
        dy = +1.0 if case else -1.0
        ax.arrow(x, y, dx, dy, width=0.03, length_includes_head=True, head_width=0.15, head_length=0.25, shape="right", edgecolor="#000000", facecolor="#000000")
    # gravity
    ax.text(0.5, 1.0, "$g$", ha="center", va="center")
    x = 0.3
    y = 1.25
    dx = 0.
    dy = -0.5
    ax.arrow(x, y, dx, dy, width=0.03, length_includes_head=True, head_width=0.15, head_length=0.25, shape="full", edgecolor="#000000", facecolor="#000000")
    keywords = {
            "aspect": "equal",
            "xlabel": "",
            "ylabel": "",
            "xlim": [0, lx],
            "ylim": [0, ly],
            "xticks": [],
            "yticks": [],
    }
    ax.set(**keywords)
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/schematic.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

if __name__ == "__main__":
    matplotlib.rcParams["font.size"] = 20
    main()

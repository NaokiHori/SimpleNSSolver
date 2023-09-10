import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import patches

def general(ax, glisize, gljsize):
    # vertical lines
    xs = np.arange(1, glisize, 1)
    ys = [0., gljsize]
    ax.vlines(xs, ys[0], ys[1], colors="#000000", linestyles="--")
    # horizontal lines
    xs = [0., glisize]
    ys = np.arange(1, gljsize, 1)
    ax.hlines(ys, xs[0], xs[1], colors="#000000", linestyles="--")
    keywords = {
            "aspect": "equal",
            "xlabel": "",
            "ylabel": "",
            "xlim": [0, glisize],
            "ylim": [0, gljsize],
    }
    ax.set(**keywords)

def main():
    glisize =  7
    gljsize = 10
    myisize = glisize
    myjsize = [3, 3, 4]
    offsets = [0, 3, 6]
    fig = pyplot.figure(figsize=(7.0, 3.5))
    ax121 = fig.add_subplot(121)
    ax122 = fig.add_subplot(122)
    general(ax121, glisize, gljsize)
    general(ax122, glisize, gljsize)
    ax121.set_title("global domain")
    ax122.set_title("local domain")
    ax121.set_xticks(ticks=[0.5, glisize-0.5], labels=["1", "glsizes[0]"])
    ax121.set_yticks(ticks=[0.5, gljsize-0.5], labels=["1", "glsizes[1]"])
    ax122.set_xticks(ticks=[0.5, myisize-0.5], labels=["1", "mysizes[0]"])
    ax122.set_yticks(
            ticks=[
                offsets[0]+0.5, offsets[0]+myjsize[0]-0.5,
                offsets[1]+0.5, offsets[1]+myjsize[1]-0.5,
                offsets[2]+0.5, offsets[2]+myjsize[2]-0.5,
            ],
            labels=[
                "1", "mysizes[1]",
                "1", "mysizes[1]",
                "1", "mysizes[1]",
            ]
    )
    ax121.add_patch(patches.Rectangle(xy=(0, 0), width=glisize, height=gljsize, facecolor="#888888", edgecolor="#888888", zorder=-10, alpha=0.5))
    ax122.add_patch(patches.Rectangle(xy=(0, offsets[0]), width=myisize, height=myjsize[0], facecolor="#FF0000", edgecolor="#FF0000", zorder=-10, alpha=0.5))
    ax122.add_patch(patches.Rectangle(xy=(0, offsets[1]), width=myisize, height=myjsize[1], facecolor="#0000FF", edgecolor="#0000FF", zorder=-10, alpha=0.5))
    ax122.add_patch(patches.Rectangle(xy=(0, offsets[2]), width=myisize, height=myjsize[2], facecolor="#33AA00", edgecolor="#33AA00", zorder=-10, alpha=0.5))
    abspath = os.path.realpath(__file__)
    dirname = os.path.dirname(abspath)
    filename = f"{dirname}/domain.png"
    print(filename)
    pyplot.savefig(filename)
    pyplot.close()

if __name__ == "__main__":
    matplotlib.rcParams["lines.linewidth"] = 1
    matplotlib.rcParams["axes.linewidth"] = 5
    main()


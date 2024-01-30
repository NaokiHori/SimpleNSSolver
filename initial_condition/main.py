import os
import sys
import numpy as np


def init_time(dest):
    # iterator and time
    step = np.array(0, dtype=np.uint64)
    time = np.array(0, dtype=np.float64)
    np.save(f"{dest}/step.npy", step)
    np.save(f"{dest}/time.npy", time)
    return


def init_domain(lengths, glsizes, uniformx, dest):
    # NOTE: cell face has +1 elements
    if uniformx:
        # uniform grid in x,
        #   which is advantageous to solve Poisson equation more efficiently
        xf = np.linspace(0., lengths[0], glsizes[0] + 1, endpoint=True)
    else:
        # stretched grid, clipped Chebyshev just as an example
        # number of grid points to be clipped at the edges
        nclip = 3
        # generate equidistant sequence
        xf = np.arange(0, glsizes[0] + 1, 1)
        # gather close to the boundaries
        xf = np.cos(np.pi * (xf + 1. * nclip) / (glsizes[0] + 2. * nclip))
        # make the descending order ascending
        xf *= -1.
        # normalse to force it changing from 0 to lx
        xf = lengths[0] * (xf - np.min(xf)) / (np.max(xf) - np.min(xf))
    # cell centers are located at the center
    #   of the two neighbouring cell faces,
    #   which are appended by the boundaries
    xc = 0.
    xc = np.append(xc, 0.5 * xf[:-1] + 0.5 * xf[1:])
    xc = np.append(xc, lengths[0])
    np.save(f"{dest}/xf.npy", np.array(xf, dtype=np.float64))
    np.save(f"{dest}/xc.npy", np.array(xc, dtype=np.float64))
    np.save(f"{dest}/glsizes.npy", np.array(glsizes, dtype=np.uint64))
    np.save(f"{dest}/lengths.npy", np.array(lengths, dtype=np.float64))
    return xf, xc


def init_fluid(lengths, glsizes, xf, xc, dest):
    shape = (glsizes[1], glsizes[0] + 1)
    ux = np.zeros(shape, dtype=np.float64)
    shape = (glsizes[1], glsizes[0] + 2)
    uy = np.zeros(shape, dtype=np.float64)
    shape = (glsizes[1], glsizes[0] + 2)
    p = np.zeros(shape, dtype=np.float64)
    shape = (glsizes[1], glsizes[0] + 2)
    t = -0.5 + np.random.random_sample(shape)
    np.save(f"{dest}/ux.npy", ux)
    np.save(f"{dest}/uy.npy", uy)
    np.save(f"{dest}/p.npy", p)
    np.save(f"{dest}/t.npy", t)


def main():
    lengths = list()
    lengths.append(float(os.environ["lx"]))
    lengths.append(float(os.environ["ly"]))
    glsizes = list()
    glsizes.append(int(os.environ["glisize"]))
    glsizes.append(int(os.environ["gljsize"]))
    uniformx = os.environ["uniformx"]
    if "True" == uniformx:
        uniformx = True
    elif "true" == uniformx:
        uniformx = True
    else:
        uniformx = False
    dest = sys.argv[1]
    # sanitise
    ndims = len(lengths)
    assert 2 == ndims
    # init and save
    init_time(dest)
    xf, xc = init_domain(lengths, glsizes, uniformx, dest)
    init_fluid(lengths, glsizes, xf, xc, dest)


main()

import sys
import numpy as np
import pyvista
from pyvista import themes

def load(dname):
    time = np.load(f"{dname}/time.npy")
    lengths = np.load(f"{dname}/lengths.npy")
    glsizes = np.load(f"{dname}/glsizes.npy")
    xc = np.load(f"{dname}/xc.npy")
    lx = lengths[0]
    ly = lengths[1]
    lz = lengths[2]
    ny = glsizes[1]
    nz = glsizes[2]
    dy = ly / ny
    dz = lz / nz
    yc = np.linspace(0.5 * dy, ly - 0.5 * dy, ny)
    zc = np.linspace(0.5 * dz, lz - 0.5 * dz, nz)
    t  = np.load(f"{dname}/t.npy")
    return time, lx, ly, lz, xc, yc, zc, t

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    idname = argv[1]
    ofname = argv[2]
    time, lx, ly, lz, x, y, z, t = load(idname)
    grid = pyvista.RectilinearGrid(x, y, z)
    t = np.ravel(t)
    my_theme = themes.Theme()
    my_theme.lighting = True
    my_theme.show_edges = False
    my_theme.background = "#000000"
    my_theme.window_size = [1200, 800]
    thresholds = [-0.1, +0.1]
    colours = ["#0000FF", "#FF0000"]
    isosurfaces = 1
    plotter = pyvista.Plotter(theme=my_theme, off_screen=True)
    for threshold, colour in zip(thresholds, colours):
        contour = grid.contour(isosurfaces=isosurfaces, scalars=t, rng=[threshold, threshold])
        contour = contour.smooth_taubin(n_iter=30, pass_band=0.05)
        plotter.add_mesh(contour, color=colour, opacity=1.)
    plotter.camera_position = 'xy'
    plotter.camera.roll = 90
    plotter.camera.azimuth = 30
    plotter.camera.elevation = 10.
    plotter.camera.zoom(1.2)
    plotter.show(screenshot=ofname)
    plotter.close()


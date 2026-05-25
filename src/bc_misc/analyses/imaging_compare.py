"""
Compare imaging with xraysky vs bc-imager
"""

import numpy as np
import matplotlib.pyplot as plt
import xraysky
from astropy.io import fits
import bc_imager
from pathlib import Path

_x = plt.ion()


def py_imager(resolution=2):
    inst = xraysky.BlackCAT(detgap=1788e-6)
    imager = xraysky.FFTImager(inst, resolution=resolution)
    return imager


def events_to_pyimager_pixels(events, instrument):
    """Convert BlackCAT Level1 events to pixel coordinatesxd,yd

    Args:
        events (_type_): _description_
        instrument (_type_): _description_
    """
    pixev_dtype = np.dtype([("xd", np.float32), ("yd", np.float32)])
    result = np.zeros(len(events), dtype=pixev_dtype)
    for d, low, pixelsize, destcol, sourcecol in zip(
        (0, 1),
        instrument.detenvelope()[0],
        instrument.detpixelsize,
        ("xd", "yd"),
        ("detx", "dety"),
    ):
        result[destcol] = (events[sourcecol] - low) * 1 / pixelsize
    return result


def pyimage(events, pyimager):
    ev = events_to_pyimager_pixels(events, pyimager.instrument)
    im = pyimager.image(ev)
    return im


def main():
    topdir = Path(__file__).parent.parent.parent.parent
    fev1 = topdir / "data/image_test/simulated/simulated_events_L1.gz"
    imager = py_imager()

    d, h = fits.getdata(fev1, header=True)
    trange = np.array([40, 70]) + h["tstart"]
    dpeak = d[(trange[0] < d["time"]) & (d["time"] < trange[1])]
    im = pyimage(dpeak, imager)
    pass


if __name__ == "__main__":
    main()

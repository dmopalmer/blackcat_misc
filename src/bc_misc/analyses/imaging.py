"""
imaging.py
"""
# %%
from tools import *
from xraysky import BC_Imaging


# From Tim
BASE_PATH = top_level / "CHECKOUT_26_05_22"
LEVEL1_PATH = BASE_PATH / "level1_tmp"
PHOTON_PATH = LEVEL1_PATH / "ph_2605221058.gz"
ORIENT_PATH = LEVEL1_PATH / "or_2605221058.gz"

orientations = fits.getdata(ORIENT_PATH, 1)
orientation_start_time = np.min(orientations["TIME"])

stable_time = orientation_start_time + 150

stable_orientations = orientations[orientations["TIME"] >= stable_time]
stable_pointings = stable_orientations["POINTING_NOMINAL"]
stable_ra = np.median(stable_pointings[:,0])
stable_dec = np.median(stable_pointings[:,1])
stable_roll = np.median(stable_pointings[:,2])

with fits.open(PHOTON_PATH) as event_hdul:
    event_hdu = event_hdul[1]
    events = event_hdu.data
    event_hdu.data = events[events["TIME"] >= stable_time]

    imager = BC_Imaging(
        radecroll_inst=[
            stable_ra,
            stable_dec,
            -stable_roll,
        ]
    )
    im = imager.evtlist2image(
        event_hdu,
        True,
        Path("/tmp/ph_2605221058_image_neg_roll.fits.gz",
    )

# %%

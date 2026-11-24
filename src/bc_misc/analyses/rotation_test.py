"""
How do you get the rotation etc of the xraysky imaging to match that you get from WCS

"""


# %%

# %%
import matplotlib

matplotlib.use("MacOSX")  # Use 'osx' in some older IPython environments

from tools import *
from xraysky import BC_Imager
from astropy.wcs import WCS
import astropy.units as u
from astropy.stats import biweight_location
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import re
import warnings
from astropy.wcs import FITSFixedWarning

# Suppress only Astropy FITS WCS warnings
warnings.filterwarnings("ignore", category=FITSFixedWarning)
try:
    sourcedir = Path(__file__).parent
except NameError:
    sourcedir = Path(os.getcwd())

outdir = Path("/tmp/images")
outdir.mkdir(exist_ok=True)


imager = BC_Imager(
    resolution=2,
    overwrite=True,
    balance=calibration["balance_boxes"],
    hide_frame=True,
)


if True:
    ph_all = sorted(
        set(
            [
                line.split()[0].split("_")[0]
                for line in sourcedir.joinpath("found_scox1.txt").open("rt").readlines()
                if "ScoX1+" in line
            ]
        )
    )
else:
    ph_all = sorted(
        top_dir.rglob("payload_checkout/CHECKOUT_26_05_*/level1_tmp/ph_2*")
    ) + sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_06_*/level1_tmp/ph_2*"))

# %%
f = Path(ph_all[-1])

obs = obs_photons(f)

counts = obs["stable_photons"]
radecroll = obs["radecroll"]
# counts,radecroll = calibrate_data(obs['stable_photons'], obs['radecroll'])
obsname = f"{obs['obsname']}"
erange = np.percentile(counts["ENERGY"], (5, 95))
imager.set_radecroll(radecroll)
imhdu = imager.evtlist2image(
    counts=counts,
    header=obs["header"],
    add_wcs=True,
    outfile=None,
)

wcs = WCS(imhdu.header)
height, width = imhdu.data.shape
scoij = wcs.world_to_pixel(scox1_coords)
fov_scox1 = 0 < scoij[0] < width and 0 < scoij[1] < height
# print("finding peaks")
peaks = imager.imager.findpeaks(imhdu.data, minsigma=5)
# print(f"Found {len(peaks)} peaks")
# Can get MANY peaks with local significance due to partial detector population
peaks = peaks[peaks["signif_global"] > 2]
# ax.plot(width-scoij[0], height-scoij[1], 'v', markerfacecolor='none', markeredgecolor=('k', 0.75))
# ax.plot(width/2, height/2, '+', markerfacecolor='none', markeredgecolor=('k', 0.75))
dt = np.ptp(counts["time"])
dets = sorted(set(counts["detid"]))
rate = len(counts) / max(dt, 0.1)

peak = peaks[0]

pixoff = peak["ijpeak"][::-1] - scoij
# alignment is handled by 'calibrate_data()' # - alignment_corr/impixsize
radec = wcs.pixel_to_world([peak["ijpeak"][1]], [peak["ijpeak"][0]])
ra, dec = [float(v.deg.item()) for v in (radec.ra, radec.dec)]
print(
    f"{obsname} ji = {peak['ijpeak'][1]:5.1f},{peak['ijpeak'][0]:5.1f} = ({ra:8.3f}, {dec:+8.3f}) ra,dec = {peak['radecpeak'][0]:6.2f},{peak['radecpeak'][1]:6.2f}  {peak['signif_local']:5.1f} σ for {dt:.0f} s at {rate:.0f} cps {len(dets):1} dets"
)
print(
    f"{obsname} ji = {peak['ijpeak'][1]:5.1f},{peak['ijpeak'][0]:5.1f} = ScoX1+ {pixoff[0]:5.1f}, {pixoff[1]:5.1f} ra,dec = {peak['radecpeak'][0]:6.2f},{peak['radecpeak'][1]:6.2f}  {peak['signif_local']:5.1f} σ for {dt:.0f} s at {rate:.0f} cps {len(dets):1} dets"
)

# Position based on xraysky conversion
nomij = wcs.world_to_pixel(
    SkyCoord(
        ra=peak["radecpeak"][0] * u.degree,
        dec=peak["radecpeak"][1] * u.degree,
        frame="icrs",
    )
)
print(f"Round trip ij = {nomij[1].item():5.1f},{nomij[0].item():5.1f},")

# %%

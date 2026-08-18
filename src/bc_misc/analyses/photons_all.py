"""
DO lots of plots of observations
"""
#%%
import matplotlib
matplotlib.use('MacOSX')  # Use 'osx' in some older IPython environments

from tools import *
from xraysky import BC_Imager
from astropy.wcs import WCS
import astropy.units as u
from astropy.stats import biweight_location
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import re
import warnings
from astropy.wcs import FITSFixedWarning
import blackcat

# Suppress only Astropy FITS WCS warnings
warnings.filterwarnings('ignore', category=FITSFixedWarning)
try:
    sourcedir = Path(__file__).parent
except NameError:
    sourcedir = Path(os.getcwd())
    
outdir = Path("/tmp/photon_stacks")
outdir.mkdir(exist_ok=True)


imager = BC_Imager(
    resolution=2,
    overwrite = True,
    balance=calibration["balance_boxes"],
    hide_frame=True,
)


if False:
    ph_all = sorted(set([line.split()[0].split("_")[0] for line in sourcedir.joinpath('found_scox1.txt').open('rt').readlines() if 'ScoX1+' in line]))
else:
    ph_all = sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_08_10/level1_tmp/ph_2*"))[::-1]
  
plot=True
fitsfiles=False

#%%

for f in ph_all:
    f=Path(f)
    try:
        obs = obs_data(f)
        fig, d, h = plot_photons_isl4(obs['obsname'], obs['l1_file'])
        fig.savefig(outdir.joinpath(obs['obsname']+"_phstack.png"))
    except Exception as e:
        print(f"Failed on {obs['obsname']} {e}")
    plt.close('all')
# %%
# Diagonals
obs = obs_data('2606212204')
ph_0 = fits.getdata(obs['l0_file'])
ph_1 = blackcat.convert_events_level0_to_level1(obs['l0_file'], obs['obsname'])
pass

# %%

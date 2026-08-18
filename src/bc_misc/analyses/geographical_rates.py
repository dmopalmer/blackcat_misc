"""
Show photons and rates
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

# Suppress only Astropy FITS WCS warnings
warnings.filterwarnings('ignore', category=FITSFixedWarning)
try:
    sourcedir = Path(__file__).parent
except NameError:
    sourcedir = Path(os.getcwd())
    
outdir = Path("/tmp/geo_images")
outdir.mkdir(exist_ok=True)


# %%


for daydir in sorted(top_dir.rglob(f"payload_checkout/CHECKOUT_*")):
    plt.close('all')
    try:
        datasource = bcd.DataSource(daydir.rglob('**/HOUSEKEEPING*'))
        daystring = f'{datasource.data('lapped')[0][-1]['t'].item():%Y-%m-%d}'
        fig,ax = plot_geographic_rates(datasource)
        fig.suptitle(f"Good count rates for {datasource.data('lapped')[0][-1]['t'].item():%Y-%m-%d}")
        fig.savefig(outdir.joinpath(f"geo_{daystring}.png"))
    except:
        pass
# %%

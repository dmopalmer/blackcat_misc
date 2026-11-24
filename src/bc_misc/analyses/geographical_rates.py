"""
Show photons and rates
"""

# %%
import matplotlib

matplotlib.use("MacOSX")  # Use 'osx' in some older IPython environments

from bc_misc.analyses.tools import plot_pointing_fovs, plot_rates_vs_particles
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

outdir = Path("/tmp/geo_images")
outdir.mkdir(exist_ok=True)


# %%


for daydir in sorted(top_dir.rglob(f"payload_checkout/CHECKOUT_*"))[::-1]:
    plt.close("all")
    try:
        datasource = bcd.DataSource(daydir.rglob("**/HOUSEKEEPING*"))
        daystring = f"{datasource.data('lapped')[0][-1]['t'].item():%Y-%m-%d}"
        fig, ax = plot_geomagnetic_rates(datasource)
        fig.suptitle(
            f"Good count rates for {datasource.data('lapped')[0][-1]['t'].item():%Y-%m-%d}"
        )
        fig.savefig(outdir.joinpath(f"geo_{daystring}.png"))
        fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12, 12))
        for ax, (particle, energy) in zip(
            np.ravel(axes), [("e", 50), ("e", 1000), ("p", 10e3), ("p", 100e3)]
        ):
            fig, ax = plot_rates_vs_particles(
                datasource, particle=particle, e_keV_low=energy, fig=fig, ax=ax
            )
            ax.set(label=f"{particle} E > {energy / 1e3} MeV", xlim=[0.9, 1.1e3])
        fig.suptitle(f"{daystring} Counts vs predicted particles")
        fig.savefig(outdir.joinpath(f"particle_{daystring}.png"))

    except Exception as e:
        print(daydir.name, e)
        break
        pass
# %%

daydir = sorted(top_dir.rglob(f"payload_checkout/CHECKOUT_26_10_??"))[-1]
datasource = bcd.DataSource(daydir.rglob("**/HOUSEKEEPING*"))
plot_pointing_fovs(datasource, maxpointings=10)
# %%

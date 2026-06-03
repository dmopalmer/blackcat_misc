"""
Image lots of data
"""
#%%
import matplotlib
matplotlib.use('MacOSX')  # Use 'osx' in some older IPython environments

from tools import *
from xraysky import BC_Imager
from astropy.wcs import WCS
import astropy.units as u

import warnings
from astropy.wcs import FITSFixedWarning

# Suppress only Astropy FITS WCS warnings
warnings.filterwarnings('ignore', category=FITSFixedWarning)

sourcedir = Path(__file__).parent
outdir = Path("/tmp/images")
outdir.mkdir(exist_ok=True)


imager = BC_Imager(
    resolution=2,
    overwrite = True,
    balance=balance_boxes
)


if True:
    ph_all = sorted(set([line.split()[0] for line in sourcedir.joinpath('found_scox1.txt').open('rt').readlines()]))
else:
    ph_all = sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_05_[23]*/level1_tmp/ph_2*"))
    
#%%
for f in ph_all:
    f=Path(f)
    try:
        obs = obs_photons(f)
        obsname = obs['obsname']
        imager.set_radecroll(obs['radecroll'])
        imhdu = imager.evtlist2image(
                counts=obs['stable_photons'], header=obs['header'],
                add_wcs=True,
                outfile=outdir.joinpath(f"{obsname}_pos_roll.fits.gz"),
            )
        wcs = WCS(imhdu.header)
        height,width = imhdu.data.shape
        fig = plt.figure(figsize=(15,9))
        ax = fig.add_subplot(1,1,1, projection=wcs)
        ax.coords['ra'].set_format_unit(u.degree)
        ax.imshow(imhdu.data, origin='lower')
        peaks = imager.imager.findpeaks(imhdu.data)
        scoij = wcs.world_to_pixel(scox1_coords)
        ax.plot(scoij[0], scoij[1], 's', markerfacecolor='none', markeredgecolor=('k', 0.75))
        ax.plot(width-scoij[0], height-scoij[1], 'v', markerfacecolor='none', markeredgecolor=('k', 0.75))
        ax.plot(width/2, height/2, '+', markerfacecolor='none', markeredgecolor=('k', 0.75))
        if len(peaks):
            for peak in peaks:
                print(f"{obsname} ij = {peak['jipeak'][1]:5.1f},{peak['jipeak'][0]:5.1f} ra,dec = {peak['radecpeak'][0]:6.2f},{peak['radecpeak'][1]:6.2f}  {peak['signif_local']:.1f} σ")
            ax.plot(peaks['jipeak'][:,1], peaks['jipeak'][:,0], 'o', markersize=50, markerfacecolor='none', markeredgecolor='k')
        ax.grid(color='black', linestyle='dotted', alpha=0.5)
        ax.set(title=f'{obsname} roll = {obs['radecroll'][2]:.2f}')
        fig.tight_layout()
        fig.savefig(outdir.joinpath(f"{obsname}_full.png"))
        plt.pause(2)
        # plt.close('all')
    except Exception as e:
        print(f'{f.name} : {e}')
        break
pass
# %%

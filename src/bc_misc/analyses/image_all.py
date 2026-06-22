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
    
outdir = Path("/tmp/images")
outdir.mkdir(exist_ok=True)


imager = BC_Imager(
    resolution=2,
    overwrite = True,
    balance=calibration["balance_boxes"],
    hide_frame=True,
)


if True:
    ph_all = sorted(set([line.split()[0].split("_")[0] for line in sourcedir.joinpath('found_scox1.txt').open('rt').readlines()]))
else:
    ph_all = (sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_05_*/level1_tmp/ph_2*")) +
        sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_06_*/level1_tmp/ph_2*")))
    

#%%
plot=True
fitsfiles=False
# Corrective offsets to add to detector detx, dety
detpixsize = 40e-6
impixsize = imager.resolution * detpixsize

# alignment_corr = calibration['alignment_corr']

pixoffs = []
for f in ph_all:
    f=Path(f)
    try:
        obs = obs_photons(f)
        
        if not 'stable_photons' in obs:
            continue
        # counts,radecroll = calibrate_data(obs['stable_photons'], obs['radecroll'])
        counts,radecroll = obs['stable_photons'], obs['radecroll']
        obsname = f"{obs['obsname']}_alldets"
        if len(counts) < 10:
            continue
        imager.set_radecroll(radecroll)
        imhdu = imager.evtlist2image(
                counts=counts, header=obs['header'],
                add_wcs=True,
                outfile=outdir.joinpath(f"{obsname}_pos_roll.fits.gz") if fitsfiles else None,
            )
        wcs = WCS(imhdu.header)
        height,width = imhdu.data.shape
        # print("finding peaks")
        peaks = imager.imager.findpeaks(imhdu.data, minsigma=5)
        # print(f"Found {len(peaks)} peaks")
        # Can get MANY peaks with local significance due to partial detector population
        peaks = peaks[peaks['signif_global'] > 2]
        scoij = wcs.world_to_pixel(scox1_coords)
        # ax.plot(width-scoij[0], height-scoij[1], 'v', markerfacecolor='none', markeredgecolor=('k', 0.75))
        # ax.plot(width/2, height/2, '+', markerfacecolor='none', markeredgecolor=('k', 0.75))
        dt = np.ptp(counts['time'])
        dets = sorted(set(counts['detid']))
        rate = len(counts)/max(dt,0.1)
        
        if len(peaks):
            for peak in peaks:
                pixoff = peak['ijpeak'][::-1] - scoij
                # alignment is handled by 'calibrate_data()' # - alignment_corr/impixsize
                if np.hypot(*pixoff) > 50:
                    continue
                print(f"{obsname} ij = {peak['ijpeak'][1]:5.1f},{peak['ijpeak'][0]:5.1f} = ScoX1+ {pixoff[0]:5.1f}, {pixoff[1]:5.1f} ra,dec = {peak['radecpeak'][0]:6.2f},{peak['radecpeak'][1]:6.2f}  {peak['signif_local']:5.1f} σ for {dt:.0f} s at {rate:.0f} cps {len(dets):1} dets")
                pixoffs.append(np.concatenate((scoij,pixoff)))

        if plot:
            plt.close('all')
            fig = plt.figure(figsize=(15,9))
            ax = fig.add_subplot(1,1,1, projection=wcs)
            ax.imshow(imhdu.data, origin='lower')
            # ax.plot(scoij[0], scoij[1], 's', markerfacecolor='none', markeredgecolor=('k', 0.75))
            ax.coords['ra'].set_format_unit(u.degree)
            ax.plot(peaks['ijpeak'][:,1], peaks['ijpeak'][:,0], 'o', markersize=15, markerfacecolor='none', markeredgecolor='k')
            ax.grid(color='black', linestyle='dotted', alpha=0.5)
            if len(peaks):
                ny,nx = imhdu.data.shape
                vfrac = 0.2
                hfrac = vfrac * ny/nx
                inax = ax.inset_axes(bounds = [vfrac/4, hfrac/4, hfrac, vfrac])
                inax.imshow(imhdu.data, origin='lower', interpolation='none')
                inax.axvline(scoij[0], color='#FF000080')
                inax.axhline(scoij[1], color='#FF000080')
                inax.set(xlim=[scoij[0]-20, scoij[0]+20], ylim=[scoij[1]-20, scoij[1]+20], xticks=[], yticks=[])
                ax.indicate_inset_zoom(inax, edgecolor="black")

            ax.set(title=f'OBSID 20{obs['obsname']} for {dt:.0f} s, {rate:.0f} cps, >5 keV')
            fig.tight_layout()
            fig.savefig(outdir.joinpath(f"{obsname}_full.png"))
            plt.pause(2)
    except Exception as e:
        print(f'{f.name} : {e}')

pixoffs = np.asarray(pixoffs)
#%%

# %%
# Significance plot

dur_kappa_ndets = []

for line in sourcedir.joinpath('found_scox1.txt').open('rt').readlines():
    try:
        spl = [v for v_ in re.split(r'[= ,]', line) if (v := v_.strip())]
        kappa, dur, rate, ndets = [float(spl[ind]) for ind in [11,14, 17, 19]]
        dur_kappa_ndets.append([dur, kappa, ndets])
    except :
        pass
    
dur_kappa_ndets = np.array(dur_kappa_ndets)

fig,ax = plt.subplots()
for n in range(4,0,-1):
    dur,kappa,nd = dur_kappa_ndets[dur_kappa_ndets[:,2] == n].T
    ax.plot(dur, kappa, '/.|vs'[n], label=f"{n} dets")
    ax.legend()
ax.set(title="Sco X-1 significance >5 keV vs duration", xlabel="Duration (s)", ylabel="Image peak significance (σ)")
fig.tight_layout()
fig.savefig("/tmp/significance.png")

fig,ax = plt.subplots()
for n in range(4,0,-1):
    dur,kappa,nd = dur_kappa_ndets[dur_kappa_ndets[:,2] == n].T
    ax.plot(dur, kappa/np.sqrt(dur), r'/.|vs'[n], label=f"{n} dets")
    ax.legend()
ax.set(title="Sco X-1 significance >5 keV per sqrt(duration)", xlabel="Duration (s)", ylabel=r"Significance rate ($\sigma s^{-{1/2}}$)")
fig.tight_layout()
fig.savefig("/tmp/significance_rate.png")



# %%

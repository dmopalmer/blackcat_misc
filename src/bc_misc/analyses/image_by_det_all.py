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
    resolution=1,
    overwrite = True,
    balance=calibration["balance_boxes"],
    hide_frame=True,
)


if True:
    ph_all = sorted(set([line.split()[0].split("_")[0] for line in sourcedir.joinpath('found_scox1.txt').open('rt').readlines()]))
else:
    ph_all = (sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_05_[23]*/level1_tmp/ph_2*")) +
        sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_06_*/level1_tmp/ph_2*")))
    
#%%
plot=True
# Corrective offsets to add to detector detx, dety
detpixsize = 40e-6
impixsize = imager.resolution * detpixsize

# The gross offset which is probably due to instrument alignment of ~0.47°
# comparable to the ~0.5° offset found between the two star camera nominals.
# alignment_corr = np.array([-0.00115047, -0.00039763])
# with 2026-06-07 code, resolution =1 measured alignment corection is [30.70622999 17.69624394] um
# And with resolution=4 the value is consistent
alignment_corr = np.array([ -1119.8e-6,   -379.9e-6])

# How much to add to each detx, dety
# Adding to detx moves image spot to the right, adding to dety moves image spot down.
corr_by_det =\
    np.array([
    [    32.9e-6,    189.0e-6],
    [     3.1e-6,   -167.7e-6],
    [  -232.2e-6,     96.2e-6],
    [   195.9e-6,   -117.6e-6],
    ])


# Error of peak location [i,j] - WCS(Sco X-1)
pixoffs = {d:[] for d in range(4)}

print(f"Assumed alignment correction = {alignment_corr * 1e6} microns = {np.rad2deg(alignment_corr / imager.instrument.flength)} degrees")
print(f"Per-detector correction = \n{corr_by_det * 1e6} microns")

for f in ph_all:
    f=Path(f)
    try:
        obs = obs_photons(f)
        if not 'stable_photons' in obs:
            continue
        for det,corr in zip(range(4), corr_by_det):
            obsname = f"{obs['obsname']}_d{det}"
            counts = obs['stable_photons']
            counts = counts[counts['DETID'] == det]
            if len(counts) < 10:
                continue

            counts['DETX'] += corr[0]
            counts['DETY'] += corr[1]
            imager.imager.balance[det] = [[counts['DETX'].min(), counts['DETY'].min()],[counts['DETX'].max(), counts['DETY'].max()]]
            imager.set_radecroll(obs['radecroll'])
            imhdu = imager.evtlist2image(
                    counts=counts, header=obs['header'],
                    add_wcs=True,
                    outfile=outdir.joinpath(f"{obsname}_pos_roll.fits.gz") if plot else None,
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
            if len(peaks):
                for peak in peaks:
                    pixoff = peak['jipeak'][::-1] - scoij - alignment_corr/impixsize
                    if np.hypot(*pixoff) > 50:
                        continue
                    print(f"{obsname} ij = {peak['jipeak'][1]:5.1f},{peak['jipeak'][0]:5.1f} = ScoX1+ {pixoff[0]:5.1f}, {pixoff[1]:5.1f} ra,dec = {peak['radecpeak'][0]:6.2f},{peak['radecpeak'][1]:6.2f}  {peak['signif_local']:.1f} σ")
                    pixoffs[det].append(pixoff)
            if plot:
                fig = plt.figure(figsize=(15,9))
                ax = fig.add_subplot(1,1,1, projection=wcs)
                ax.imshow(imhdu.data, origin='lower')
                ax.plot(scoij[0], scoij[1], 's', markerfacecolor='none', markeredgecolor=('k', 0.75))
                ax.coords['ra'].set_format_unit(u.degree)
                ax.plot(peaks['jipeak'][:,1], peaks['jipeak'][:,0], 'o', markersize=15, markerfacecolor='none', markeredgecolor='k')
                ax.grid(color='black', linestyle='dotted', alpha=0.5)
                ax.set(title=f'{obsname} roll = {obs['radecroll'][2]:.2f}')
                fig.tight_layout()
                fig.savefig(outdir.joinpath(f"{obsname}_full.png"))
                plt.pause(2)
        plt.close('all')
    except Exception as e:
        print(f'{f.name} : {e}')

# Per-detector measured correction
# If the image shows the peak as in +i direction from ScoX-1, (offs is positive), then add the corresponding meters to detX
# and subtract the corresponding meters from detY
imeas_im_offset_by_det_m = []
for i in range(4):
    pixoffs[i] = np.array(pixoffs[i])
    imeas_im_offset_by_det_m.append(biweight_location(pixoffs[i], axis=0) * impixsize) 
    
imeas_im_offset_by_det_m = np.array(imeas_im_offset_by_det_m)
meas_alignment_corr = np.mean(imeas_im_offset_by_det_m, axis=0)
imeas_im_offset_by_det_m -= meas_alignment_corr

new_alignment_corr = meas_alignment_corr + alignment_corr
new_corr_by_det = imeas_im_offset_by_det_m *[+1,-1] + corr_by_det

print(f"Measured alignment correction = {meas_alignment_corr * 1e6} um + previous correction = {alignment_corr * 1e6}um gives")
print(f"alignment_corr = np.array([{new_alignment_corr[0]*1e6:8.1f}e-6, {new_alignment_corr[1]*1e6:8.1f}e-6])")

print(f"TODO: convert this to a quaternion")
print(f"TODO: Include instrument twist")

print(f"Per-detector image offset =\n{imeas_im_offset_by_det_m}")
print(f"Added [+,-] to current correction = \n{corr_by_det}")
print(f"Gives a new correction of\ncorr_by_det =\\\n    np.array([")
for d in range(4):
    print(f"    [{new_corr_by_det[d,0]*1e6:8.1f}e-6, {new_corr_by_det[d,1]*1e6:8.1f}e-6],")
print("    ])")
fig,axes = plt.subplots(2,2, figsize=[10,10], sharex=True, sharey=True)
axes = axes.ravel()[[1,2,0,3]]
for det,ax in enumerate(axes):
    ax.plot(*(pixoffs[det].T*impixsize*1e6), '.')
    ax.text(0,0,f'{det}')
    ax.plot(*imeas_im_offset_by_det_m[det,:]*1e6, 'r+', markersize=15)
    ax.plot(*(imeas_im_offset_by_det_m[det,:]+meas_alignment_corr)*1e6, 'rs', markersize=15, markerfacecolor='none')
pass
# %%

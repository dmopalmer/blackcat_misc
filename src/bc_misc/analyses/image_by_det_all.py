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
    balance=balance_boxes
)


if True:
    ph_all = sorted(set([line.split()[0].split("_")[0] for line in sourcedir.joinpath('found_scox1.txt').open('rt').readlines()]))
else:
    ph_all = (sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_05_[23]*/level1_tmp/ph_2*")) +
        sorted(top_dir.rglob("payload_checkout/CHECKOUT_26_06_*/level1_tmp/ph_2*")))
    
#%%
plot=False
# Corrective offsets to add to detector detx, dety
detpixsize = 40e-6
impixsize = imager.resolution * detpixsize

# The gross offset which is probably due to instrument alignment of ~0.47°
# comparable to the ~0.5° offset found between the two star camera nominals.
# alignment_corr = np.array(  [0,0])  # 2605201842_d0 ij = 1155.0,598.0 = ScoX1+ -13.7,  -5.0 ra,dec = 215.31,-35.29  10.3 σ
alignment_corr = np.array([-0.00115047, -0.00039763])  # 
# How much to add to each detx, dety
# after iter 1
# corr_by_det =\
#     np.array([
#     [    20.2e-6,   -182.5e-6],
#     [    28.4e-6,    148.0e-6],
#     [  -230.3e-6,    -64.1e-6],
#     [   181.7e-6,     98.6e-6],
#     ])
# after iter 2
# corr_by_det =\
#     np.array([
#     [    13.6e-6,   -203.4e-6],
#     [    48.0e-6,    184.2e-6],
#     [  -236.5e-6,    -99.3e-6],
#     [   174.9e-6,    118.5e-6],
#     ])
# After iter 3 (resolution 4)
# corr_by_det =\
#     np.array([
#     [    41.0e-6,   -183.0e-6],
#     [   -34.7e-6,    173.2e-6],
#     [  -208.6e-6,   -113.1e-6],
#     [   202.2e-6,    122.9e-6],
#     ])
# After iter 4 (resolution now 2)
corr_by_det =\
    np.array([
    [    22.6e-6,   -183.4e-6],
    [    16.7e-6,    176.2e-6],
    [  -222.4e-6,    -87.4e-6],
    [   182.9e-6,     94.6e-6],
    ])
# After 5 iterations:
corr_by_det =\
    np.array([
    [    56.0e-6,   -198.3e-6],
    [   -22.2e-6,    166.8e-6],
    [  -253.3e-6,    -82.5e-6],
    [   219.3e-6,    114.0e-6],
    ])
# After 6 iterations
corr_by_det =\
    np.array([
    [    23.9e-6,   -204.9e-6],
    [     2.6e-6,    157.4e-6],
    [  -212.9e-6,    -80.9e-6],
    [   186.2e-6,    128.4e-6],
    ])
# After 7 iterations, the last at resolution=1
corr_by_det =\
    np.array([
    [    42.8e-6,   -177.0e-6],
    [   -11.5e-6,    180.1e-6],
    [  -222.7e-6,    -88.0e-6],
    [   191.1e-6,     84.9e-6],
    ])
# After 8 iterations and 200 detections
corr_by_det =\
    np.array([
    [    49.2e-6,   -185.8e-6],
    [   -11.0e-6,    149.4e-6],
    [  -229.5e-6,    -76.7e-6],
    [   191.1e-6,    113.2e-6],
    ])
"""
Per-detector residual correction =
[[-1.76133614e-05 -7.03470001e-06]
 [ 1.99810275e-05 -8.49028884e-06]
 [-1.07802612e-05 -1.50402480e-06]
 [ 8.41259508e-06  1.70290136e-05]]
Added to current correction = 
[[ 4.280e-05 -1.770e-04]
 [-1.150e-05  1.801e-04]
 [-2.227e-04 -8.800e-05]
 [ 1.911e-04  8.490e-05]]
Gives a new correction of
corr_by_det =\
    np.array([
    [    25.2e-6,   -184.0e-6],
    [     8.5e-6,    171.6e-6],
    [  -233.5e-6,    -89.5e-6],
    [   199.5e-6,    101.9e-6],
    ])
"""
# After 8 iterations
corr_by_det =\
    np.array([
    [    25.2e-6,   -184.0e-6],
    [     8.5e-6,    171.6e-6],
    [  -233.5e-6,    -89.5e-6],
    [   199.5e-6,    101.9e-6],
    ])
# After 9 iterations
corr_by_det =\
    np.array([
    [    55.9e-6,   -186.1e-6],
    [     1.0e-6,    167.5e-6],
    [  -230.8e-6,    -86.6e-6],
    [   173.7e-6,    105.3e-6],
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
            # WHY NEGATIVE?
            counts['DETY'] -= corr[1]
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
det_meas_corr_m = []
for i in range(4):
    pixoffs[i] = np.array(pixoffs[i])
    det_meas_corr_m.append(biweight_location(pixoffs[i], axis=0) * impixsize) 
    
det_meas_corr_m = np.array(det_meas_corr_m)
meas_alignment_corr = np.mean(det_meas_corr_m, axis=0)
det_meas_corr_m -= meas_alignment_corr

new_corr_by_det = det_meas_corr_m + corr_by_det

print(f"Measured alignment correction = {meas_alignment_corr * 1e6} u + previous {alignment_corr * 1e6}u = {meas_alignment_corr + alignment_corr}")

print(f"Per-detector residual correction =\n{det_meas_corr_m}")
print(f"Added to current correction = \n{corr_by_det}")
print(f"Gives a new correction of\ncorr_by_det =\\\n    np.array([")
for d in range(4):
    print(f"    [{new_corr_by_det[d,0]*1e6:8.1f}e-6, {new_corr_by_det[d,1]*1e6:8.1f}e-6],")
print("    ])")
fig,axes = plt.subplots(2,2, figsize=[10,10], sharex=True, sharey=True)
axes = axes.ravel()[[1,2,0,3]]
for det,ax in enumerate(axes):
    ax.plot(*(pixoffs[det].T*impixsize*1e6), '.')
    ax.text(0,0,f'{det}')
    ax.plot(*det_meas_corr_m[det,:]*1e6, 'r+', markersize=15)
    ax.plot(*(det_meas_corr_m[det,:]+meas_alignment_corr)*1e6, 'rs', markersize=15, markerfacecolor='none')
pass
# %%

"""
imaging.py
"""
# %%
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *

from tools import *
from xraysky import BC_Imager


# From Tim
BASE_PATH = top_dir / "payload_checkout/CHECKOUT_26_05_22"
LEVEL1_PATH = BASE_PATH / "level1_tmp"
ph_file = LEVEL1_PATH / "ph_2605221058.gz"
ori_file = LEVEL1_PATH / "or_2605221058.gz"

orientations = fits.getdata(ori_file, 1)
orientation_start_time = np.min(orientations["TIME"])

stable_time = orientation_start_time + 150

stable_orientations = orientations[orientations["TIME"] >= stable_time]
stable_pointings = stable_orientations["POINTING_NOMINAL"]
stable_ra = np.median(stable_pointings[:,0])
stable_dec = np.median(stable_pointings[:,1])
stable_roll = np.median(stable_pointings[:,2])

with fits.open(ph_file) as event_hdul:
    counts, header = fits.getdata(ph_file, header=True)
    counts = counts[counts["TIME"] >= stable_time]

    balance = []
    for det in range(4):
        d = counts[counts['DETID'] == det]
        if len(d) > 10:
            balance.append([[d['DETX'].min(), d['DETY'].min()], [d['DETX'].max(), d['DETY'].max()]])
    balance = np.array(balance)
    imager = BC_Imager(
        resolution=2,
        radecroll_inst=[
            stable_ra,
            stable_dec,
            -stable_roll,
        ],
        overwrite = True,
        balance=balance
    )
    imhdu = imager.evtlist2image(
        counts=counts, header=header,
        add_wcs=True,
        outfile=Path("/tmp/ph_2605221058_image_neg_roll.fits.gz"),
    )
    plt.imshow(imhdu.data)
    pass
#%%
imager = BC_Imager(
    resolution=1,
    radecroll_inst=[
        stable_ra,
        stable_dec,
        -stable_roll,
    ],
    overwrite = True,
    balance=balance_boxes
)
fig,axes = plt.subplots(2,2, sharex=True, sharey=True, figsize=(12,8))
axes = np.ravel(axes)[[1,2,0,3]]
for det,ax in enumerate(axes):
    imhdu = imager.evtlist2image(counts=counts[counts['detid']==det], header=header)
    im=ax.imshow(imhdu.data)
    ax.set(title=f"Det={det} peak {imhdu.data.max():.0f} ± {imhdu.data.std():.0f}")
fig.suptitle(ph_file.stem)
fig.tight_layout()
# plt.colorbar(im=im,axes=axes)
# %%

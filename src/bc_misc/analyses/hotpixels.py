#%%
"""
Use the detector calibration readout of hot pixels to
make an fpga threshold map that excludes the worst pixels.

"""



from tools import *

f_detcal=Path('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_03_27/raw_level0/mar26cal_1774623516.gz')
f_gains = Path('/Users/palmer/repo/blackcat/bc_misc/src/bc_misc/analyses/gainsMar26.gz')
f_gains_new=Path('/Users/palmer/repo/blackcat/bc_misc/src/bc_misc/analyses/gainsMar31.gz')
fits.info(f_detcal)
fits.info(f_gains)
# %%
d_detcal=np.stack([fits.getdata(f_detcal, ("hotness",det)) for det in range(4)])
d_gains_fpgathresh=np.stack([fits.getdata(f_gains, ("fpgathresh",det)) for det in range(4)])
d_gains_softwarethresh=np.stack([fits.getdata(f_gains, ("softwarethresh",det)) for det in range(4)])
# %%
%matplotlib ipympl

images = []
plt.close('all')
fig,axes = plt.subplots(3,4, figsize=(12,8), sharex=True, sharey=True)
axes=np.ravel(axes)
for i,ax in enumerate(axes[:4]):
    images.append(ax.imshow(d_detcal[i]))
for i,ax in enumerate(axes[4:8]):
    images.append(ax.imshow(d_gains_fpgathresh[i]))
for i,ax in enumerate(axes[8:12]):
    images.append(ax.imshow(d_gains_softwarethresh[i]))

ax.set(aspect=1)
# fig.colorbar(images[-1], ax=axes,orientation='vertical')
fig.suptitle(f"Hotness of pixels from {f_detcal.name}")
fig.tight_layout()

# plt.close('all')
# figt,axest = plt.subplots(2,2, figsize=(12,12), sharex=True, sharey=True)
# axest=np.ravel(axest)
# for i,ax in enumerate(axest):
#     images.append(ax.imshow(d_gains_fpgathresh[i]))
# ax.set(aspect=1)
# fig.colorbar(images[0], ax=axes,orientation='vertical')
# fig.suptitle("Hotness of pixels from {}")


# %%

plt.close('all')
sortvals=d_detcal.ravel().copy()
sortvals.sort()
plt.plot(sortvals)
# %%
plt.close('all')
plt.plot(sortvals,np.cumsum(sortvals)/np.sum(sortvals))
# %%
d_hotthresh = np.where(d_detcal >= 252, 65535, 0)
new_fpga_thresholds = np.max(np.stack((d_gains_fpgathresh, d_gains_softwarethresh, d_hotthresh), axis=0),axis=0)

# %%
plt.close('all')
fig,axes = plt.subplots(4,4, figsize=(12,12), sharex=True, sharey=True)
axes=np.ravel(axes)
for i,ax in enumerate(axes[:4]):
    images.append(ax.imshow(d_detcal[i]))
for i,ax in enumerate(axes[4:8]):
    images.append(ax.imshow(d_gains_fpgathresh[i]))
for i,ax in enumerate(axes[8:12]):
    images.append(ax.imshow(d_gains_softwarethresh[i]))
for i,ax in enumerate(axes[12:16]):
    images.append(ax.imshow(new_fpga_thresholds[i]))

ax.set(aspect=1)
# fig.colorbar(images[-1], ax=axes,orientation='vertical')
fig.suptitle(f"Hotness of pixels from {f_detcal.name}")
fig.tight_layout()

# %%

"""
Create new thresholds file
"""

hdulist = []

for hdu in fits.open(f_gains):
    hdu = hdu.copy()
    if hdu.header['EXTNAME'].lower() == "fpgathresh":
        n_orig = np.count_nonzero(hdu.data >= 3000)
        detnum = int(hdu.header['EXTVER'])
        hdu.data = new_fpga_thresholds[detnum]
        hdu.header.add_comment("Set FPGA threshold based on hot pixels and software threshold 2026-03-30")
        n_final = np.count_nonzero(hdu.data >= 3000)
        print(f"Detector {detnum} {n_orig} -> {n_final} pixels supressed in FPGA")
    hdulist.append(hdu)

hdul = fits.HDUList(hdulist)
hdul.writeto(f_gains_new)


# %%

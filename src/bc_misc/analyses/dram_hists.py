# %%
"""
Histograms from DRAM dumps
"""

# from tools import *
from mylab import *
f=Path('/Users/palmer/Downloads/dram_nomApr18te-with-histograms.fits')


hists = np.concatenate([fits.getdata(f, ('dethist', i)) for i in range(4)])

plt.plot(hists.sum(axis=0))

# %%
plt.close('all')
for det in range(4):
    hist = np.mean(fits.getdata(f, ('dethist', det)), axis=0)
    plt.plot(hist, label =f"det {det}")

plt.legend()
plt.gca().set(yscale='log', title=f"{f.name} histograms", ylabel="Pixels per readout, averaged", xlabel="Raw ADC value")
plt.savefig(f'/tmp/{f.name}.png')
# %%
plt.close('all')
spacing=5
det=3
d=fits.getdata(f, ('dethist', det))
for i,row in enumerate(d):
    plt.plot(row * spacing ** -i)
plt.gca().set(yscale='log', title=f"det = {det}  {f.name} histograms")

# %%
hrandom = fits.getheader('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_04_22/raw_level0/ph_TH1940Apr21qk.gz',1)
framerate = hrandom['TS_RATTK']/hrandom['TS_FRMPR']
framerate

# %%
plt.close('all')
for det in range(4):
    hist = np.mean(fits.getdata(f, ('dethist', det)), axis=0)
    cumhist = np.cumsum(hist[::-1])[::-1]
    plt.plot(cumhist * framerate, label =f"det {det}")

plt.legend()
plt.gca().set(yscale='log', title=f"cumulative {f.name} rates", ylabel="Pixels/second >ADC", xlabel="Raw ADC value", ylim=[1e5,cumhist.max()*framerate*2])
plt.savefig(f'/tmp/cum_{f.name}.png')

# %%

"""
Distributions of hots and biases
"""
#%%
import matplotlib
matplotlib.use('MacOSX')  # Use 'osx' in some older IPython environments

from tools import *
_ion = plt.ion()
f=Path('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_06_18/raw_level0/dc_Jun18bk.gz')

plt.figure() 
for i in range(4):
    d=fits.getdata(f, 2*i+1)
    d=np.ravel(d)
    d.sort()
    plt.plot(d,np.arange(len(d)))
plt.gca().set(title=f"{Path(f).stem} Bias values", xlabel="Pixel bias", ylabel="$N_{pixels} >$")
plt.savefig("/tmp/bias.png")

plt.figure() 
for i in range(4):
    d=fits.getdata(f, 2*i+2)
    d=np.ravel(d)
    d.sort()
    plt.plot(d,np.arange(len(d)))
plt.gca().set(title=f"{Path(f).stem} Hotness values", xlabel="Pixel hotness", ylabel="$N_{pixels} >$")
plt.savefig("/tmp/hotness.png")


# %%

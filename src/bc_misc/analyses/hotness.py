#%%
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *
import blackcat
from collections import Counter
dir_latest = sorted(top_dir.glob("payload_checkout/CHECKOUT_26*"))[-1]

print(f"Latest data is {dir_latest.name}")
outdir = new_outdir()

# %%
plt.close('all')
plt.figure(figsize=[10,10])
dc_files = sorted(dir_latest.glob('raw_level0/dc_*'))
f=dc_files[-2]
dbias=stack_extension(f, 'bias')
dhot = stack_extension(f, 'hotness')
plt.imshow(dhot, origin='lower')
plt.colorbar()

# %%

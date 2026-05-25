# %%
from tools import *

fnames = get_ph_files("260326*", level=1, no_x=False)
fnames
# %%
#% matplotlib widget
for f in fnames:
    d = fits.getdata(f)
    plt.figure()
    plt.plot(d['time'], d['rawx'] + 2000 * d['detid'].astype(int), ',')
    plt.suptitle(f.name)
    plt.savefig(outdir().joinpath(f"{f.name}_bytime.png"))
# %%

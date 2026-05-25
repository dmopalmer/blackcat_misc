#%%
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *
import blackcat
from collections import Counter
dir_latest = sorted(top_dir.glob("payload_checkout/CHECKOUT_26*"))[-1]

print(f"Latest data is {dir_latest.name}")
outdir = new_outdir()

# ALL starts in April, and is only the most recent 20
if False:
    ph1_list = sorted(top_dir.glob("payload_checkout/CHECKOUT_26_04_*/level1*/ph_*"))[-20:]
else:
    # Or just all today
    ph1_list = sorted(dir_latest.glob("level1*/ph_*"))
ph1_byname = {fname.stem.split('_')[1] : fname for fname in ph1_list}
try:
    ph1_byname.pop('TOO')
    # startup has the detectors in the wrong places
    # ph1_byname.pop('startup')
    for name in list(ph1_byname):
        if '2603161733x' <= name and name <= '2603171926x':
            ph1_byname.pop(name)
except Exception as e:
    # print(name, e)
    pass

# %%
plt.close('all')


for name, fname in ph1_byname.items():
    fig,d,h = plot_photons_isl4(name, fname)
    d_sub19 = d[d['ENERGY'] < 19]
    print(f"{name:20} {len(d)} {100 * len(d_sub19)/len(d):5.1f}% below 19 keV {detsplit(d, apply=len)} per detector in {np.ptp(d['TIME']):.1f} s")
    fig.savefig(outdir.joinpath(f"{name}_photons.png"))
    # c.update(zip(d['detid'],d['rawx'], d['rawy']))
    pass
# print(f"Most common: {c.most_common(20)}")

# %%
plt.close('all')
plt.plot(d['TIME'],',')
detsplit(d, apply=len)


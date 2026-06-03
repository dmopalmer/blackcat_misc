#%%
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *
import blackcat
from collections import Counter
dir_latest = sorted(top_dir.glob("payload_checkout/CHECKOUT_26*"))[-1]

ph_list = sorted(top_dir.glob('payload_checkout/CHECKOUT_26_05_25/level1_tmp/ph_*.gz'))
print(len(ph_list))
# %%
plt.close('all')
fig,axes = plt.subplots(nrows=2, figsize=(12,6), sharex=True)
for ph in ph_list:
    d = fits.getdata(ph)
    tedges,rates = ph_rate(d)
    axes[0].stairs(rates, tedges, color='b')
axes[0].set(yscale='log', ylim=[0.5,None])
# %%

# FIND LOOPS
recent_ph = sorted(top_dir.glob('payload_checkout/CHECKOUT_26_05_23/level1_tmp/ph_*.gz'))
for ph in recent_ph:
    d=fits.getdata(ph)
    if loopcheck(d):
        print(f'Loop in {ph}')
        # break
    else:
        print(f'{ph.name} OK {len(d):6d} / {np.ptp(d['TIME']):.1f} s = {len(d) /np.ptp(d['TIME']):6.0f} cps in {len(set(d['detid']))} detectors')
# %%
dloop = fits.getdata('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_05_23/level1_tmp/ph_2605191543.gz')
dlp = dloop[(dloop['framenum'] == 9) & (dloop['detid'] == 3)]
d3 = dloop[dloop['detid'] == 3]
print(len(dlp))
# %%
for f in Path('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_05_23/level1_tmp/').glob('ph_*.gz'):
    plot_photons_isl4(f.name, f)
# %%


"""
Look at when voltages fail to turn on

2026-05-15T06:30 is one case where det0 voltages don't turn on, and dets 3, 4 vs2mon have strange behavior
2026-05-14T12:30
"""

#%%

import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *


ds = bcd.DataSource(
    list(bcd.topdir.rglob("**/CHECKOUT_26_05_1[3456789]*/**/raw_level0/**/HOUSEKEEPING*"))
)

lines=ds.messages(fmt=True)
Path("messages.txt").open("w").write("\n".join(lines) + "\n")

regs = ds.fpga_registers(namecolumn='value')
# %%

voltnames = [name.split('_')[0] for name in ds.names('v*_0')]

nvolts = len(voltnames)
vdata = ds.data('vcds_0')[0]

plt.close('all')
fig, axes = plt.subplots(nrows=nvolts, ncols=1, figsize=(12,7), sharex=True)
all_axes = np.ravel(axes)

for vname,ax in zip(voltnames, all_axes):
    for det,color in enumerate('rgbc'):
        ax.plot('t', f'{vname}_{det}', '-', color=color, data=vdata)
        if det==0:
            ax.set(ylabel=vname)
fig.tight_layout()
# %%

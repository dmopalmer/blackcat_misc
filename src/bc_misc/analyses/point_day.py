#%%
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *


ds = bcd.DataSource(
    list(bcd.topdir.rglob("**/CHECKOUT_26_06_18/**/raw_level0/**/HOUSEKEEPING*"))
)


# %%
rolls = ds.data('point_roll')[0]
plt.plot(rolls['t'], rolls['point_roll'],'.')
# %%

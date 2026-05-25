#%%
from tools import *

daystring="26_03_31"
orientation_file=Path(f'/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_{daystring}/raw_level0/HOUSEKEEPING_ORIENTATION.fits.gz')
fits.info(orientation_file)

# %%
d=fits.getdata(orientation_file, 'ORIENTATION_TLM')
d.names
# %%
%matplotlib ipympl
plt.plot('time', 'POINT_roll', '.', data=d)

# %%

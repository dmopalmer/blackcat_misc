from tools import *
from mylab import *
import blackcat


f0=Path('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_04_13/raw_level0/ph_2604102351.gz')

d0,h0 = fits.getdata(f0, header=True)

d1,h10,h11,_ = blackcat.convert_events_level0_to_level1(f0, observation_id='20'+f0.name[3:13], overwrite=True)
pass
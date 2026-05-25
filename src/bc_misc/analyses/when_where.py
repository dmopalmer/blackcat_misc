# %%
"""
When where
"""

from tools import *

tstamps = [1776452366,
            1776484123,
            1776489515,
            1776539539,
            1776559911,]

sat = get_satellite()

for ts in tstamps:
    when = datetime.fromtimestamp(ts)
    pos = satellite_location(when, sat)[0]
    print(ts, when.isoformat()[:19], pos.latitude.degrees, "N, ", pos.longitude.degrees, "E")
# %%


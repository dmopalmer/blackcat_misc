"""
When is BlackCAT in eclipse?

"""
#%%
from tools import *
from mylab import *
import numpy as np
import skyfield
import skyfield.api as sfapi
from skyfield import almanac
from datetime import timedelta
import xraysky


eph = skyfield.load('de421.bsp')

tle_20260619 = [
        "0 BLACKCAT",
        "1 67369U 26004G   26168.62496213  .00001273  00000-0  14572-3 0  9994",
        "2 67369  97.7980 167.0599 0008293 133.9176 226.2726 14.87661774 23353"
]

sat = sfapi.EarthSatellite(tle_20260619[1], tle_20260619[2], tle_20260619[0])

def is_sunlit(t):
    return sat.at(t).is_sunlit(eph)
# check on a 1 minute cadence
is_sunlit.step_days=1/1440

trange = xraysky.sftime([parsedate("2026-03-01"), parsedate("2026-08-31")])

t_transition,is_sun_ = almanac.find_discrete(trange[0], trange[1], is_sunlit)
#%%
# Fail if we start eclipsed
fig = plt.figure()
assert(is_sun_[0] == 0)
t = t_transition.utc_datetime()
tstart,dur = zip(*[(t_in, (t_out - t_in).total_seconds()/60) for t_in,t_out in zip(t[::2], t[1::2])])
with open('/tmp/eclipses.txt', "w") as f:
    for tstart_, dur_ in zip(tstart,dur):
        print(f"{tstart_:%Y-%m-%dT%H:%M:%S} {dur_:5.2f}", file=f)
plt.plot(tstart, dur, '.')
ax = plt.gca()
ax.set(title="BlackCAT Eclipses", ylabel="Minutes of eclipse per orbit", xlabel ="Date (2026)")
fig.autofmt_xdate()
fig.tight_layout()
plt.savefig('/tmp/eclipses.png')

# %%

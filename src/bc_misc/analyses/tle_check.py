"""
Compare GPS data to TLEs

"""

from mylab import *
import numpy as np
import skyfield
import skyfield.api as sfapi
from datetime import timedelta

gpslogfile = Path("/Users/palmer/Downloads/20260204_GPS_Log_content.log")
tlefile = Path(__file__).parent.parent.parent.parent.joinpath(
    "data/orbital/gps_v_st.txt"
)
sats = sfapi.load.tle_file(str(tlefile))

# Conversion from gps_Seconds_of_week to utc:
"""Example Lines
$GPGSV,3,3,10,14,11,217,34,04,06,003,41*7A
#BESTXYZA,COM1,0,81.0,FINESTEERING,2402,265220.000,0200400c,44cf,16410;SOL_COMPUTED,SINGLE,3917099.1487,5488322.4121,-1837626.9099,4.1151,5.6238,3.5977,SOL_COMPUTED,DOPPLER_VELOCITY,122.5236,-2500.9525,-7215.1684,0.4919,0.6722,0.4300,"",0.000,0.000,0.000,27,24,24,0,0,02,01,11*106f6616
$GPGGA,014002.00,1520.0391,S,05429.0384,E,1,24,0.5,612090.04,M,-21.40,M,,*6F
#BESTPOSA,COM1,0,81.0,FINESTEERING,2402,265220.000,0200400c,cdba,16410;SOL_COMPUTED,SINGLE,-15.33398463219,54.48397270322,612090.0437,-21.4000,WGS84,3.1198,3.2600,6.4144,"",0.000,0.000,27,24,24,0,00,02,01,11*822be753
$GPRMC,014002.00,A,1520.0390779,S,05429.0383622,E,14845.710,191.7,210126,0.0,E,A*2C
$GPVTG,191.730,T,191.730,M,14845.710,N,27494.256,K,A*24
$GPGSA,M,3,08,03,06,16,04,09,07,32,26,17,14,27,0.8,0.5,0.6*36
$GPGSV,3,1,10,01,62,218,,03,49,013,42,02,43,172,,08,37,132,36*76
"""
# From $GPRMC and BESTXYZA lines
t0 = parse_as_utc("2026-01-21T01:40:02") - timedelta(seconds=265220) + timedelta(days=7)
week_secs = 86400.0 * 7
# BESTXYZA,COM1,0,81.0,FINESTEERING,2402,265220.000,0200400c,44cf,16410;SOL_COMPUTED,SINGLE,3917099.1487,5488322.4121,-1837626.9099,4.1151,5.6238,3.5977,SOL_COMPUTED,DOPPLER_VELOCITY,122.5236,-2500.9525,-7215.1684,0.4919,0.6722,0.4300,"",0.000,0.000,0.000,27,24,24,0,0,02,01,11*106f6616
#  0        1  2  3     4            5    6          7        8    9                  10     11           12             13
with Path(gpslogfile).open() as f:
    txyz = np.array(
        [
            [l[6], l[11], l[12], l[13]]
            for l in [
                l.split(",")
                for l in Path(gpslogfile).open().readlines()
                if "FINESTEERING" in l and "SOL_COMPUTED,SINGLE" in l
            ]
            if l[0] == "#BESTXYZA" and l[4] == "FINESTEERING"
        ]
    ).astype(float)
    # Last ~60 elements or so are unreliable
    txyz = txyz[:-70]

ts = skyfield.load.timescale(builtin=True)
t = ts.from_datetime(t0) + np.unwrap(txyz[:, 0], period=week_secs) / 86400.0

mindist = 1e9
bestsat = None
meandists = []
for sat in sats:
    pos = sat.at(t).itrf_xyz().m.transpose()
    dist = np.sqrt(np.sum((pos - txyz[:, 1:]) ** 2, axis=1))
    meandist = np.mean(dist)
    meandists.append(meandist)
    if True or meandist < 1e6:
        print(f"{sat.name:15} Mean distance: {meandist / 1000:10.3f}")
        plt.plot(t.utc_datetime(), dist / 1000.0, ".", label=sat.name)
    if meandist < mindist:
        mindist = meandist
        bestsat = sat
plt.legend()
print(bestsat)

pass

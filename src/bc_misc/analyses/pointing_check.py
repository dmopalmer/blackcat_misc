"Checking pointing with quaternions"

"""
Example from the outstanding questions tracker
https://nanoavionics.sharepoint.com/:x:/r/sites/BLACKCAT052/_layouts/15/Doc.aspx?sourcedoc=%7BB4F6C249-A369-4B0E-86F4-8C54DF6A272C%7D&file=BlackCAT%20Question%20Tracker.xlsx&action=default&mobileredirect=true

Star Target  ECI J2000
	
Pointing Target
	
ECI Target Quaternion J2000
	
Performance


Betelgeuse 
RA: 89.1521 deg
DEC: 7.4110 deg
	
17234
Payload Pointing Option #2 (-X to Star Vector, -Y to Nadir)
(Star Pointing J2000)
	
Meas. Start:
2026-02-07 11:24:00 
q(w,x,y,z) = (-0.489, -0.503, 0.443 ,  0.558)
Meas. Mid:
2026-02-07 11:29:00 
q(w,x,y,z) = (-0.422, -0.561, 0.510, 0.498)
Meas. End:
2026-02-07 11:34:00 
q(w,x,y,z) = (-0.342, -0.613, 0.572,  0.442)



	
Control Error [deg]
CE mean: 0.01882
CE 1 sigma: 0.02207
CE 2 sigma: 0.03455
CE 3 sigma: 0.04324
Attitude Drift [deg/s]
wCE mean: 0.008023
wCE 1 sigma: 0.01009
wCE 2 sigma: 0.02004
wCE 3 sigma: 0.03141

"""


"""
Conclusions:
The last quaternion probably has a typo (the norm of the quaternion is 1.0153 and pointing directions are off by of order a degree).
But for the other two pointings my software agrees that according to the quaternions:
The instrument boresight (spacecraft –X) points at Betelgeuse (~0.07° error).
The instrument boresight does not mispoint by the ~0.5 degrees offset you would get by confusing J2000 and current-equinox coordinates.
The Earth direction (Nadir) is as near as possible (10s of degrees) to the Spacecraft -Y Axis while still pointing -X at Betelgeuse
More specifically, the Nadir is in the spacecraft X-Y plane (within 0.02 degrees)

Therefore we have verified that the PSU understanding of the quaternions matches that of KNA, and this question can be closed.

"""
from mylab import *
import xraysky
import astropy
import skyfield.api as sfapi
import skyfield
from quaternion import quaternion
import quaternion as quat
from skyfield.positionlib import ICRF


radec_star = [89.1521, 7.4110]
ts = sfapi.Loader(os.path.expanduser("~/.config/skyfield")).timescale(builtin=True)

satellite = sfapi.EarthSatellite(
    "1 67369U 26004G   26040.94764658  .00003294  00000-0  37099-3 0  9993",
    "2 67369  97.7966  41.5893 0006825 191.9283 168.1773 14.87297113  4371",
    name="BlackCAT",
)

times = ts.from_datetimes(
    [
        parse_as_utc(t)
        for t in ["2026-02-07 11:24:00", "2026-02-07 11:29:00", "2026-02-07 11:34:00"]
    ]
)
measquats = [
    quaternion(-0.489, -0.503, 0.443, 0.558),
    quaternion(-0.422, -0.561, 0.510, 0.498),
    quaternion(-0.342, -0.613, 0.572, 0.442),  # This one has a bad norm
]

inst = xraysky.BlackCAT()

satlocs = satellite.at(times)
nadir_vectors = -satlocs.position.au
nadir_pos = ICRF(nadir_vectors, t=times)

for q, loc, nadir_pos in zip(measquats, satlocs, nadir_pos):
    inst.orientation = xraysky.Orientation(q_sc=q)
    nadir_ra, nadir_dec, nadir_dist = nadir_pos.radec()
    theta, phi = inst.orientation.radec2thetaphi(nadir_ra.degrees, nadir_dec.degrees)
    print(
        f"ra,dec,roll = {inst.orientation.radecroll()}  {q.norm()=:.5} Earth ra,dec = {nadir_ra.degrees:.2}, {nadir_dec.degrees:.2} theta (dist) = {np.rad2deg(theta)} phi (posang) = {np.rad2deg(phi)}"
    )


pass

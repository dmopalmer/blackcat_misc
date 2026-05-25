# %%
from tools import *
import shutil

old_to_right_mapping = [2, 0, 3, 1]

f_gainsMar13 = f = (
    "/Users/palmer/Downloads/2026-03-11T20_57_22Z_science-fixes1/assets/gains.gz"
)
f_output = "./gainsMar23.gz"

hdulist = []

for hdu in fits.open(f_gainsMar13):
    hdu = hdu.copy()
    hdu.header["EXTVER"] = old_to_right_mapping[hdu.header["extver"]]
    hdulist.append(hdu)

hdul = fits.HDUList(hdulist)
hdul.writeto(f_output)

rightver = "foo"

# shutil.copy(f_gainsMar13, f_output)

for i in range(16):
    h = fits.getheader(f_output, i)
    try:
        # rightver = old_to_right_mapping[h['extver']]
        print(
            f"{i:2} {h['extname']:15} is {h['extver']} should be {rightver} {h.get('DETS_SN', 'No SN')}"
        )
    except Exception as e:
        print(f"{i:2} {e}")
# %%

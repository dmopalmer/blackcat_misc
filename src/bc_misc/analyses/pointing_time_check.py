"""
Compare pointing from the orientation housekeeping (clocked by the CFE clock)
with pointing from the OR files (clocked by the GPS)

"""

# %%
import matplotlib as mpl

mpl.use("module://ipympl.backend_nbagg")
import bc_misc
from bc_misc.analyses.tools import *
from mylab import *
import blackcat
import blackcat_data as bcd

day = parse_as_utc("2026-11-21 12:00")

daydir = top_dir.joinpath(f"payload_checkout/CHECKOUT_{day:%y_%m_%d}")

l1_orfiles = sorted(daydir.glob("level1_tmp/or*"))
datasource = bcd.DataSource(daydir.rglob("**/HOUSEKEEPING*"))
if True:
    l1_orfiles = sorted(
        top_dir.glob(f"payload_checkout/CHECKOUT_26_11_*/level1_tmp/or*")
    )
    datasource = bcd.DataSource(
        top_dir.glob(f"payload_checkout/CHECKOUT_26_11_*/**/HOUSEKEEPING*")
    )


dpt = datasource.data("point_ra")[0]

# %%
plt.close("all")
fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(12, 12))
(ax_ra, ax_dec, ax_roll) = axes
for field, ax in zip(["point_ra", "point_dec", "point_roll"], axes):
    ax.plot("t", field, ".-r", data=dpt)
# %%

for f in l1_orfiles:
    d = fits.getdata(f)
    t = ts_to_time(d["time"])
    for coord, ax in zip(d["POINTING_NOMINAL"].T, axes):
        # plot the line
        ax.plot(t, coord, "-b")
        # mark its endpoints
        ax.plot(t[[0, -1]], coord[[0, -1]], "|k")
fig.suptitle(
    "Do the red points on the plot (coord vs CFE Time) match the blue curve (coord vs GPS time)?"
)
fig.tight_layout()
# %%

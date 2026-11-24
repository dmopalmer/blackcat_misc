"""
How fast is the LTC FPGA clock ticking compared ot GPS PPS,
and what do the residuals look like?

"""

# %%
import matplotlib as mpl

mpl.use("module://ipympl.backend_nbagg")
import bc_misc
from bc_misc.analyses.tools import *
from mylab import *
import blackcat
import blackcat_data as bcd

day = parse_as_utc("2026-11-17 12:00")

daydir = top_dir.joinpath(f"payload_checkout/CHECKOUT_{day:%y_%m_%d}")

l1_orfiles = sorted(daydir.glob("level1_tmp/or*"))
datasource = bcd.DataSource(daydir.rglob("**/HOUSEKEEPING*"))
d_timectl = datasource.data("pps_age")[0]

orheaders = Headers(l1_orfiles)

times = orheaders.value("TS_TIME", as_time=True)
seconds = orheaders.value("TS_TIME")
ticks = orheaders.value("TS_PPSTK")
ijump = np.ravel(np.argwhere(np.diff(ticks) < 1))[-1] + 3
times, seconds, ticks = [v[ijump:] for v in (times, seconds, ticks)]
print(f"{times[0].item().isoformat()}")

# %%
plt.close("all")
fig, axes = plt.subplots(nrows=3, sharex=True, figsize=[10, 10])


slope, intercept = np.polyfit(seconds, ticks, 1)
y_fit = slope * seconds + intercept
residuals_ticks = ticks - y_fit
residuals = residuals_ticks / slope

axes[0].plot(times, ticks, ".")
axes[0].set(ylabel="Ticks")
axes[1].plot(times, residuals_ticks / slope)
axes[1].set(ylabel="Residual (s)")

w = (residuals > (-0.5 * slope)) & (residuals_ticks < (0.5 * slope))
slope_w, intercept_w = np.polyfit(seconds[w], ticks[w], 1)
y_fit_w = slope_w * seconds[w] + intercept_w
residuals_w = ticks[w] - y_fit_w

axes[2].plot(times[w], residuals_w / slope_w)


fig.suptitle(f"Ticks fit to {slope_w} ticks/s")


# %%

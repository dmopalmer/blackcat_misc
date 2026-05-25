# %%
# % matplotlib widget
from tools import *

# Photon file that shows a byte shift
badname = "2603241324x"
# %%

fl0 = sorted(top_dir.rglob(f"raw_level0/*{badname}*"))[-1]
fl1 = sorted(top_dir.rglob(f"level1*/*{badname}*"))[-1]
d0 = fits.getdata(fl0)
d1 = fits.getdata(fl1)
# %%
plt.close("all")
plt.plot("detx", "dety", ",", data=d1)
plt.gca().set(aspect=1)

# %%
plt.close("all")
plt.plot("time", "rawy", ",", data=d1)
# %%
plt.close("all")
plt.plot(d1["time"], ",")
plt.gca().set(xlabel="Record number", ylabel="time", title=fl1.name)
plt.savefig("/tmp/time_problem.png")
# %%
column = "bias"
plt.close("all")
plt.plot(d1[column], ",")
plt.gca().set(xlabel="Record number", ylabel=column, title=f"{fl1.name} {column}")
# plt.savefig(f"/tmp/{column}_problem.png")

# %%
for column in d0.names:
    plt.close("all")
    plt.plot(d0[column], ",")
    plt.gca().set(
        xlabel="Record number", ylabel=column, title=f"d0 {fl1.name} {column}"
    )
    plt.savefig(f"/tmp/d0_{column}_problem.png")

for col_left, col_right in zip(d0.names, d0.names[1:]):
    plt.close("all")
    plt.plot(((d0[col_left] & 0xFF) << 8) | ((d0[col_right] & 0xFF00) >> 8), ",")
    plt.gca().set(
        xlabel="Record number",
        ylabel=column,
        title=f"d0 {fl1.name} {col_left} smush {col_right}",
    )
    plt.savefig(f"/tmp/d0_{col_left}_{col_right}_problem.png")


# %%

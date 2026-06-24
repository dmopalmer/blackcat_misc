# %%
"""
Look at the new registers that count how often we fail-back to global threshold

This should have added the following changes - 
added counter of processed pixels with threshold and fallback.
- # of pixels with successful local thresholding
registers Xdc,Xdb
- # of pixels with fallback global thresholding
registers Xda,Xd9

"""
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *


ds = bcd.DataSource(
    list(bcd.topdir.rglob("**/CHECKOUT_26_06_2*/**/raw_level0/**/HOUSEKEEPING*"))
)

lines=ds.messages(fmt=True, supress=False)
Path("messages.txt").open("w").write("\n".join(lines) + "\n")

regs = ds.fpga_registers(namecolumn='value')
#%%

if False:
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(12,8), sharex=True)
    axes = np.ravel(axes)
    detbases = [0x700, 0x600, 0x500, 0x400]

    for off,_,color in ((0xd9, axes[-1], 'r'), (0xdb, axes[0], 'b')):
        for detbase,marker,ax in zip(detbases, 'osD*', axes):
            data = regs[detbase + off]
            ax.plot(data['t'], data['value'], linestyle='none', color=color, marker=marker)
    ax.set(ylabel="Counter low 32bits")



if False:
    plt.close('all')
    plt.figure()
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(12,8), sharex=True)
    for detbase, marker,ax in zip(detbases, 'osD*', np.ravel(axes)):
        # Red for row, cyan for column
        for off,color,label in [(0xdf, 'r', 'row'), (0xde,'c', 'column')]:
            i = detbase + off
            data = regs[i]
            print(f"reg {i:x} {bin(data['value'].min())} - {bin(data['value'].max())}")
            stop,start = np.divmod(data['value'], 1024)
            ax.plot(data['t'], start, '-', color=color)
            ax.plot(data['t'], stop, ':', color=color)
            # ax.plot(data['t'], stop + start/65536, '--', color=color)

if False:
    for det in range(4):
        first_col=100 + 10 * det
        last_col=350 + 10 * det
        first_row=150 + 10 * det
        last_row=400 + 10 * det
        print(f"{last_col:b} {first_col:b} {last_row:b} {first_row:b}")

# %%

plt.close('all')

fig, axes = plt.subplots(nrows=13, ncols=1, figsize=(12,12), sharex=True)
all_axes = np.ravel(axes)
detbases = [0x700, 0x600, 0x500, 0x400]
dout = Path("/tmp/fpgaregs")
dout.mkdir(exist_ok=True)
fig.suptitle("e3 Write DRAM address; e4 wrap count;   db local and global pixel count")
for off,axes,color,desc in ((0xe3, all_axes[:4], 'r','dram address'), (0xe2, all_axes[4:8], 'b', 'wrap count'), (0xdb, all_axes[8:], 'b', 'pixel count'), (0xd9, all_axes[8:], 'r', 'pixel_count')):
    for detbase,marker,ax in zip(detbases, 'osD*', axes):
        detnum = detbase + off
        data = regs[detnum]
        with dout.joinpath(f"fpga_x{detnum}.csv").open("wt") as fout:
            print(f'time,detnum', file=fout)
            for row in data:
                print(f"{row['t']},{row['value']}", file=fout)
        if color == 'r':
            data = data[data['value'] != 0]
        ax.plot(data['t'], data['value'], color=color, marker=marker)
        ax.set(ylabel=f"{desc} {detnum:x}")
laps = ds.data('lapped')[0]
laps = laps[laps['t'] > data['t'][0]]
all_axes[12].plot(laps['t'],laps['lapped'])
all_axes[12].set(ylabel="Lapped")
# ax.set(ylabel="Counter low 32bits")
fig.tight_layout()
pass

# %%
laps = ds.data('lapped')

# %%

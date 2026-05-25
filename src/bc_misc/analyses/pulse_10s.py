from tools import *

# THIS OBSERVATION PRODUCES BURSTS OF PHOTONS EVERY 10 s
plt.close('all')
f2=Path('/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_05_21/level1_tmp/')/'ph_2605202022.gz'
d = fits.getdata(f2)
if f2.exists():
    plt.figure()
    plot_photons_isl4(f2.name, f2)
    plt.savefig(f"/tmp/counts_{f2.name}.png")
    plt.figure()
    plt.plot(d['framenum'], d['bias'], '.')
    plt.gca().set(xlabel=f"Frame number in {f2.name}", ylabel="Frame bias (ADU histogram median)")
    plt.savefig(f"/tmp/bias_{f2.name}.png")
else:
    print(f"no {f2}")
c = Counter(d['time'] - d[0]['time'])
print([(float(t % 10 - 3.45), float(t), n) for t,n in sorted(c.most_common(20)) if n > 10])
# %%

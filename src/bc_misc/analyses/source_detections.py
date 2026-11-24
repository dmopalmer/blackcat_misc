"""
On-board source detections
"""

# %%
import matplotlib as mpl

mpl.use("module://ipympl.backend_nbagg")
import bc_misc
from bc_misc.analyses.tools import *
from mylab import *
import blackcat
import blackcat_data as bcd
from collections import Counter

cat = SourceCatalog()
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


sid = datasource.data("source_id")[0]
counts = Counter(sid["source_id"])


def angs(angs_):
    return np.array(
        [float(str.strip(a.decode(), "'")) for a in angs_], dtype=np.float64
    )


# %%
plt.close('all')
outdir = Path('/tmp/found_sources')
outdir.mkdir(exist_ok = True)
cat = SourceCatalog()
hsz = 1.0
sourcenames = []
sources = []
for sourceid, count in counts.most_common(30):
    if count < 5:
        break
    if sourceid == 0:
        continue
    source = cat.by_id(sourceid)
    sources.append(source)
    sourcenames.append(f"{source['source_name']} ({count})")

print(f"Found:\n{'\n'.join(sourcenames)}")


for source in sources:
    sourceid = source['source_id']
    w = np.ravel(np.argwhere(sid['source_id'] == sourceid))
    matches = sid[w]
    ra,dec = [angs(matches[f'ra_dec({i})']) for i in range(2)]
    ra_all,dec_all = [angs(sid[f'ra_dec({i})']) for i in range(2)]
    ra_all_s, dec_all_s = cat.catalog['ra_dec'].T
    ra_s, dec_s = source['ra_dec']
    sec = 1/np.cos(np.deg2rad(source['ra_dec'][1]))
    fig = plt.figure(figsize=(7.5, 10.5))
    ax_top = fig.add_axes([0.2, 0.25, 0.7, 0.7])
    ax_bottom = fig.add_axes(([0.2, 0.1, 0.7, 0.15]))
    fig.suptitle(source['source_name'])
    ax_top.plot(ra_all, dec_all, 'g.')
    ax_top.plot(ra, dec, 'b.')
    ax_top.plot(ra_s, dec_s, 'r+')
    ax_top.set(title=source['source_name'], xlim=[ra_s + sec*hsz, ra_s - sec*hsz], ylim=[dec_s-hsz, dec_s + hsz], aspect = sec)
    ax_bottom.plot(matches['t'], matches['strength']/matches['local_sigma'], 'b.')
    ax_bottom.tick_params(axis="x", labelrotation=45)
    # fig.autofmt_xdate()
    fig.canvas.draw_idle()
    fig.savefig(outdir.joinpath(f'{source['source_name']}_found.png'))
    


# %%

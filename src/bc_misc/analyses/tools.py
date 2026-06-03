from mylab import *
from typing import Optional
from numpy.typing import NDArray
from functools import cache
import mission_planner
import xraysky
import blackcat_data as bcd
from collections import Counter
from astropy.coordinates import SkyCoord

top_dir = Path("/Volumes/data_x8/blackcat/soc_archive_mirror")

BC_EPOCH_SECS = (datetime(2025,1,1) - datetime(1970,1,1)).total_seconds()

@cache
def new_outdir():
    result = Path(f"/tmp/bc_image_check_plots/plots_{datetime.now(tz=utc):%Y%m%d_%H%M}")
    result.mkdir(parents=True)
    return result


def describe(table):
    for colname in table.names:
        c = table[colname]
        print(
            f"{colname:20}  {np.min(c):20} - {np.max(c):20} {np.mean(c):20} ± {np.std(c):20}"
        )


class Headers:
    def __init__(self, filelist, ext=1):
        self.filelist = []
        self.headers = []
        for f in sorted(filelist):
            try:
                self.headers.append(fits.getheader(f, ext=ext))
                self.filelist.append(f)
            except Exception as e:
                print(f"Skipping {f.name} : {e}")

    def value(self, key, start_diff=False, as_time=False):
        values = np.asarray([header.get(key, np.nan) for header in self.headers])
        if as_time:
            print(type(values[0]))
            values = np.array(
                [np.datetime64(int(value * 1e6), "us") for value in values]
            )
        if start_diff:
            values = values - values[0]
        return values

    def __getitem__(self, key):
        if isinstance(key, str):
            values = self.value(key)
        else:
            values = self.headers[key]
        return values


def where_time_good(d):
    # After 2026 Jan 1.0
    w = d["time"] > 1767250800
    if np.count_nonzero(w):
        return d[w]
    else:
        print("All data timestamped before 2025")
        return d


def get_ph_files(name: str, level: Optional[int] = None, no_x=False) -> list[Path]:
    result = []
    if level is None or level == 0:
        result += sorted(top_dir.rglob(f"raw_level0/ph_{name}*"))
    if level is None or level == 1:
        result += sorted(top_dir.rglob(f"level1_tmp/ph_{name}*"))
    if no_x:
        result = [fname for fname in result if fname.name[13] != "x"]
    return result


def ts_to_time(values):
    """Timestamp either 1970 or 2025 based converted to numpy array time

    Args:
        values (_type_): _description_

    Returns:
        _type_: _description_
    """
    values = np.atleast_1d(values)
    if len(values) == 0:
        return values
    # Check for BlackCAT Epoch
    if values.max() < 3e8:
        values = values + BC_EPOCH_SECS
    return np.array(
        [np.datetime64(int(value * 1e6), "us") for value in values]
    )


def get_satellite() -> sfapi.EarthSatellite:
    from mission_planner.spacetrack_client import get_bc_tle
    tle = get_bc_tle()
    sat = sfapi.EarthSatellite(tle[-2], tle[-1])
    return sat


def satellite_location(t:datetime, satellite = False):
    """Return satellite geographical position (and optionally satellite)

    Args:
        t (datetime): _description_
        satellite (bool, optional): _description_. Defaults to False.
    """
    sft = xraysky.sftime(t)
    sat = get_satellite().at(sft)
    position = sat.subpoint()
    if satellite:
        return position, sat
    else:
        return position

def stack_extension(f, extname:str, mosaic=True) -> NDArray:
    """Read the 4 extver (detectors) of a given extname into an array

    Args:
        f (_type_): _description_
        extension (str): _description_
        mosaic (bool, optional): _description_. Defaults to True.

    Returns:
        NDArray: _description_
    """
    v0,v1,v2,v3 = [fits.getdata(f, (extname, i)) for i in range(4)]
    if mosaic:
        # FIXME: do the rotations and put in the right place for origin='lower'
        # This is probably wrong
        result = np.vstack([np.hstack([v1,v3[::-1,:].T]), np.hstack([v2[:,::-1].T,v0[::-1,::-1]])])
    else:
        result = np.stack([v0,v1,v2,v3], axis=0)
    return result

def plot_orientation(ordata, obsname=None, ax=None):
    if ax is None:
        plt.figure()
        ax_ = plt
    else:
        ax_ = ax
    t = ts_to_time(ordata['time'])
    for i,(name, offset, color) in enumerate((('RA', 0, 'r'), ('Dec + 180', 180, 'g'), ('Roll', 0, 'b'))):
        ax_.plot(t, ordata['POINTING'][:,i] + offset, label=name, color=color)
    if ax is None: ax = ax_.gca()
    if obsname is None:
        obsname=f"{t[0].item():%Y-%m-%d %H:%M:%S} - {t[-1].item():%H:%M:%S}"
    ax.set(title=f"{obsname} Orientation", ylim=[0, 360], xlabel = f"Time on {t[0].item():%Y-%m-%d}")
    ax.legend()
    return ax




def plot_detpos(data, name=None, ax=None):
    if ax is None:
        plt.figure()
        ax_ = plt
    else:
        ax_ = ax
    ax_.plot(data['DETX'], data['DETY'], ',')
    if ax is None: ax = ax_.gca()
    if name is None:
        name=f"{ts_to_time(data[0]['time']).item():%Y-%m-%d %H:%M:%S} - {ts_to_time(data[-1]['time']).item():%H:%M:%S}"
    ax.set(aspect=1, title=f"{name} {len(data)} counts in {np.ptp(data['TIME']):.2f} s", xlim=[-0.024, 0.024], ylim=[-0.024, 0.024])
    return ax


def detsplit(data, field='DETID', matches=range(4), apply=None):
    result = []
    for match in matches:
        r = data[data[field] == match]
        if apply is not None:
            r = apply(r)
        result.append(r)
    return result


def plot_photons_isl4(name,fname, fig=None):
    """Plot photons including Island4

    Full detector plane, per/detector energy histogram, position vs time, Island4 vs time

    Args:
        name (_type_): _description_
        fname (_type_): _description_
        fig (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    # plt.close('all')
    # Positions of detectors in dety,detx indices
    detgrid = [[1,1],[0,0],[1,0],[0,1]]
    # def plot_photons(fname, name=None, fig=None):
    d,h = fits.getdata(fname, header=True)
    if fig is None:
        fig = plt.figure(figsize=[12,12])
    fig.suptitle(name)
    gs = plt.GridSpec(4,4,figure=fig)
    # fig.suptitle(f"{name}\n{detsplit(d, apply=len)} total, {detsplit(d[d['ENERGY'] < 19], apply=len)} below 19 keV in {np.ptp(d['TIME']):.1f} s")
    fig.suptitle(f"{name}  ({len(d)} photons in {np.ptp(d['TIME']):.1f} s -> {len(d)/np.ptp(d['TIME']):.1f} cps")
    ax_dph = fig.add_subplot(gs[0:2, 0:2])
    ax_dph.plot(d['detx'], d['dety'], 'k,')
    ax_dph.set(aspect=1, xlim=[-0.024,0.024], ylim=[-0.024, 0.024])
    
    for i,dd in enumerate(detsplit(d)):
        if len(dd) > 0:
            # xave, yave = [np.mean(dd[f'det{c}']) for c in ('x','y')]
            # ax_dph.text(xave, yave, f"SP{i}", color="r", alpha=1)
            axhist = fig.add_subplot(gs[1-detgrid[i][0], 2+detgrid[i][1]])
            axhist.hist(dd['energy'], bins=np.logspace(np.log10(0.5),np.log10(25),num=100),histtype='step')
            axhist.set(xscale='log',yscale='log',xlim=[1,21],ylim=[0.5,1e5])
    ax_ybytime = fig.add_subplot(gs[2,:])
    for dd,c in zip((d[d['island4'] > 2000], d[d['island4'] <= 2000]), 'rb'):
        if len(dd) > 0:
            ax_ybytime.plot(ts_to_time(dd['time']),dd['rawy']+1600*dd['detid'].astype(int), c+',')
    ax_ybytime.set(ylabel="rawY + K*detid")
    ax_I4bytime = fig.add_subplot(gs[3,:], sharex=ax_ybytime)
    ax_I4bytime.plot(ts_to_time(d['time']),d['island4']+2500*d['detid'].astype(int), 'b,')
    ax_I4bytime.set(ylabel="island4 + K*detid")
    fig.tight_layout()
    fig.show()
    return fig, d, h

def ph_rate(data, binsize:float=1) -> tuple[NDArray[np.datetime64], NDArray[float]]:
    """Photon rate from photon data
    
    use `plt.stairs(edges, values,...)` to do histogram-like with edges 1 longer than values

    Args:
        data (_type_): _description_
        binsize (float, optional): _description_. Defaults to 1.
        edges (bool, optional): _description_. Defaults to True.

    Returns:
        tuple[NDArray[np.datetime64], NDArray[np.float]]: _description_
    """
    # photons must be sorted monotonically
    assert not np.any(data['time'][:-1] > data['time'][1:])
    # Start to work in posix seconds
    tlow,thigh = data[[0,-1]]['time']
    dt = (thigh - tlow)
    nbins = int(np.ceil(dt/binsize))
    bins = np.linspace(tlow, thigh, nbins+1)
    counts, tbins = np.histogram(data['time'], bins=bins)
    rates = counts / (bins[1] - bins[0])
    tedges = ts_to_time(tbins)
    return (tedges, rates)
    
def obsdata(obsname: str|Path):
    if isinstance(obsname, Path):
        obsname = str(obsname)
    # foo/path/ph_blah.gz -> blah
    obsname = obsname.split('/')[-1]
    obsname = obsname.split('.')[0]
    obsname = obsname.split('_')[-1]
    result = dict(obsname = obsname)
    for name, pat in (('l1_file', f'level1_tmp/ph_{obsname}*'),
                      ('l0_file', f'raw_level0/ph_{obsname}*'),
                      ('ori_file', f'level1_tmp/or_{obsname}*')):
        try:
            f = next(top_dir.rglob(pat))
            result[name] = f
        except:
            result[name] = None
    return result


def obs_photons(obsname: str|Path, pointing_tolerance=0.25):
    result = obsdata(obsname)
    try:
        photons,ph_header = fits.getdata(result['l1_file'], header=True)
        # Add a 't' column of np.datetime64
        # photons['t'] = ts_to_time(photons['time'])
        result['photons'] = photons
        result['header'] = ph_header
    except:
        pass
    try:
        ori = fits.getdata(result['ori_file'])
        # Add a 't' column of np.datetime64
        # ori['t'] = ts_to_time(ori['time'])
        
        result['orientation'] = ori
        radecroll = ori['pointing']
        for c in (0,2):
            # Remove the wraps so that 359->360->361 for RA and roll
            radecroll[:,c] = np.unwrap(radecroll[:,c], period=360)
        radecroll_med = np.median(radecroll, axis=0)
        # Convert output to ranges [0,360),[-90,+90],[0,360)

        result['radecroll'] = (radecroll_med + [0,180,0]) % 360 - [0,180,0]
        # Offset from wrapped median values, sum of squares in degrees
        offset_sq = np.sum((radecroll - radecroll_med)**2, axis=1)
        # Map of where the pointing is off by more than the tolerance
        off_point_indices = np.ravel(np.argwhere(offset_sq > pointing_tolerance**2))
        zero_lengths = np.diff(off_point_indices, append=len(offset_sq))
        ii_max_zero_length = np.argmax(zero_lengths)
        # First ori point within tolerance
        i_start = off_point_indices[ii_max_zero_length] + 1
        # Number of within-tolerance points in the run
        maxlen = max(np.max(zero_lengths)-1,0)
        if maxlen == 0:
            # guard kicks to exception
            raise RuntimeError("No stable times")
        assert np.all(offset_sq[i_start:i_start+maxlen] < 1.01*pointing_tolerance**2)
        stable_times = ori['time'][[i_start, i_start+maxlen-1]]
        stable_photons = photons[(stable_times[0] < photons['time']) & (photons['time'] < stable_times[1])]
        result['stable_times'] = ts_to_time(stable_times)
        result['stable_photons'] = stable_photons
    except:
        pass
    
    return result
    
    

def plot_stack(data, name=None, cols=['dety','bias','energy','ISLAND4']):
    fig,ax = plt.subplots(len(cols),1,sharex=True, figsize=(8,2.5*len(cols)))
    fig.suptitle(name)
    for col,ax in zip(cols,np.ravel(ax)):
        ax.plot(d['time'],d[col],',')
        ax.set(ylabel=col)
    fig.tight_layout()
    return fig

def loopcheck(data):
    loops = []
    # First assert fail if frame number jumps backwards
    assert not np.any(np.diff(data['framenum']) < 0)
    # Try looking for dupes
    ctr = Counter([(int(d['framenum']), int(d['detid']), int(d['rawx']), int(d['rawy'])) for d in data])
    (vals, count) = ctr.most_common()[0]
    if count > 1:
        print(f"{count} copies of framenum, detid, rawx, rawy = {vals}")
        return True

balance_boxes = np.array([
    [[ 0.00210067,  0.00210067], [ 0.02206067,  0.02206067]],
    [[-0.02206067, -0.02206067],[-0.00210067, -0.00210067]],
    [[-0.02202067,  0.00210067],[-0.00210067,  0.02206067]],
    [[ 0.00214067, -0.02206067],[ 0.02206067, -0.00210067]]], dtype=np.float32)

scox1_coords  = SkyCoord(ra='244.9794552787600d', dec='-15.6402826851500d', frame='icrs')

# %%

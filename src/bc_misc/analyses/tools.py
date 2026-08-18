from mylab import *
from typing import Optional
from numpy.typing import NDArray
from functools import cache
import mission_planner
import xraysky
import blackcat_data as bcd
from collections import Counter
from collections.abc import Sequence
from typing import Tuple
from astropy.coordinates import SkyCoord
from scipy.spatial.transform import Rotation
from blackcat import update_calibration
from astropy.coordinates import SkyCoord
from astropy import units as u

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


class HotMonitor:
    """Accumulate photons, then 
    """
    __RECENT_HOT_FILE = Path(__file__).parent.joinpath('../../../data/gains/hot_pixels_20260630.gz')
    def __init__(self, hotfile: Optional[Path|bool] = True):
        if hotfile is True:
            hotfile = self.__RECENT_HOT_FILE
        if hotfile:
            # hots_from_file = 1 if hot
            # IMPROVEME also use gains/threshold file
            self.hots_from_file = np.array([fits.getdata(hotfile, ('HOTNESS', i)) for i in range(4)])
        else:
            self.hots_from_file = None
        self.hot_counters = np.zeros((4,550,550), dtype=np.int32)
        
    def _id_py_px(self, counts) -> tuple[NDArray, NDArray, NDArray]:
        """detid, pixel_x, pixel_y
        """
        id = counts['DETIC']
        py  = counts['RAWY']
        px  = counts['RAWX']
        
    def add_hots(self, counts):
        id,py,px = self._id_py_px(counts)
        np.add.at(self.hot_counters, (id, py, px), 1)
        
    def suggest_threshold(self) -> int:
        """Suggest a hot pixel threshold
        
        Remove the zeros and the file-hots from the counters
        Take the median of what remains.
        Remove the values less than 10% of the median
        Take the new median
        Threshold at 3 sigma over twice the median
        

        Returns:
            int: Threshold  hot = (count > result)
        """
        h = np.ravel(self.hot_counters)
        if self.hots_from_file is not None:
            h = h[0 == np.ravel(self.hots_from_file)]
        h = h[h != 0]
        med = np.median(h)
        if med > 10:
            med = np.median(h[h > med//10])
        return int(2*med + 3 * np.sqrt(med))
        
        
    def hot_map(self, threshold: Optional[int] = None, ignore_file = False) -> NDArray[bool]:
        if threshold is None:
            threshold = self.suggest_threshold()
        result = self.hot_counter > threshold
        if not ignore_file:
            result = np.logical_and(result, self.hots_from_file)
        return result
    
    def split_counts(self, counts, threshold: Optional[int] = None, ignore_file = False) -> tuple[NDArray, NDArray]:
        """Split counts into those from cool and hot pixels

        Args:
            counts (_type_): _description_
            threshold (Optional[int], optional): _description_. Defaults to None.
            ignore_file (bool, optional): _description_. Defaults to False.

        Returns:
            tuple[NDArray, NDArray]: _description_
        """
        hm = self.hot_map(threshold=threshold, ignore_file = ignore_file)
        id,py,px = self._id_py_px(counts)
        hot = hot_map[id, py, px]
        hots = counts[hot]
        cools = counts[np.logical_not(hot)]
        return cools, hots
        
        
            
        

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

@cache
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

def stable_pointings(datasource:bcd.DataSource, max_change:float=1):
    """Time intervals when pointing is stable 
    
    Return type:
        dtype([('t_start', '<M8[us]), ('t_end', '<M8[us]'), ('point_dec', '<f8'), ('point_ra', '<f8'), ('point_roll', '<f8'))])
    
    Args:
        datasource (bcd.DataSource): _description_
    Returns:
        NDArray[]
    """
    d = datasource.data('point_ra')[0]
    out_dtype = np.dtype([('t_start', '<M8[us]'), ('t_end', '<M8[us]'), 
                        ('point_dec', '<f8'), ('point_ra', '<f8'), ('point_roll', '<f8')])
    # Yields an array of
    # dtype([('t', '<M8[us]'), ('adcs_state', 'S8'), ('num_attitudes', '<i8'), ('num_cmds', '<i8'), ('num_msg_errors', '<i8'), ('num_positions', '<i8'), ('pointing_state', 'S7'), ('point_dec', '<f8'), ('point_ra', '<f8'), ('point_roll', '<f8'), ('time', '<f8')])
    # Use numpy to find all runs where successive pointing changes (addded in quadrature)
    # are below max_change
    if d is None or len(d) == 0:
        return np.empty(0, dtype=out_dtype)
    
    # Calculate successive differences for coordinates
    diff_ra = np.diff(d['point_ra'])
    diff_dec = np.diff(d['point_dec'])
    diff_roll = np.diff(d['point_roll'])
    
    # Calculate change added in quadrature
    quad_change = np.sqrt(diff_ra**2 + diff_dec**2 + diff_roll**2)
    
    # Create boolean mask where changes are within acceptable bounds
    # Note: len(stable_mask) is len(d) - 1
    stable_mask = quad_change < max_change
    
    # Find transitions (pad to capture boundaries correctly)
    padded = np.concatenate(([False], stable_mask, [False]))
    diff_padded = np.diff(padded.astype(int))
    
    # Start indices are where diff is 1, End indices are where diff is -1
    start_indices = np.where(diff_padded == 1)[0]
    end_indices = np.where(diff_padded == -1)[0]
    
    results = []
    for start, end in zip(start_indices, end_indices):
        # Slice the continuous run array segment
        segment = d[start:end + 1]
        
        results.append((
            segment['t'][0],              # t_start
            segment['t'][-1],             # t_end
            np.median(segment['point_dec']), # Average position over stable run
            np.median(segment['point_ra']),
            np.median(segment['point_roll'])
        ))
        
    return np.array(results, dtype=out_dtype)

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

def rate_from_count(ts, counts) -> tuple[NDArray, NDArray]:
    """Rate from cumulative counts

    Args:
        times (_type_): in seconds
        counts (_type_): cumulative
    Returns:
        datetime64, counts_per_second
    """
    rates = np.diff(counts)/np.clip(np.diff(ts), 0.01, None)
    w = ts[:-1] > 1e8
    return ts_to_time(ts[:-1][w]), rates[w]
    

def plot_rate_from_count(times, counts, ax, *args, drawstyle='steps-pre', **kwargs):
    """Plot rate from cumulative counts

    Args:
        times (_type_): in seconds will have ts_to_time applied
        counts (_type_): cumulative
        ax (_type_): axes to plot on
        **kwargs: passed to plot
    """
    rates = np.diff(counts)/np.clip(np.diff(times), 0.01, None)
    w = times[:-1] > 1e8
    # w = w & np.roll(w,1) & np.roll(w,-1)
    ax.plot(ts_to_time(times[:-1][w]), rates[w], *args, drawstyle=drawstyle, **kwargs)

def binspace(low, high, approxwidth):
    """Generate bins approximately the requested width

    Args:
        low (_type_): _description_
        high (_type_): _description_
        approxwidth (_type_): _description_
    returns:
        (bins, binwidth)
    """
    dt = high-low
    if dt == 0:
        return np.array([low, low+approxwidth]), approxwidth
    nbins = max(np.round(dt/approxwidth).astype(int), 1)
    return np.linspace(low, high, nbins+1, endpoint=True), dt/nbins

def plot_photons_and_rates(l1files, datasource, fig=None):
    if fig is None:
        fig = plt.figure(figsize=[12,12])
    gs = plt.GridSpec(nrows=6, ncols=1, figure=fig)
    # height space
    gs.update(
        left=0.08,    # 5% margin from left edge
        right=0.95,   # 5% margin from right edge
        top=0.95,     # 5% margin from top edge
        bottom=0.05,  # 5% margin from bottom edge
        hspace=0.05   # Keep the subplots close together vertically
    )
    ax_ybytime = fig.add_subplot(gs[0,0])
    # A plot for each detector with 0 at the bottom
    ax_det = [fig.add_subplot(gs[4-i, 0], sharex=ax_ybytime) for i in range(4)]
    ax_rings = fig.add_subplot(gs[5,0], sharex=ax_ybytime)
    for fname in l1files:
        d,h = fits.getdata(fname, header=True)
        times = ts_to_time(d['time'])
        ax_ybytime.plot(times, d['rawy']+1600*d['detid'].astype(int), 'b,')
        ax_ybytime.plot(times[[0,-1]], [6400,7000], 'r', linewidth=0.25, linestyle='-')
        bins, binwidth = binspace(d['time'][0], d['time'][-1], 1)
        for detid in range(4):
            dd = d[d['detid'] == detid]
            h, e = np.histogram(dd['time'], bins=bins)
            ax_det[detid].plot(ts_to_time(bins[:-1]), h/binwidth, 'b-', drawstyle='steps-pre')
    
    for detid, ax in enumerate(ax_det):
        for name_,style,label in [('good_events', 'b:', 'good'), ('hot_events', 'r-', 'hot'), 
                                  ('adc_limit_events', 'g-', 'adc_lim'), 
                                  ('fpga_headwrite', 'c-', 'DRAM written'), 
                                  ('fpga_globalthresh', 'm-', 'Global Thresh'), 
                                  ('fpga_pixthresh', 'y-', 'Pixel thresh')]:
            name = f"{name_}({detid})"
            values = datasource.data(name)[0]
            plot_rate_from_count(values['time'], values[name], ax, style, label=label)
        ax.set(ylim=[1,None], ylabel=f"det {detid} rate (Hz)", yscale="log")
        # if detid != 0:
        #     ax.set(xticklabels = [])
    ax_det[-1].legend(ncols=5, loc="lower center", frameon=False, bbox_to_anchor=(0.5,0.75))
    
    ringvalues = datasource.data('good-max_index')[0]
    for name,style in [('good', 'g-'), ('sus', 'b-'), ('bad', 'c:'), ('cr', 'r')]:
        plot_rate_from_count(ringvalues['time'], ringvalues[f'{name}-max_index'], ax_rings, style, label=name)
    ax_rings.set(ylim=[1,None], ylabel="ring buffers rate (Hz)", yscale="log")
    ax_rings.legend(ncols=4,loc="lower center", frameon=False, bbox_to_anchor=(0.5,0.75))
    fig.tight_layout()
    return fig

def plot_geographic_rates(datasource, name='good', sat=None, fig=None, ax=None) -> tuple[plt.Figure, plt.Axes]:
    import cartopy.crs as ccrs
    from xraysky import sftime
    from skyfield.api import wgs84
    if ax is None:
        if fig is None:
            fig = plt.figure(figsize=(13,7))
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.coastlines()
    name_rate = name+"-max_index_rate"
    data = datasource.data(name_rate)[0]
    if sat is None:
        sat = get_satellite()
    sft = sftime(data['t'].astype('object'))
    geocentric = sat.at(sft)
    lats, lons = wgs84.latlon_of(geocentric)
    # 2-6 = 1e2 - 1e4 counts/s
    log_rate = np.log10(np.clip(data[name_rate], 1e2, 1e4))
    marker_size = log_rate*10
    scatter = ax.scatter(lons.degrees, lats.degrees, c=log_rate, s=3, cmap='jet', vmin=2, vmax=4)
    fig.colorbar(scatter, label="Rate (log10)")
    return fig, ax
    
def plot_pointing_info(datasource: bcd.DataSource):
    pointings = stable_pointings(datasource)
    inst = xraysky.BlackCAT()
    tan_corners = inst.fov()[0]
    imager = xraysky.BC_Imager()
    
    

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
            axhist.set(xscale='log',yscale='log',xlim=[0.5,23],ylim=[0.5,1e5])
    ax_ybytime = fig.add_subplot(gs[2,:])

    try:
        for dd,c in zip((d[d['island4'] > 2000], d[d['island4'] <= 2000]), 'rb'):
            if len(dd) > 0:
                ax_ybytime.plot(ts_to_time(dd['time']),dd['rawy']+1600*dd['detid'].astype(int), c+',')
        ax_I4bytime = fig.add_subplot(gs[3,:], sharex=ax_ybytime)
        ax_I4bytime.plot(ts_to_time(d['time']),d['island4']+2500*d['detid'].astype(int), 'b,')
        ax_I4bytime.set(ylabel="island4 + K*detid")
    except KeyError:
        dd = d
        if len(dd) > 0:
            ax_ybytime.plot(ts_to_time(dd['time']),dd['rawy']+1600*dd['detid'].astype(int), 'b,')
    ax_ybytime.set(ylabel="rawY + K*detid")
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
    
def obs_data(obsname: str|Path):
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
    result = obs_data(obsname)
    try:
        photons,ph_header = fits.getdata(result['l1_file'], header=True)
        # Add a 't' column of np.datetime64
        # photons['t'] = ts_to_time(photons['time'])
        photons,ph_header = update_calibration(photons, ph_header)
        result['photons'] = photons
        result['header'] = ph_header
    except:
        pass
    try:
        ori = fits.getdata(result['ori_file'])
        ori,_ = update_calibration(ori)
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
        try:
            ax.plot(data['time'],data[col],',')
            ax.set(ylabel=col)
        except KeyError:
            pass
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


calibration = dict(
    balance_boxes = np.array([
        [[ 0.00210067,  0.00210067], [ 0.02206067,  0.02206067]],
        [[-0.02206067, -0.02206067],[-0.00210067, -0.00210067]],
        [[-0.02202067,  0.00210067],[-0.00210067,  0.02206067]],
        [[ 0.00214067, -0.02206067],[ 0.02206067, -0.00210067]]], dtype=np.float32),

    alignment_corr = np.array([ -1119.8e-6,   -379.9e-6]),

    # How much to add to each detx, dety
    # Adding to detx moves image spot to the right, adding to dety moves image spot down.
    corr_by_det =\
        np.array([
        [    32.9e-6,    189.0e-6],
        [     3.1e-6,   -167.7e-6],
        [  -232.2e-6,     96.2e-6],
        [   195.9e-6,   -117.6e-6],
        ])

)


def manual_calibrate_data(photons:NDArray, radecroll:Sequence[float], calibration = calibration) -> Tuple[NDArray, Sequence[float]]:
    """Data and pointing adjusted for calibration

    Args:
        photons (NDArray): _description_
        radecroll (Sequence[float]): _description_

    Returns:
        Tuple[NDArray, Sequence[float]]: _description_
    """
    photons = photons.copy()
    cbd = calibration['corr_by_det']
    photons['DETX'] += cbd[photons['DETID'], 0]
    photons['DETY'] += cbd[photons['DETID'], 1]
    # Pointing alignment:
    # This is the nominal way to get from (im_tanx, im_tany, +1) to (sc_x, sc_y, sc_z)
    nom_inst2sc = np.array(
        [[0., 0., -1.],
         [0., 1.,  0.],
         [1., 0.,  0.]])
    nom_flength = 0.154
    # Angle offsets
    align_x_rad, align_y_rad = calibration['alignment_corr']/nom_flength
    twist_rad = 0.0
    # Rotate around im_y (to move image in x) then around im_x (to move in y) then around z (to twist)
    # WARNING: I did not derive the signs from first principles
    # For small angles, ignore non-commutation
    rot = Rotation.from_euler('yxz', [align_x_rad, align_y_rad, twist_rad], degrees=False)
    inst2sc_matrix = nom_inst2sc @ rot.as_matrix()
    sc2inst_matrix = inst2sc_matrix.T
    orientation_nom = xraysky.Orientation(radecroll_inst = radecroll)
    orientation_corr = xraysky.Orientation(q_sc = orientation_nom.q_sc, sc2inst_matrix = sc2inst_matrix)
    
    new_radecroll = orientation_corr.radecroll()
    
    return (photons, new_radecroll)

scox1_coords  = SkyCoord(ra='244.9794552787600d', dec='-15.6402826851500d', frame='icrs')

def ltan(sat, hours=True):
    """
    Local time of Ascending Node of a satellite.

    Makes the most sense for an SSO satellite.


    :param sat: EarthSatellite
    :param hours: if True, result is in hours, else degrees
    :return: LTAN in hours or degrees
    """
    # gmst at gmt midnight in hours
    epoch_utc = sat.epoch.utc
    epoch_hours = epoch_utc.hour + epoch_utc.minute / 60 + epoch_utc.second / 3600
    # Good enough approximation
    gmst_midnight = sat.epoch.gmst - epoch_hours
    ltan_h = ((np.rad2deg(sat.model.nodeo) * 24 / 360.0) - gmst_midnight) % 24
    if hours:
        return ltan_h
    else:
        return 15 * ltan_h


# %%

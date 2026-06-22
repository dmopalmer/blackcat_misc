#%%
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')
from tools import *
import blackcat

dir_latest = sorted(top_dir.glob("payload_checkout/CHECKOUT_26*"))[-1]

print(f"Latest data is {dir_latest.name}")
outdir = new_outdir()


if False:
    ph0_list_all = sorted(top_dir.glob("payload_checkout/CHECKOUT*/raw_level0/ph_*.gz"))
    ph1_list_all = sorted(top_dir.glob("payload_checkout/CHECKOUT*/level1*/*"))
else:
    # Just June
    ph0_list_all = sorted(top_dir.glob("payload_checkout/CHECKOUT_26_06*/raw_level0/ph_*.gz"))
    ph1_list_all = sorted(top_dir.glob("payload_checkout/CHECKOUT_26_06*/level1*/*"))
    

headers0 = Headers(ph0_list_all)
headers0[-1]
#%%

"""
TS_TIME =      1.773344518E+09 / Time for TS_PPSTK                              
TS_UTERR=      1.014830638E-02 / Error in Q7S clock                             
TS_PPSTK=           3495442155 / FPGA tick at PPS                               
TS_TSRC =                    1 / Source of TS_TIME. (1 if PPS)                  
TS_RATTK=            100000000 / Tick rate per second                           
TS_FRMTK=            940297114 / FPGA tick at frame                             
TS_FRMNM=               204760 / Number of frame                                
TS_FRMPR=               656540 / Period of frame in ticks                       
SL_START=     -6.310000000E+02 / Selection start time relative to base          
SL_BASE =                    1 / Selection start time base                      
SL_DUR  =      6.010000000E+02 / Selection duration                             
SL_DET  =                   15 / Selected detector mask                         
SL_GRL  =                 8191 / Selected grades: [63..=0]                      
SL_GRH  =                    0 / Selected grades: [99..=64]                     
SL_BUF  =                    7 / Selected event buffers                         
SL_ELOW =                  250 / Selection low energy limit                     
SL_EHIGH=                25000 / Selection high energy limit                    
"""
#%%
# %matplotlib widget
plt.close('all')
times = headers0.value("TS_TIME", as_time=True)

# If the plot slopes up, increase this.  If it slopes down, decrease,
# by 2.778e-12 * (residual rate in ticks/hour)
# For:
# 2026-03-25 -22.37e-6
# 2026-03-25 Has stair-steps up (1 second jumps) between 03-27 12:00 - and 03-29 12:00
# Then a ~30 second (really a lot more given the data gap) betwee
# 2026-03-29 has stair-steps about 350e-6 = 30 seconds/day

rate_err = -22.37e-6
dt = headers0.value("TS_TIME", start_diff=True)
dtick_pps = headers0.value("TS_PPSTK", start_diff=True)
tick_residual = ((dtick_pps - dt * 1e8*(1+rate_err)) + 2**31)% 2**32 -2**31
plt.plot(times, tick_residual, '.-')
plt.gca().set(xlabel=f"Time", ylabel="Residual (Ticks) for Time", title=f"Time match from LDP headers0 at clock rate {rate_err*1e6:+.2f} PPM")
plt.gcf().autofmt_xdate()
plt.tight_layout()
plt.savefig("/tmp/tick_residual.png")


# %%
plt.close('all')
plt.plot(times, headers0['TS_UTERR'].astype(float),'.-')
plt.gcf().autofmt_xdate()
plt.gca().set(xlabel=f"Time", ylabel="TS_UTERR (s)", title=f"TS_UTERR should be PPS vs Q7S clock error {times[0].astype(datetime):%Y-%m-%d}")
plt.tight_layout()
plt.savefig(outdir.joinpath("UTERR.png"))

#%%
plt.close('all')
f_sysinfo = dir_latest.joinpath('raw_level0/HOUSEKEEPING_SYSINFO.fits.gz')
fits.info(f_sysinfo)
regvalues = fits.getdata(f_sysinfo, "FPGA_REGISTER_VALUES_TLM")
# print(regvalues)
d_sysinfo= fits.getdata(f_sysinfo, "SYSINFO_TLM")
d_sysinfo
"""
'BOOT_CHIP, BOOT_COPY, CFE_TIME, DIE_TEMP, FC_TIME, 
FPGA_REGISTERS-ADDRESS(0), FPGA_REGISTERS-ADDRESS(1), FPGA_REGISTERS-ADDRESS(2), FPGA_REGISTERS-ADDRESS(3), FPGA_REGISTERS-ADDRESS(4), FPGA_REGISTERS-ADDRESS(5), FPGA_REGISTERS-ADDRESS(6), FPGA_REGISTERS-ADDRESS(7), 
FPGA_REGISTERS-VALUE(0), FPGA_REGISTERS-VALUE(1), FPGA_REGISTERS-VALUE(2), FPGA_REGISTERS-VALUE(3), FPGA_REGISTERS-VALUE(4), FPGA_REGISTERS-VALUE(5), FPGA_REGISTERS-VALUE(6), FPGA_REGISTERS-VALUE(7), 
FREE_RAM, LNX_TIME, LOAD_15MIN, LOAD_1MIN, LOAD_5MIN, 
LTC_PPS, LTC_TIME, NUM_CMDS, NUM_FC_PRECISE_TIME, NUM_FC_TIME_OTHER, NUM_MSG_ERRORS, 
RTC_TEMP, SEQUENCE, TOTAL_RAM, UPTIME, TIME'
"""
d = where_time_good(d_sysinfo)
t = ts_to_time(d['time'])


t0 = d['TIME'][0]
tick0 = d['LTC_PPS'][0]
tick_delta = np.diff(d['LTC_PPS'])
# Ticks/s to nearest count
tick_rate = int((1+rate_err) * 1e8)
tick_delta_sec = np.round(tick_delta / tick_rate).astype(np.int64)
tick_sec = np.cumsum(tick_delta_sec)
tick_err = d['LTC_PPS'] - tick0
tick_err[1:] -= tick_delta_sec * tick_rate


# tick_err = (np.array((d['TIME'] - t0) * tick_rate - tick0).astype(np.int64) & 0xffff_ffff).astype(np.int32)
plt.figure()

plt.plot(t, d['lnx_time'] - d['fc_time'], '.', label="Linux time - FC time")
plt.plot(t, d['lnx_time'] - d['time'], '.', label="Linux time - CFE time")
plt.plot(t, tick_err/tick_rate, '.', label="CFE time - tick time")

plt.suptitle("Linux time vs FC time")
plt.gca().set(ylabel="Linux - FC time (s)")
plt.gcf().autofmt_xdate()
plt.legend()
plt.savefig(outdir.joinpath("linux_vs_fc_time.png"))

# %%
plt.close('all')
file_ori_0=dir_latest.joinpath('raw_level0/HOUSEKEEPING_ORIENTATION.fits.gz')
d_change=where_time_good(fits.getdata(file_ori_0,1))

imaxtime = np.argmax(d_change['TIME'])
offtime = d_change['TIME-SEC'][imaxtime] - d_change['TIME'][imaxtime]
if abs(offtime) < 1000:
    offtime = 0.0
times_change=ts_to_time(d_change['TIME-SEC'])

print(f"Time offset = {offtime} = {ts_to_time(offtime)[0].astype(datetime):%Y-%m-%dT%H:%M:%S}")

d_pointing=where_time_good(fits.getdata(file_ori_0,2))


times_point = ts_to_time(d_pointing['TIME'] + offtime)
plt.plot(times_point, d_pointing["point_ra"], 'r.', label="RA˚")
plt.plot(times_point, d_pointing["point_dec"]+180, 'g.', label="dec+180˚")
plt.plot(times_point, d_pointing["point_roll"], 'b.', label="roll")

ax1 = plt.gca()

angle_change = np.sqrt(np.sum(
                        [np.diff(d_pointing[f"point_{coo}"])**2 for coo in ('ra', 'dec', 'roll')],
                        axis=0)) * 60
axr=ax1.twinx()
axr.plot(times_point[1:], angle_change, 'c+', label="Delta angle")
axr.set(ylabel="Angle change (arcmin)", ylim=[0,5])     
plt.gcf().autofmt_xdate()

# states=np.concatenate([d_change[0]['STATE_OLD']], d_change['STATE_NEW'], 


for i in range(len(d_change)-1):
    if d_change['STATE_NEW'][i] != 'STABLE' and d_change['STATE_NEW'][i+1] == 'STABLE' :
        ax1.axvspan(times_change[i], times_change[i+1], color='red', alpha=0.05)

# ax0=ax1.twinx()
# for val,char in [('STABLE', '>'), ('SLEWING', '<')]:
#     w = d_change['STATE_NEW'] == val
#     ax0.plot(times_change[w], d_change['TIME'][w], 'c'+char, label="Transitions")
# ax0.gcf().autofmt_xdate()
ax1.set(xlabel=f"Time based on ORIENTATION_STATE_CHANGE", ylabel="Angles", 
        title=f"Orientation vs time  {times_change[0].astype(datetime):%Y-%m-%d}")

ax1.legend()

# plt.gcf().autofmt_xdate()
plt.tight_layout()
plt.savefig(outdir.joinpath("Pointing.png"))
plt.figure()
plt.plot(times_change, d_change['TIME-SEC'] - d_change['TIME'] - offtime, ".")
plt.gcf().autofmt_xdate()
plt.gca().set(title=f"Orientation change time - CFS clock - offset\nfor offset {offtime:.2f} = {ts_to_time(offtime)[0].astype(datetime):%Y-%m-%dT%H:%M:%S}")


# %%
plt.close('all')
photfiles = sorted(dir_latest.glob('level*1*/ph*.gz'))
photfiles0 = sorted(dir_latest.glob('*level0/ph*.gz*'))


if len(photfiles) == 0 and len(photfiles0) != 0:
    photfiles = []
    for infile in photfiles0:
        obsid = '20'+infile.name[3:13]
        d,h = fits.getdata(infile, header=True)
        if len(d) == 0:
            continue
        try:
            obsid = '20'+infile.name[3:13]
            outfile = outdir.joinpath(f"ph_{obsid}_L1.gz")
            result = blackcat.convert_events_level0_to_level1(infile, observation_id=obsid, outfile=outfile, overwrite=True)
            print(obsid, len(result[0]))
            if len(result[0]):
                photfiles.append(result[-1])
        except Exception as e:
            print(obsid, e)
    
for fname in photfiles:
    d,h= fits.getdata(fname, header=True)
    print(f"{fname.name}: {len(d)} photons")
    if len(d) == 0:
        continue
    ax=plot_detpos(d, fname.name)
    plt.savefig(outdir.joinpath(f"{fname.name}_detpos.png"))



# %%
plt.close('all')
file_phomon = dir_latest.joinpath('raw_level0/HOUSEKEEPING_PHOTON_MON.fits.gz')
fits.info(file_phomon)
d_phmon = where_time_good(fits.getdata(file_phomon))
FLAG=0xffffffff
offtime=0
time_phmon = ts_to_time(d_phmon['TIME'] + offtime)
for grade in ['GOOD','SUS','BAD','CR']:
    c = np.array(d_phmon[grade+'-MAX_INDEX'])
    c[c == FLAG] = 0
    rate = np.diff(c) / 60
    plt.plot(time_phmon, c, '.', label=grade)
plt.legend()
plt.gca().set(title="Photon monitor", ylabel="Cumulative photons", xlabel="date")
plt.gcf().autofmt_xdate()
plt.savefig(outdir.joinpath("photon_monitor.png"))


# %%
for f in sorted(dir_latest.glob('raw_level0/*.fits.gz')):
    # print(f)
    for ext in range(1, 100):
        try:
            d,h=fits.getdata(f, ext, header=True)
        except Exception as e:
            # print(e)
            break
        try:
            t = d['TIME']
            print(f"{f.stem:40} {ext:3} {len(t):5} {t.min():20.0f} {t.max():20.0f}")
        except Exception as e:
            # print(e)
            continue
# %%
plt.close('all')

f='/Volumes/data_x8/blackcat/soc_archive_mirror/payload_checkout/CHECKOUT_26_03_19/raw_level0/HOUSEKEEPING_HEALTH.fits.gz'
fits.info(f)
d=fits.getdata(f,1)
plt.plot(d['time'],'.')
plt.savefig('/tmp/hk_time.png')


# %%

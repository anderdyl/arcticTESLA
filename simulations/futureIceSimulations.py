import xarray as xr
from scipy.io.matlab.mio5_params import mat_struct
import scipy.io as sio
from datetime import datetime, timedelta, date
import numpy as np
from time_operations import xds2datetime as x2d
from time_operations import xds_reindex_daily as xr_daily
from time_operations import xds_common_dates_daily as xcd_daily
import pickle
from dateutil.relativedelta import relativedelta
from alr import ALR_WRP
import sys
import subprocess


if 'darwin' in sys.platform:
    print('Running \'caffeinate\' on MacOSX to prevent the system from sleeping')
    subprocess.Popen('caffeinate')

def ReadMatfile(p_mfile):
    'Parse .mat file to nested python dictionaries'

    def RecursiveMatExplorer(mstruct_data):
        # Recursive function to extrat mat_struct nested contents

        if isinstance(mstruct_data, mat_struct):
            # mstruct_data is a matlab structure object, go deeper
            d_rc = {}
            for fn in mstruct_data._fieldnames:
                d_rc[fn] = RecursiveMatExplorer(getattr(mstruct_data, fn))
            return d_rc

        else:
            # mstruct_data is a numpy.ndarray, return value
            return mstruct_data

    # base matlab data will be in a dict
    mdata = sio.loadmat(p_mfile, squeeze_me=True, struct_as_record=False)
    mdata_keys = [x for x in mdata.keys() if x not in
                  ['__header__','__version__','__globals__']]

    # use recursive function
    dout = {}
    for k in mdata_keys:
        dout[k] = RecursiveMatExplorer(mdata[k])
    return dout

def dateDay2datetimeDate(d_vec):
   '''
   Returns datetime list from a datevec matrix
   d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
   '''
   return [date(int(d[0]), int(d[1]), int(d[2])) for d in d_vec]

def dateDay2datetime(d_vec):
   '''
   Returns datetime list from a datevec matrix
   d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
   '''
   return [date(int(d[0]), int(d[1]), int(d[2])) for d in d_vec]




import pickle
# with open(r"iceData25ClustersMediumishMin150withDayNumRG.pickle", "rb") as input_file:
# with open(r"iceData36ClustersDayNumRGAdjusted.pickle", "rb") as input_file:
with open(r"iceData36ClustersDayNumRGFinalized.pickle", "rb") as input_file:
    # with open(r"iceData36ClustersMin60.pickle", "rb") as input_file:
    iceDWTs = pickle.load(input_file)

iceTime = iceDWTs['dayTime']
#iceTime = [np.array((iceDWTs['year'][i],iceDWTs['month'][i],iceDWTs['day'][i])) for i in range(len(iceDWTs['year']))]
iceDateTimes = iceTime# dateDay2datetime(iceTime)
# bmus = iceDWTs['bmus_corrected']
bmus = iceDWTs['bmus_final']
bmus = bmus[151:]+1
bmus_dates = iceDateTimes[151:]


xds_KMA_fit = xr.Dataset(
    {
        'bmus':(('time',), bmus),
    },
    coords = {'time': bmus_dates}#[datetime(r[0],r[1],r[2]) for r in timeDWTs]
)


# Loading historical Arctic Temperatures for 1979 to 2022
with open(r"predictedArcticTemp.pickle", "rb") as input_file:
   arcticTemperatures = pickle.load(input_file)

timeTemp = arcticTemperatures['futureDate']
arcticTemp = np.array(arcticTemperatures['futureSims'][0])

xds_Temp_fit = xr.Dataset(
    {
        'temp': (('time',), arcticTemp),
    },
    coords = {'time': timeTemp}
)

# reindex to daily data after 1979-01-01 (avoid NaN)
xds_Temp_fit = xr_daily(xds_Temp_fit, datetime(1979, 6, 1),datetime(2022,5,31))

# FUTURE ARCTIC AIR TEMPERATURE
futureArcticTime = arcticTemperatures['futureDate']
futureArcticTemp = arcticTemperatures['futureSims']

# futureArcticTemp[0] = futureArcticTemp[1]

#
# ###### LOADING HISTORICAL MJO AND CHOPPING TO 1979 to 2021
#
# dataMJO = ReadMatfile('/media/dylananderson/Elements/NC_climate/mjo_australia_2021.mat')
#
# yearMonth = np.vstack((dataMJO['year'],dataMJO['month']))
# Dates = np.vstack((yearMonth,dataMJO['day']))
# mjoHistDates = np.array([datetime(r[0],r[1],r[2]) for r in Dates.T])
#
# mjoIndex = np.where((mjoHistDates >= datetime(1979,6,1)) & (mjoHistDates <= datetime(2021,5,31)))
# rmm1Historical = dataMJO['rmm1'][mjoIndex]
# rmm2Historical = dataMJO['rmm2'][mjoIndex]
# mjoDates = mjoHistDates[mjoIndex]
# ### TODO: Change mjoDates to a list?
# xds_MJO_fit = xr.Dataset(
#     {
#         'rmm1': (('time',), rmm1Historical),
#         'rmm2': (('time',), rmm2Historical),
#     },
#     coords = {'time': mjoDates}
# )
# # reindex to daily data after 1979-01-01 (avoid NaN)
# xds_MJO_fit = xr_daily(xds_MJO_fit, datetime(1979, 6, 1),datetime(2021,5,31))
#
#
#
# ###### LOADING AN MJO SIMULATION AND CHOPPING TO the length time of 1979 to 2021
# with open(r"mjoSimulations.pickle", "rb") as input_file:
#    simMJOs = pickle.load(input_file)
# mjoBMUSSim = simMJOs['evbmus_sim']
# mjoRMM1Sim = simMJOs['rmm1Sims']
# mjoRMM2Sim = simMJOs['rmm2Sims']
# mjoDatesSim = simMJOs['dates_sim']

#
# ####### LOADING HISTORICAL SWTS AND CHOPPING TO JUNE TO MAY
#
# with open(r"awtPCs.pickle", "rb") as input_file:
#     historicalAWTs = pickle.load(input_file)
#
# dailyPC1 = historicalAWTs['dailyPC1']
# dailyPC2 = historicalAWTs['dailyPC2']
# dailyPC3 = historicalAWTs['dailyPC3']
# dailyDates = historicalAWTs['dailyDates']
# awt_bmus = historicalAWTs['awt_bmus']
# annualTime = historicalAWTs['annualTime']
# dailyAWT = historicalAWTs['dailyAWT']
#
# # AWT: PCs (Generated with copula simulation. Annual data, parse to daily)
# xds_PCs_fit = xr.Dataset(
#     {
#         'PC1': (('time',), dailyPC1),
#         'PC2': (('time',), dailyPC2),
#         #'PC3': (('time',), dailyPC3),
#     },
#     coords = {'time': [datetime(r[0],r[1],r[2]) for r in dailyDates]}
# )
# # reindex annual data to daily data
# xds_PCs_fit = xr_daily(xds_PCs_fit, datetime(1979,6,1),datetime(2021,5,31))
#
#
#
# ###### LOADING A SIMULATED SWT AND ALIGNING WITH JUNE/MAY and daily values
# # simulated seasonal
# with open(r"awtSimulations.pickle", "rb") as input_file:
#    simSWTs = pickle.load(input_file)
# swtBMUS = simSWTs['evbmus_sim']
# swtPC1 = simSWTs['pc1Sims']
# swtPC2 = simSWTs['pc2Sims']
# swtPC3 = simSWTs['pc3Sims']
# #swtPC4 = simSWTs['pc4Sims']
# swtDatesSim = simSWTs['dates_sim']
#
#
# import matplotlib.pyplot as plt
# plt.figure()
# plt.plot(swtDatesSim,swtPC1[0])
#
#
#
#
# d1 = datetime(2021, 6, 1)
# dt = datetime(2021, 6, 1)
# end = datetime(2121, 6, 2)
# # step = datetime.timedelta(months=1)
# step = relativedelta(days=1)
# simDailyTime = []
# while dt < end:
#     simDailyTime.append(dt)
#     dt += step
# simDailyDatesMatrix = np.array([[r.year,r.month,r.day] for r in simDailyTime])
#



# --------------------------------------
# Mount covariates matrix

# available data:
# model fit: xds_KMA_fit, xds_MJO_fit, xds_PCs_fit
# model sim: xds_MJO_sim, xds_PCs_sim

# covariates: FIT
# d_covars_fit = xcd_daily([xds_PCs_fit, xds_Temp_fit, xds_KMA_fit])
d_covars_fit = xcd_daily([xds_Temp_fit, xds_KMA_fit])

# # # # # PCs covar
# cov_PCs = xds_PCs_fit.sel(time=slice(d_covars_fit[0],d_covars_fit[-1]))
# cov_1 = cov_PCs.PC1.values.reshape(-1,1)
# cov_2 = cov_PCs.PC2.values.reshape(-1,1)
# # cov_3 = cov_PCs.PC3.values.reshape(-1,1)

# MJO covars
cov_Temp = xds_Temp_fit.sel(time=slice(d_covars_fit[0],d_covars_fit[-1]))
cov_1 = cov_Temp.temp.values.reshape(-1,1)

# join covars and norm.
# cov_T = np.hstack((cov_1, cov_2, cov_3))#, cov_4))
#cov_T = np.hstack((cov_1))
cov_T = cov_1
cov_T_mean = np.mean(cov_T,axis=0)
cov_T_std = np.std(cov_T,axis=0)
#cov_T_std = np.array(cov_T_std[0])
# multCovT = np.array([0.31804979/0.31804979, 0.16031134/0.31804979, 0.12182678/0.31804979, 0.09111769/0.31804979, 1, 1])
# multCovT = np.array([0.4019148/0.4019148, 0.11355852/0.4019148, 0.10510168/0.4019148, 1])
multCovT = np.array([1])
# multCovT = np.array([0.4019148/0.4019148, 0.11355852/0.4019148, 1])

covTNorm = np.divide(np.subtract(cov_T,cov_T_mean),cov_T_std)
covTNormalize = np.multiply(covTNorm,multCovT)


i0 = d_covars_fit.index(datetime(int(xds_KMA_fit.time.dt.year[0]),int(xds_KMA_fit.time.dt.month[0]),int(xds_KMA_fit.time.dt.day[0])))

# i0 = d_covars_fit.index(x2d(xds_KMA_fit.time[0]))
cov_KMA = cov_T[i0:,:]
d_covars_fit = d_covars_fit[i0:]

# generate xarray.Dataset
# cov_names = ['PC1', 'PC2', 'Temp']
cov_names = ['Temp']
xds_cov_fit = xr.Dataset(
    {
        'cov_values': (('time','cov_names'), covTNormalize),
    },
    coords = {
        'time': d_covars_fit,
        'cov_names': cov_names,
    }
)


# use bmus inside covariate time frame
xds_bmus_fit = xds_KMA_fit.sel(
    time=slice(d_covars_fit[0], d_covars_fit[-1])
)



# Autoregressive logistic wrapper
num_clusters = 36
sim_num = 10
fit_and_save = True # False for loading
p_test_ALR = '/users/dylananderson/documents/data/pointHope/testIceALR/'

# ALR terms
d_terms_settings = {
    'mk_order'  : 2,
    'constant' : False,
    'long_term' : True,
    'seasonality': (True, [2]),
    'covariates': (True, xds_cov_fit),
}

print('ALR model fit   : {0} --- {1}'.format(
    d_covars_fit[0], d_covars_fit[-1]))
# Autoregressive logistic wrapper
ALRW = ALR_WRP(p_test_ALR)
ALRW.SetFitData(
    num_clusters, xds_bmus_fit, d_terms_settings)

ALRW.FitModel(max_iter=10000)


sim_num = 10
diffSims = 100
evbmus_sim = np.nan*np.ones((len(futureArcticTime)-1,1000))
c = 0
for simIndex in range(diffSims):

    print('working on large-scale climate simulation #{}'.format(simIndex))

    simTemp = np.array(arcticTemperatures['futureSims'][simIndex])

    xds_Temp_sim = xr.Dataset(
        {
            'temp': (('time',), simTemp),
        },
        coords={'time': futureArcticTime}
    )

    # reindex to daily data after 1979-01-01 (avoid NaN)
    xds_Temp_sim = xr_daily(xds_Temp_sim, datetime(1979, 6, 2), datetime(2050, 5, 31))
    #
    # swtBMUsim = swtBMUS[simIndex][0:100]#[0:len(awt_bmus)]
    # swtPC1sim = swtPC1[simIndex][0:100]#[0:len(awt_bmus)]
    # swtPC2sim = swtPC2[simIndex][0:100]#[0:len(awt_bmus)]
    # #swtPC3sim = swtPC3[simIndex][0:100]#[0:len(awt_bmus)]
    #
    # # trainingDates = mjoDatesSim#[datetime(r[0],r[1],r[2]) for r in dailyDates]
    # trainingDates = [datetime(r[0],r[1],r[2]) for r in simDailyDatesMatrix]
    # dailyAWTsim = np.ones((len(trainingDates),))
    # dailyPC1sim = np.ones((len(trainingDates),))
    # dailyPC2sim = np.ones((len(trainingDates),))
    # #dailyPC3sim = np.ones((len(trainingDates),))
    #
    # dailyDatesSWTyear = np.array([r[0] for r in simDailyDatesMatrix])
    # dailyDatesSWTmonth = np.array([r[1] for r in simDailyDatesMatrix])
    # dailyDatesSWTday = np.array([r[2] for r in simDailyDatesMatrix])
    # normPC1 = swtPC1sim
    # normPC2 = swtPC2sim
    # #normPC3 = swtPC3sim
    #
    # for i in range(len(swtBMUsim)):
    #     sSeason = np.where((simDailyDatesMatrix[:, 0] == swtDatesSim[i].year) & (
    #                 simDailyDatesMatrix[:, 1] == swtDatesSim[i].month) & (simDailyDatesMatrix[:, 2] == 1))
    #     ssSeason = np.where((simDailyDatesMatrix[:, 0] == swtDatesSim[i].year + 1) & (
    #                 simDailyDatesMatrix[:, 1] == swtDatesSim[i].month) & (simDailyDatesMatrix[:, 2] == 1))
    #
    #     dailyAWTsim[sSeason[0][0]:ssSeason[0][0] + 1] = swtBMUsim[i] * dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]
    #     dailyPC1sim[sSeason[0][0]:ssSeason[0][0] + 1] = normPC1[i] * np.ones(
    #         len(dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]), )
    #     dailyPC2sim[sSeason[0][0]:ssSeason[0][0] + 1] = normPC2[i] * np.ones(
    #         len(dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]), )
    #     #dailyPC3sim[sSeason[0][0]:ssSeason[0][0] + 1] = normPC3[i] * np.ones(
    #         #len(dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]), )
    #
    # xds_PCs_sim = xr.Dataset(
    #     {
    #         'PC1': (('time',), dailyPC1sim),
    #         'PC2': (('time',), dailyPC2sim),
    #         #'PC3': (('time',), dailyPC3sim),
    #     },
    #     # coords = {'time': [datetime(r.year,r.month,r.day) for r in swtDatesSim]}
    #     coords = {'time': [datetime(r[0], r[1], r[2]) for r in simDailyDatesMatrix]}
    # )
    #
    # # reindex annual data to daily data
    # xds_PCs_sim = xr_daily(xds_PCs_sim,datetime(2021,6,2),datetime(2050,5,31))
    # # xds_PCs_sim = xr_daily(xds_PCs_sim,datetime(1979,6,1),datetime(2021,5,31))
    #


    d_covars_sim = xcd_daily([xds_Temp_sim])#,xds_PCs_sim])
    # cov_PCs_sim = xds_PCs_sim.sel(time=slice(d_covars_sim[0],d_covars_sim[-1]))
    # cov_1_sim = cov_PCs_sim.PC1.values.reshape(-1,1)
    # cov_2_sim = cov_PCs_sim.PC2.values.reshape(-1,1)
    # #cov_3_sim = cov_PCs_sim.PC3.values.reshape(-1,1)
    # #cov_4_sim = cov_PCs_sim.PC4.values.reshape(-1,1)
    cov_Temp_sim = xds_Temp_sim.sel(time=slice(d_covars_sim[0],d_covars_sim[-1]))
    cov_1_sim = cov_Temp_sim.temp.values.reshape(-1,1)
    #cov_6_sim = cov_MJOs_sim.rmm2.values.reshape(-1,1)
    # cov_T_sim = np.hstack((cov_1_sim, cov_2_sim, cov_3_sim))#, cov_5_sim, cov_6_sim))
    cov_T_sim = cov_1_sim#, cov_5_sim, cov_6_sim))

    covTSimNorm = np.divide(np.subtract(cov_T_sim,np.mean(cov_T_sim,axis=0)),np.std(cov_T_sim,axis=0))
    # covTSimNorm = np.divide(np.subtract(cov_T_sim,cov_T_mean),cov_T_std)
    covTSimNormalize = np.multiply(covTSimNorm,multCovT)


    # generate xarray.Dataset
    xds_cov_sim = xr.Dataset(
        {
            'cov_values': (('time','cov_names'), covTSimNormalize),
        },
        coords = {
            'time': d_covars_sim,
            'cov_names': cov_names,
        }
    )



    # --------------------------------------
    # Autoregressive Logistic Regression

    # available data:
    # model fit: xds_KMA_fit, xds_cov_sim, num_clusters
    # model sim: xds_cov_sim, sim_num, sim_years


    #p_report = op.join(p_data, 'r_{0}'.format(name_test))

    #ALRW.Report_Fit() #'/media/dylananderson/Elements/NC_climate/testALR/r_Test', terms_fit==False)


    # # ALR model simulations
    # sim_years = 100
    # # start simulation at PCs available data
    # d1 = x2d(xds_cov_sim.time[0])
    # d2 = datetime(d1.year+sim_years, d1.month, d1.day)
    # dates_sim = [d1 + timedelta(days=i) for i in range((d2-d1).days+1)]
    # dates_sim = dates_sim[0:-2]
    # ALR model simulations
    sim_years = 71
    # start simulation at PCs available data
    # d1 = x2d(xds_cov_sim.time[0])
    d1 = datetime(int(xds_cov_sim.time.dt.year[0]), int(xds_cov_sim.time.dt.month[0]), int(xds_cov_sim.time.dt.day[0]))

    d2 = datetime(d1.year+sim_years, d1.month, d1.day)
    dates_sim = [d1 + timedelta(days=i) for i in range((d2-d1).days+1)]
    dates_sim = dates_sim[0:-2]
    # print some info

    # print('ALR model sim   : {0} --- {1}'.format(
    #     dates_sim[0], dates_sim[-1]))
    print('ALR model sim   : {0} --- {1}'.format(
        d_covars_fit[0], d_covars_fit[-1]))

    # launch simulation
    xds_ALR = ALRW.Simulate(
        sim_num, dates_sim, xds_cov_sim)

    # dates_sim = dates_sim[0:-2]


    # Save results for matlab plot
    evbmus_sim[:,c:c+10] = xds_ALR.evbmus_sims.values
    c = c + 10
    # evbmus_probcum = xds_ALR.evbmus_probcum.values




# p_mat_output = ('/media/dylananderson/Elements/NC_climate/testALR/testFutureAnnual_y{0}s{1}.h5'.format(
#         sim_years, sim_num))
# import h5py
# with h5py.File(p_mat_output, 'w') as hf:
#     hf['bmusim'] = evbmus_sim
#     # hf['probcum'] = evbmus_probcum
#     hf['dates'] = np.vstack(
#         ([d.year for d in dates_sim],
#         [d.month for d in dates_sim],
#         [d.day for d in dates_sim])).T
#

samplesPickle = 'ice36FinalFutureSimulations1000.pickle'
outputSamples = {}
outputSamples['evbmus_sim'] = evbmus_sim
# outputSamples['evbmus_probcum'] = evbmus_probcum
outputSamples['sim_years'] = sim_num
outputSamples['dates_sim'] = dates_sim
# outputSamples['awtBMUsim'] = swtBMUS
# outputSamples['awtPC1sim'] = swtPC1
# outputSamples['awtPC2sim'] = swtPC2
outputSamples['futureTemp'] = futureArcticTemp
outputSamples['futureArcticTime'] = futureArcticTime

# outputSamples['awtPC3sim'] = swtPC3
# outputSamples['mjoRMM1Sim'] = mjoRMM1Sim
# outputSamples['mjoRMM2Sim'] = mjoRMM2Sim

with open(samplesPickle,'wb') as f:
    pickle.dump(outputSamples, f)







def GenOneYearDaily(yy=1981, month_ini=1):
   'returns one generic year in a list of datetimes. Daily resolution'

   dp1 = datetime(yy, month_ini, 2)
   dp2 = dp1 + timedelta(days=364)

   return [dp1 + timedelta(days=i) for i in range((dp2 - dp1).days)]


def GenOneSeasonDaily(yy=1981, month_ini=1):
   'returns one generic year in a list of datetimes. Daily resolution'

   dp1 = datetime(yy, month_ini, 1)
   dp2 = dp1 + timedelta(3*365/12)

   return [dp1 + timedelta(days=i) for i in range((dp2 - dp1).days)]



import matplotlib.pyplot as plt



bmus_dates_months = np.array([d.month for d in dates_sim])
bmus_dates_days = np.array([d.day for d in dates_sim])


# generate perpetual year list
list_pyear = GenOneYearDaily(month_ini=6)
m_plot = np.zeros((num_clusters, len(list_pyear))) * np.nan
num_sim=1
# sort data
for i, dpy in enumerate(list_pyear):
   _, s = np.where(
      [(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)]
   )
   b = evbmus_sim[s,:]
   # b = bmus[s]
   b = b.flatten()

   for j in range(num_clusters):
      _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!
      # _, bb = np.where([(j == b)])  # j+1 starts at 1 bmus value!

      m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)

import matplotlib.cm as cm
dwtcolors = cm.rainbow(np.linspace(0, 1, num_clusters))
# tccolors = np.flipud(cm.autumn(np.linspace(0,1,21)))
# dwtcolors = np.vstack((etcolors,tccolors[1:,:]))




fig = plt.figure(figsize=(10,4))
ax = plt.subplot2grid((1,1),(0,0))
# plot stacked bars
bottom_val = np.zeros(m_plot[1, :].shape)
for r in range(num_clusters):
   row_val = m_plot[r, :]
   ax.bar(list_pyear, row_val, bottom=bottom_val,width=1, color=np.array([dwtcolors[r]]))
   # store bottom
   bottom_val += row_val

import matplotlib.dates as mdates
# customize  axis
months = mdates.MonthLocator()
monthsFmt = mdates.DateFormatter('%b')
ax.set_xlim(list_pyear[0], list_pyear[-1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.set_ylim(0, 1000)
ax.set_ylabel('Probability')




# generate perpetual year list
list_pyear = GenOneYearDaily(month_ini=6)
m_plot = np.zeros((num_clusters, len(list_pyear))) * np.nan
num_sim=1
# sort data
for i, dpy in enumerate(list_pyear):
   _, s = np.where(
      [(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)]
   )
   b = evbmus_sim[s,2]
   # b = bmus[s]
   b = b.flatten()

   for j in range(num_clusters):
      _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!
      # _, bb = np.where([(j == b)])  # j+1 starts at 1 bmus value!

      m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)

import matplotlib.cm as cm
# etcolors = cm.viridis(np.linspace(0, 1, 70-20))
# tccolors = np.flipud(cm.autumn(np.linspace(0,1,21)))
# dwtcolors = np.vstack((etcolors,tccolors[1:,:]))

fig = plt.figure(figsize=(10,4))
ax = plt.subplot2grid((1,1),(0,0))
# plot stacked bars
bottom_val = np.zeros(m_plot[1, :].shape)
for r in range(num_clusters):
   row_val = m_plot[r, :]
   ax.bar(list_pyear, row_val, bottom=bottom_val,width=1, color=np.array([dwtcolors[r]]))
   # store bottom
   bottom_val += row_val

import matplotlib.dates as mdates
# customize  axis
months = mdates.MonthLocator()
monthsFmt = mdates.DateFormatter('%b')
ax.set_xlim(list_pyear[0], list_pyear[-1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
# ax.set_ylim(0, 100)
ax.set_ylim(0, 1)
ax.set_ylabel('Probability')



bmus_dates_months = np.array([d.month for d in iceDateTimes[151:-214]])
bmus_dates_days = np.array([d.day for d in iceDateTimes[151:-214]])


# generate perpetual year list
list_pyear = GenOneYearDaily(month_ini=6)
m_plot = np.zeros((num_clusters, len(list_pyear))) * np.nan
num_sim=1
# sort data
for i, dpy in enumerate(list_pyear):
   _, s = np.where([(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)])
   # b = evbmus_sim[s,:]
   b = bmus[s]
   b = b.flatten()

   for j in range(num_clusters):
      _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!
      # _, bb = np.where([(j == b)])  # j+1 starts at 1 bmus value!

      m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)

# plt.style.use('dark_background')
fig = plt.figure(figsize=(10,4))
ax = plt.subplot2grid((1,1),(0,0))
# plot stacked bars
bottom_val = np.zeros(m_plot[1, :].shape)
for r in range(num_clusters):
   row_val = m_plot[r, :]
   ax.bar(
      list_pyear, row_val, bottom=bottom_val,
      width=1, color=np.array([dwtcolors[r]]))

   # store bottom
   bottom_val += row_val

import matplotlib.dates as mdates

# customize  axis
months = mdates.MonthLocator()
monthsFmt = mdates.DateFormatter('%b')

ax.set_xlim(list_pyear[0], list_pyear[-1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.set_ylim(0, 1)
ax.set_ylabel('Probability')

#
# from matplotlib import gridspec
#
#
# #
# # Lets get complicated...
# # a grid, 8 x 4 for the 8 SWTs and the 4 seasons?
# # generate perpetual seasonal list
# fig = plt.figure()
# gs = gridspec.GridSpec(2, 1, wspace=0.1, hspace=0.15)
#
# # onlyDo = [0,1,2,3,4,5]
# onlyDo = [0,5]
#
# c = 0
# for m in np.unique(onlyDo):
#
#     list_pSeason = GenOneYearDaily(month_ini=6)
#     m_plot = np.zeros((70, len(list_pSeason))) * np.nan
#     num_clusters=70
#     num_sim=1
#     # Identify on the data that occurs in this awt
#     awtIndex = np.where(dailyAWT[0:-2] == m)
#     subsetMonths = bmus_dates_months[awtIndex]
#     subsetDays = bmus_dates_days[awtIndex]
#     subsetBmus = bmus[awtIndex]
#     #subsetBmus = evbmus_sim[awtIndex[0], :]
#     # sort data
#
#     for i, dpy in enumerate(list_pSeason):
#         _, s = np.where([(subsetMonths == dpy.month) & (subsetDays == dpy.day)])
#         # _, s = np.where([(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)])
#         # b = evbmus_sim[s,:]
#         # b = bmus[s]
#         b = subsetBmus[s]
#         b = b.flatten()
#
#         for j in range(num_clusters):
#             _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!
#
#             m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)
#
#     ax = plt.subplot(gs[c])
#     # plot stacked bars
#     bottom_val = np.zeros(m_plot[1, :].shape)
#     for r in range(num_clusters):
#         row_val = m_plot[r, :]
#         ax.bar(list_pSeason, row_val, bottom=bottom_val,width=1, color=np.array([dwtcolors[r]]))
#
#         # store bottom
#         bottom_val += row_val
#     # customize  axis
#     months = mdates.MonthLocator()
#     monthsFmt = mdates.DateFormatter('%b')
#     ax.set_xlim(list_pSeason[0], list_pSeason[-1])
#     ax.xaxis.set_major_locator(months)
#     ax.xaxis.set_major_formatter(monthsFmt)
#     #ax.set_ylim(0, 100)
#     ax.set_ylabel('')
#     c = c + 1
#

# dailyMWT = dailyMWT[0:-2]
#
# # evbmus_sim = evbmus_sim - 1
# # bmus = bmus + 1
# fig = plt.figure()
# gs = gridspec.GridSpec(9, 4, wspace=0.1, hspace=0.15)
# c = 0
# for awt in np.unique(awt_bmus):
#
#     ind = np.where((dailyMWT == awt))[0][:]
#     timeSubDays = bmus_dates_days[ind]
#     timeSubMonths = bmus_dates_months[ind]
#     a = evbmus_sim[ind,:]
#
#     monthsIni = [3,6,9,12]
#     for m in monthsIni:
#
#         list_pSeason = GenOneSeasonDaily(month_ini=m)
#         m_plot = np.zeros((70, len(list_pSeason))) * np.nan
#         num_clusters=70
#         num_sim=1
#         # sort data
#         for i, dpy in enumerate(list_pSeason):
#             _, s = np.where([(timeSubMonths == dpy.month) & (timeSubDays == dpy.day)])
#             b = a[s,:]
#             # b = bmus[s]
#             b = b.flatten()
#             if len(b) > 0:
#                 for j in range(num_clusters):
#                     _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!
#                     # _, bb = np.where([(j == b)])  # j starts at 0 bmus value!
#
#                     m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)
#
#
#         # indNan = np.where(np.isnan(m_plot))[0][:]
#         # if len(indNan) > 0:
#         #     m_plot[indNan] = np.ones((len(indNan),))
#         #m_plot = m_plot[1:,:]
#         ax = plt.subplot(gs[c])
#         # plot stacked bars
#         bottom_val = np.zeros(m_plot[1, :].shape)
#         for r in range(num_clusters):
#             row_val = m_plot[r, :]
#             indNan = np.where(np.isnan(row_val))[0][:]
#             if len(indNan) > 0:
#                 row_val[indNan] = 0
#             ax.bar(list_pSeason, row_val, bottom=bottom_val,width=1, color=np.array([dwtcolors[r]]))
#
#             # store bottom
#             bottom_val += row_val
#         # customize  axis
#         months = mdates.MonthLocator()
#         monthsFmt = mdates.DateFormatter('%b')
#         ax.set_xlim(list_pSeason[0], list_pSeason[-1])
#         ax.xaxis.set_major_locator(months)
#         ax.xaxis.set_major_formatter(monthsFmt)
#         ax.set_ylim(0, 100)
#         ax.set_ylabel('')
#         c = c + 1
#

# sim_num=100

# Lets make a plot comparing probabilities in sim vs. historical
probH = np.nan*np.ones((num_clusters,))
probS = np.nan * np.ones((1000,num_clusters))

for h in np.unique(bmus):
    findH = np.where((bmus == h))[0][:]
    probH[int(h-1)] = len(findH)/len(bmus)
    # probH[int(h)] = len(findH)/len(bmus)

    for s in range(1000):
        findS = np.where((evbmus_sim[0:15706,s] == h))[0][:]
        probS[s,int(h-1)] = len(findS)/len(bmus)
        # probS[s,int(h)] = len(findS)/len(bmus)



plt.figure()
# plt.plot(probH,np.mean(probS,axis=0),'.')
# plt.plot([0,0.03],[0,0.03],'.--')
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
for i in range(num_clusters):
    temp = probS[:,i]
    temp2 = probH[i]
    box1 = ax.boxplot(temp,positions=[temp2],widths=.0006,notch=True,patch_artist=True,showfliers=False)
    plt.setp(box1['boxes'],color=dwtcolors[i])
    plt.setp(box1['means'],color=dwtcolors[i])
    plt.setp(box1['fliers'],color=dwtcolors[i])
    plt.setp(box1['whiskers'],color=dwtcolors[i])
    plt.setp(box1['caps'],color=dwtcolors[i])
    plt.setp(box1['medians'],color=dwtcolors[i],linewidth=0)

    #box1['boxes'].set(facecolor=dwtcolors[i])
    #plt.set(box1['fliers'],markeredgecolor=dwtcolors[i])
ax.plot([0,0.05],[0,0.05],'k.--', zorder=10)
plt.xlim([0,0.15])
plt.ylim([0,0.15])
plt.xticks([0,0.01,0.02,0.03, 0.04,0.05], ['0','0.01','0.02','0.03','0.04','0.05'])
plt.xlabel('Historical Probability')
plt.ylabel('Simulated Probability')
plt.title('Validation of ALR DWT Simulations')






timeAsArray = np.array(iceDateTimes[151:])
plt.figure()
ax1 = plt.subplot2grid((2,1),(0,0))
for qq in range(len(np.unique(bmus))):
    getBMUS = np.where((bmus == qq+1))
    temp = bmus[getBMUS]
    tempTime = timeAsArray[getBMUS]
    ax1.plot(np.array(iceDateTimes[151:])[getBMUS[0]],qq*np.ones((len(temp),)),'.',color=dwtcolors[qq])

ax1.set_xlim([futureArcticTime[0], futureArcticTime[-1]])
timeAsArray = np.array(futureArcticTime)

simIce = 105
ax2 = plt.subplot2grid((2,1),(1,0))
for qq in range(len(np.unique(bmus))):
    getBMUS = np.where((evbmus_sim[:,simIce] == qq+1))
    temp = evbmus_sim[getBMUS,simIce]
    tempTime = timeAsArray[getBMUS]
    ax2.plot(np.array(futureArcticTime)[getBMUS[0]],qq*np.ones((len(temp[0]),)),'.',color=dwtcolors[qq])
ax2.set_xlim([futureArcticTime[0], futureArcticTime[-1]])











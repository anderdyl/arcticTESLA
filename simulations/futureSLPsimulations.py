import xarray as xr
from scipy.io.matlab.mio5_params import mat_struct
import scipy.io as sio
from datetime import datetime, timedelta, date
import numpy as np
from functions.time_operations import xds2datetime as x2d
from functions.time_operations import xds_reindex_daily as xr_daily
from functions.time_operations import xds_common_dates_daily as xcd_daily
import pickle
from dateutil.relativedelta import relativedelta
from functions.alr import ALR_WRP


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
with open(r"dwts.pickle", "rb") as input_file:
    slpDWTs = pickle.load(input_file)


timeDWTs = slpDWTs['SLPtime']
bmus = slpDWTs['bmus_corrected']
bmus_dates = dateDay2datetimeDate(timeDWTs)
bmus_dates_months = np.array([d.month for d in bmus_dates])
bmus_dates_days = np.array([d.day for d in bmus_dates])

bmus = bmus[120:]+1
timeDWTs = timeDWTs[120:]
bmus_dates = bmus_dates[120:]

xds_KMA_fit = xr.Dataset(
    {
        'bmus':(('time',), bmus),
    },
    coords = {'time': [datetime(int(r[0]),int(r[1]),int(r[2])) for r in timeDWTs]}
)



# Loading historical Arctic Temperatures for 1979 to 2022
with open(r"predictedArcticTemps.pickle", "rb") as input_file:
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
xds_Temp_fit = xr_daily(xds_Temp_fit, datetime(1979, 6, 1),datetime(2021,5,31))

# FUTURE ARCTIC AIR TEMPERATURE
futureArcticTime = arcticTemperatures['futureDate']
futureArcticTemp = arcticTemperatures['futureSims']

with open(r"arcticMJO.pickle", "rb") as input_file:
   mjoData = pickle.load(input_file)


timeMJO = mjoData['mjoTime']

bmusMJO = mjoData['bmus']
bmus_datesMJO = timeMJO#dateDay2datetimeDate(timeMJO)
bmus_dates_yearsMJO = np.array([d.year for d in bmus_datesMJO])
bmus_dates_monthsMJO = np.array([d.month for d in bmus_datesMJO])
bmus_dates_daysMJO = np.array([d.day for d in bmus_datesMJO])

xds_MJO_fit = xr.Dataset(
    {
        'rmm1': (('time',), mjoData['mjoRmm1']),
        'rmm2': (('time',), mjoData['mjoRmm2']),
    },
    coords = {'time': [datetime(bmus_dates_yearsMJO[r],bmus_dates_monthsMJO[r],bmus_dates_daysMJO[r]) for r in range(len(bmus_dates_daysMJO))]}
    #coords = {'time': timeMJO}
)
# reindex to daily data after 1979-01-01 (avoid NaN)
xds_MJO_fit = xr_daily(xds_MJO_fit, datetime(1979, 6, 1),datetime(2023,5,31))


#### AWT FROM ENSO SSTs
with open(r"pastAWTs.pickle", "rb") as input_file:
   historicalAWTs = pickle.load(input_file)
awtClusters = historicalAWTs['clusters']
awtPredictor = historicalAWTs['predictor']

awtBmus = awtClusters.bmus.values
pc1Annual = awtClusters.PCs[:,0]
pc2Annual = awtClusters.PCs[:,1]
pc3Annual = awtClusters.PCs[:,2]

awtVariance = awtPredictor['variance'].values
nPercent = awtVariance / np.sum(awtVariance)

dt = datetime(1880, 6, 1)
end = datetime(2023, 6, 1)
#step = datetime.timedelta(months=1)
step = relativedelta(years=1)
sstTime = []
while dt < end:
    sstTime.append(dt)
    dt += step

years = np.arange(1979,2023)
awtYears = np.arange(1880,2023)

awtDailyBmus = np.nan * np.ones(np.shape(bmus))
PC1 = np.nan * np.ones(np.shape(bmus))
PC2 = np.nan * np.ones(np.shape(bmus))
PC3 = np.nan * np.ones(np.shape(bmus))

for hh in years:
   indexDWT = np.where((np.asarray(bmus_dates) >= date(hh,6,1)) & (np.asarray(bmus_dates) <= date(hh+1,5,31)))
   indexAWT = np.where((awtYears == hh))
   awtDailyBmus[indexDWT] = awtBmus[indexAWT]*np.ones(len(indexDWT[0]))
   PC1[indexDWT] = pc1Annual[indexAWT]*np.ones(len(indexDWT[0]))
   PC2[indexDWT] = pc2Annual[indexAWT]*np.ones(len(indexDWT[0]))
   PC3[indexDWT] = pc3Annual[indexAWT]*np.ones(len(indexDWT[0]))

dailyAWT = awtDailyBmus
# AWT: PCs (Generated with copula simulation. Annual data, parse to daily)
xds_PCs_fit = xr.Dataset(
    {
        'PC1': (('time',), PC1),
        'PC2': (('time',), PC2),
        'PC3': (('time',), PC3),
    },
    coords = {'time': [datetime(int(r[0]),int(r[1]),int(r[2])) for r in timeDWTs]}
)
# reindex annual data to daily data
xds_PCs_fit = xr_daily(xds_PCs_fit, datetime(1979,6,1),datetime(2023,5,31))




###### LOADING AN MJO SIMULATION AND CHOPPING TO the length time of 1979 to 2021
with open(r"mjoSimulations.pickle", "rb") as input_file:
   simMJOs = pickle.load(input_file)
mjoBMUSSim = simMJOs['evbmus_sim']
mjoRMM1Sim = simMJOs['rmm1Sims']
mjoRMM2Sim = simMJOs['rmm2Sims']
mjoDatesSim = simMJOs['dates_sim']

###### LOADING A SIMULATED SWT AND ALIGNING WITH JUNE/MAY and daily values
# simulated seasonal
with open(r"futureAWTs.pickle", "rb") as input_file:
   simSWTs = pickle.load(input_file)
swtBMUS = simSWTs['evbmus_sim']
swtPC1 = simSWTs['pc1Sims']
swtPC2 = simSWTs['pc2Sims']
swtPC3 = simSWTs['pc3Sims']
swtDatesSim = simSWTs['dates_sim']


d1 = datetime(2022, 6, 1)
dt = datetime(2022, 6, 1)
end = datetime(2122, 6, 2)
step = relativedelta(days=1)
simDailyTime = []
while dt < end:
    simDailyTime.append(dt)
    dt += step
simDailyDatesMatrix = np.array([[r.year,r.month,r.day] for r in simDailyTime])


d_covars_fit = xcd_daily([xds_MJO_fit, xds_PCs_fit, xds_Temp_fit, xds_KMA_fit])

# PCs covar
cov_PCs = xds_PCs_fit.sel(time=slice(d_covars_fit[0],d_covars_fit[-1]))
cov_1 = np.array(cov_PCs.PC1.values.reshape(-1,1), dtype=np.float64)
cov_2 = np.array(cov_PCs.PC2.values.reshape(-1,1), dtype=np.float64)
cov_3 = np.array(cov_PCs.PC3.values.reshape(-1,1), dtype=np.float64)

# MJO covars
cov_MJO = xds_MJO_fit.sel(time=slice(d_covars_fit[0],d_covars_fit[-1]))
cov_4 = np.array(cov_MJO.rmm1.values.reshape(-1,1), dtype=np.float64)
cov_5 = np.array(cov_MJO.rmm2.values.reshape(-1,1), dtype=np.float64)

cov_Temp = xds_Temp_fit.sel(time=slice(d_covars_fit[0],d_covars_fit[-1]))
cov_6 = np.array(cov_Temp.temp.values.reshape(-1,1), dtype=np.float64)


# join covars and norm.
cov_T = np.hstack((cov_1, cov_2, cov_3, cov_4, cov_5, cov_6))


cov_T_mean = np.mean(cov_T,axis=0)
cov_T_std = np.std(np.array(cov_T, dtype=np.float64),axis=0)
multCovT = np.array([0.45632212/0.45632212, 0.10552158/0.45632212, 0.08360907/0.45632212, 1, 1, 1])

covTNorm = np.divide(np.subtract(cov_T,cov_T_mean),cov_T_std)
covTNormalize = np.multiply(covTNorm,multCovT)

# KMA related covars starting at KMA period
i0 = d_covars_fit.index(x2d(xds_KMA_fit.time[0]))
cov_KMA = cov_T[i0:,:]
d_covars_fit = d_covars_fit[i0:]

# generate xarray.Dataset
cov_names = ['PC1', 'PC2', 'PC3', 'MJO1', 'MJO2', 'TEMP']
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
num_clusters = 49
fit_and_save = True # False for loading
p_test_ALR = '/testALR/'

# ALR terms
d_terms_settings = {
    'mk_order'  : 2,
    'constant' : False,
    'long_term' : True,
    'seasonality': (True, [2, 4]),
    'covariates': (True, xds_cov_fit),
}

print('ALR model fit   : {0} --- {1}'.format(
    d_covars_fit[0], d_covars_fit[-1]))
# Autoregressive logistic wrapper
ALRW = ALR_WRP(p_test_ALR)
ALRW.SetFitData(
    num_clusters, xds_bmus_fit, d_terms_settings)

ALRW.FitModel(max_iter=10000)


diffSims = 50
sim_num = 10
evbmus_sim = np.nan*np.ones((18992,diffSims*sim_num))
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
    xds_Temp_sim = xr_daily(xds_Temp_sim, datetime(2023, 6, 1), datetime(2075, 5, 31))

    rmm1Sim = mjoRMM1Sim[simIndex]
    rmm2Sim = mjoRMM2Sim[simIndex]

    xds_MJOs_sim = xr.Dataset(
        {
            'rmm1': (('time',), rmm1Sim),
            'rmm2': (('time',), rmm2Sim),
        },
        coords = {'time': mjoDatesSim}
    )

    xds_MJOs_sim = xr_daily(xds_MJOs_sim,datetime(2023,6,1),datetime(2075,5,31))

    swtBMUsim = swtBMUS[simIndex][0:100]
    swtPC1sim = swtPC1[simIndex][0:100]
    swtPC2sim = swtPC2[simIndex][0:100]
    swtPC3sim = swtPC3[simIndex][0:100]

    trainingDates = [datetime(r[0],r[1],r[2]) for r in simDailyDatesMatrix]
    dailyAWTsim = np.ones((len(trainingDates),))
    dailyPC1sim = np.ones((len(trainingDates),))
    dailyPC2sim = np.ones((len(trainingDates),))
    dailyPC3sim = np.ones((len(trainingDates),))

    dailyDatesSWTyear = np.array([r[0] for r in simDailyDatesMatrix])
    dailyDatesSWTmonth = np.array([r[1] for r in simDailyDatesMatrix])
    dailyDatesSWTday = np.array([r[2] for r in simDailyDatesMatrix])
    normPC1 = swtPC1sim
    normPC2 = swtPC2sim
    normPC3 = swtPC3sim

    for i in range(len(swtBMUsim)):
        sSeason = np.where((simDailyDatesMatrix[:, 0] == swtDatesSim[i].year) & (
                    simDailyDatesMatrix[:, 1] == swtDatesSim[i].month) & (simDailyDatesMatrix[:, 2] == 1))
        ssSeason = np.where((simDailyDatesMatrix[:, 0] == swtDatesSim[i].year + 1) & (
                    simDailyDatesMatrix[:, 1] == swtDatesSim[i].month) & (simDailyDatesMatrix[:, 2] == 1))

        dailyAWTsim[sSeason[0][0]:ssSeason[0][0] + 1] = swtBMUsim[i] * dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]
        dailyPC1sim[sSeason[0][0]:ssSeason[0][0] + 1] = normPC1[i] * np.ones(
            len(dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]), )
        dailyPC2sim[sSeason[0][0]:ssSeason[0][0] + 1] = normPC2[i] * np.ones(
            len(dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]), )
        dailyPC3sim[sSeason[0][0]:ssSeason[0][0] + 1] = normPC3[i] * np.ones(
            len(dailyAWT[sSeason[0][0]:ssSeason[0][0] + 1]), )

    xds_PCs_sim = xr.Dataset(
        {
            'PC1': (('time',), dailyPC1sim),
            'PC2': (('time',), dailyPC2sim),
            'PC3': (('time',), dailyPC3sim),
        },
        # coords = {'time': [datetime(r.year,r.month,r.day) for r in swtDatesSim]}
        coords = {'time': [datetime(r[0], r[1], r[2]) for r in simDailyDatesMatrix]}
    )

    # reindex annual data to daily data
    xds_PCs_sim = xr_daily(xds_PCs_sim,datetime(2023,6,1),datetime(2075,5,31))
    # xds_PCs_sim = xr_daily(xds_PCs_sim,datetime(1979,6,1),datetime(2021,5,31))



    d_covars_sim = xcd_daily([xds_MJOs_sim,xds_PCs_sim,xds_Temp_sim])
    cov_PCs_sim = xds_PCs_sim.sel(time=slice(d_covars_sim[0],d_covars_sim[-1]))
    cov_1_sim = cov_PCs_sim.PC1.values.reshape(-1,1)
    cov_2_sim = cov_PCs_sim.PC2.values.reshape(-1,1)
    cov_3_sim = cov_PCs_sim.PC3.values.reshape(-1,1)
    #cov_4_sim = cov_PCs_sim.PC4.values.reshape(-1,1)
    cov_MJOs_sim = xds_MJOs_sim.sel(time=slice(d_covars_sim[0],d_covars_sim[-1]))
    cov_4_sim = cov_MJOs_sim.rmm1.values.reshape(-1,1)
    cov_5_sim = cov_MJOs_sim.rmm2.values.reshape(-1,1)

    cov_Temp_sim = xds_Temp_sim.sel(time=slice(d_covars_sim[0], d_covars_sim[-1]))
    cov_6_sim = np.array(cov_Temp_sim.temp.values.reshape(-1, 1), dtype=np.float64)

    cov_T_sim = np.hstack((cov_1_sim, cov_2_sim, cov_3_sim, cov_4_sim, cov_5_sim, cov_6_sim))

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


    # ALR model simulations
    sim_years = 52
    # start simulation at PCs available data
    d1 = x2d(xds_cov_sim.time[0])
    d2 = datetime(d1.year+sim_years, d1.month, d1.day)
    dates_sim = [d1 + timedelta(days=i) for i in range((d2-d1).days+1)]
    dates_sim = dates_sim[0:-2]
    print('ALR model sim   : {0} --- {1}'.format(
        d_covars_fit[0], d_covars_fit[-1]))
    # launch simulation
    xds_ALR = ALRW.Simulate(
        sim_num, dates_sim, xds_cov_sim)
    # Save results for matlab plot
    evbmus_sim[:,c:c+sim_num] = xds_ALR.evbmus_sims.values
    c = c + sim_num




samplesPickle = 'dwtFutureSimulations.pickle'
outputSamples = {}
outputSamples['evbmus_sim'] = evbmus_sim
# outputSamples['evbmus_probcum'] = evbmus_probcum
outputSamples['sim_years'] = sim_num
outputSamples['dates_sim'] = dates_sim
outputSamples['awtBMUsim'] = swtBMUS
outputSamples['awtPC1sim'] = swtPC1
outputSamples['awtPC2sim'] = swtPC2
outputSamples['awtPC3sim'] = swtPC3
outputSamples['mjoRMM1Sim'] = mjoRMM1Sim
outputSamples['mjoRMM2Sim'] = mjoRMM2Sim

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



bmus_dates_months = np.array([d.month for d in mjoDatesSim[0:10226]])
bmus_dates_days = np.array([d.day for d in mjoDatesSim[0:10226]])

num_clusters = 49

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
   b = b.flatten()

   for j in range(num_clusters):
      _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!

      m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)

import matplotlib.cm as cm
dwtcolors = cm.rainbow(np.linspace(0, 1, num_clusters))



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
ax.set_ylim(0, 500)
ax.set_ylabel('Probability')



# generate perpetual year list
list_pyear = GenOneYearDaily(month_ini=6)
m_plot = np.zeros((num_clusters, len(list_pyear))) * np.nan
num_sim=1
# sort data
for i, dpy in enumerate(list_pyear):
   _, s = np.where(
      [(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)])
   b = evbmus_sim[s,2]
   b = b.flatten()

   for j in range(num_clusters):
      _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!

      m_plot[j, i] = float(len(bb) / float(num_sim)) / len(s)


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
ax.set_ylim(0, 1)
ax.set_ylabel('Probability')



bmus_dates_months = np.array([d.month for d in bmus_dates])
bmus_dates_days = np.array([d.day for d in bmus_dates])


# generate perpetual year list
list_pyear = GenOneYearDaily(month_ini=6)
m_plot = np.zeros((num_clusters, len(list_pyear))) * np.nan
num_sim=1
# sort data
for i, dpy in enumerate(list_pyear):
   _, s = np.where([(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)])
   b = bmus[s]
   b = b.flatten()

   for j in range(num_clusters):
      _, bb = np.where([(j + 1 == b)])  # j+1 starts at 1 bmus value!

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



# Lets make a plot comparing probabilities in sim vs. historical
probH = np.nan*np.ones((num_clusters,))
probS = np.nan * np.ones((500,num_clusters))
for h in np.unique(bmus):
    findH = np.where((bmus == h))[0][:]
    probH[int(h-1)] = len(findH)/len(bmus)
    for s in range(500):
        findS = np.where((evbmus_sim[:,s] == h))[0][:]
        probS[s,int(h-1)] = len(findS)/len(evbmus_sim[:,s])



plt.figure()
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
for i in range(num_clusters):
    temp = probS[:,i]
    temp2 = probH[i]
    box1 = ax.boxplot(temp,positions=[temp2],widths=.0005,notch=True,patch_artist=True,showfliers=False)
    plt.setp(box1['boxes'],color=dwtcolors[i])
    plt.setp(box1['means'],color=dwtcolors[i])
    plt.setp(box1['fliers'],color=dwtcolors[i])
    plt.setp(box1['whiskers'],color=dwtcolors[i])
    plt.setp(box1['caps'],color=dwtcolors[i])
    plt.setp(box1['medians'],color=dwtcolors[i],linewidth=0)

ax.plot([0,0.06],[0,0.06],'k.--', zorder=10)
plt.xlim([0,0.06])
plt.ylim([0,0.06])
plt.xticks([0,0.02,0.04,0.06], ['0','0.02','0.04','0.06'])
plt.xlabel('Historical Probability')
plt.ylabel('Simulated Probability')
plt.title('Validation of ALR DWT Simulations')












import os
import numpy as np
import datetime
from netCDF4 import Dataset
from scipy.stats.kde import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import gridspec
import pickle
from scipy.io.matlab.mio5_params import mat_struct
from datetime import datetime, date, timedelta
import random
import itertools
import operator
import scipy.io as sio
import statsmodels.api as sm
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.interpolate import interp1d
from scipy.stats import norm, genpareto, t
from scipy.special import ndtri  # norm inv
import matplotlib.dates as mdates
from scipy.stats import gumbel_l, genextreme
from scipy.spatial import distance
import scipy.io
import h5py
import mat73
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from datetime import timedelta
import itertools
import operator

# define some constants
epoch = dt.datetime(1970, 1, 1)
matlab_to_epoch_days = 719529  # days from 1-1-0000 to 1-1-1970
matlab_to_epoch_seconds = matlab_to_epoch_days * 24 * 60 * 60

def matlab_to_datetime(matlab_date_num_seconds):
    # get number of seconds from epoch
    from_epoch = matlab_date_num_seconds - matlab_to_epoch_seconds

    # convert to python datetime
    return epoch + dt.timedelta(seconds=from_epoch)

def dateDay2datetimeDate(d_vec):
   '''
   Returns datetime list from a datevec matrix
   d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
   '''
   return [dt.date(int(d[0]), int(d[1]), int(d[2])) for d in d_vec]

def dateDay2datetime(d_vec):
   '''
   Returns datetime list from a datevec matrix
   d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
   '''
   return [dt.datetime(int(d[0]), int(d[1]), int(d[2])) for d in d_vec]

import datetime as dt
def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    return dt.datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           + timedelta(hours=int(hours)) \
           + timedelta(minutes=int(minutes)) \
           + timedelta(seconds=round(seconds)) \
           - timedelta(days=366)


def datenum_to_date(datenum):
    """
    Convert Matlab datenum into Python date.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    return dt.datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           - timedelta(days=366)
           #+ timedelta(hours=int(hours)) \
           #+ timedelta(minutes=int(minutes)) \
           #+ timedelta(seconds=round(seconds)) \




import scipy.io as sio

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
   return [datetime(int(d[0]), int(d[1]), int(d[2])) for d in d_vec]





with open(r"dwts.pickle", "rb") as input_file:
    historicalDWTs = pickle.load(input_file)

numClusters = historicalDWTs['num_clusters']
timeDWTs = historicalDWTs['SLPtime']
bmus = historicalDWTs['bmus_corrected']

bmus_dates = dateDay2datetimeDate(timeDWTs)
bmus_dates_months = np.array([d.month for d in bmus_dates])
bmus_dates_days = np.array([d.day for d in bmus_dates])

with open(r"ice.pickle", "rb") as input_file:
    historicalICEs = pickle.load(input_file)

bmusIce = historicalICEs['bmus']
timeIce = historicalICEs['dayTime']
timeArrayIce = np.array(timeIce)
areaBelow = historicalICEs['areaBelow']
bmus_corrected = bmusIce


with open(r"historicalNTR.pickle", "rb") as input_file:
    historicalWLs = pickle.load(input_file)

tNTR = historicalWLs['time']
ntr = historicalWLs['ntr1']


import pickle
with open(r"/Users/dylananderson/Documents/projects/arcticClimate/waves.pickle","rb") as input_file:
    wavesWinds = pickle.load(input_file)

endTime = wavesWinds['endTime']
startTime = wavesWinds['startTime']

import datetime as dt
from dateutil.relativedelta import relativedelta
st = dt.datetime(startTime[0], startTime[1], startTime[2])
end = dt.datetime(endTime[0],endTime[1]+1,1)
step = relativedelta(hours=1)
hourTime = []
while st < end:
    hourTime.append(st)#.strftime('%Y-%m-%d'))
    st += step


beginTime = np.where((np.asarray(hourTime) == datetime(bmus_dates[0].year,bmus_dates[0].month,bmus_dates[0].day,0,0)))
endingTime = np.where((np.asarray(hourTime) == datetime(bmus_dates[-1].year,bmus_dates[-1].month,bmus_dates[-1].day,0,0)))

wh = wavesWinds['metOcean'].Hs[beginTime[0][0]:endingTime[0][0]+24]
tp = wavesWinds['metOcean'].Tp[beginTime[0][0]:endingTime[0][0]+24]
dm = wavesWinds['metOcean'].Dm[beginTime[0][0]:endingTime[0][0]+24]
u10 = wavesWinds['metOcean'].u10[beginTime[0][0]:endingTime[0][0]+24]
v10 = wavesWinds['metOcean'].v10[beginTime[0][0]:endingTime[0][0]+24]
sst = wavesWinds['metOcean'].sst[beginTime[0][0]:endingTime[0][0]+24]
ssr = wavesWinds['metOcean'].ssr[beginTime[0][0]:endingTime[0][0]+24]
t2m = wavesWinds['metOcean'].t2m[beginTime[0][0]:endingTime[0][0]+24]

# Point Hope
waveNorm = dm - 334
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
neg2 = np.where((waveNorm < -180))
waveNorm[neg2[0]] = waveNorm[neg2[0]]+360
dmOG = dm
dm = waveNorm

time_all = np.asarray([datetime.utcfromtimestamp((qq - np.datetime64(0, 's')) / np.timedelta64(1, 's'))for qq in wavesWinds['metOcean'].timeWave])[beginTime[0][0]:endingTime[0][0]+24]


def hourlyVec2datetime(d_vec):
   '''
   Returns datetime list from a datevec matrix
   d_vec = [[y1 m1 d1 H1 M1],[y2 ,2 d2 H2 M2],..]
   '''
   return [datetime(d[0], d[1], d[2], d[3]) for d in d_vec]

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


realWavesPickle = 'realWaves.pickle'

outputReal = {}
outputReal['tWave'] = time_all
outputReal['hsCombined'] = wh
outputReal['tpCombined'] = tp
outputReal['dmCombined'] = dmOG
outputReal['waveNorm'] = waveNorm
outputReal['ntr'] = ntr
outputReal['tNTR'] = tNTR
outputReal['t2m'] = t2m

with open(realWavesPickle, 'wb') as f:
    pickle.dump(outputReal, f)





dt = datetime(bmus_dates[0].year, bmus_dates[0].month, bmus_dates[0].day)
end = datetime(bmus_dates[-1].year,bmus_dates[-1].month+1,1)
step = timedelta(days=1)
midnightTime = []
while dt < end:
    midnightTime.append(dt)#.strftime('%Y-%m-%d'))
    dt += step



grouped = [[e[0] for e in d[1]] for d in itertools.groupby(enumerate(bmus), key=operator.itemgetter(1))]

groupLength = np.asarray([len(i) for i in grouped])
bmuGroup = np.asarray([bmus[i[0]] for i in grouped])
timeGroup = [np.asarray(midnightTime)[i] for i in grouped]

startTimes = [i[0] for i in timeGroup]
endTimes = [i[-1] for i in timeGroup]

hydrosInds = np.unique(bmuGroup)

hydros = list()
c = 0
for p in range(len(np.unique(bmuGroup))):

    tempInd = p
    print('working on bmu = {}'.format(tempInd))
    index = np.where((bmuGroup==tempInd))[0][:]

    print('should have at least {} days in it'.format(len(index)))
    subLength = groupLength[index]
    m = np.ceil(np.sqrt(len(index)))
    tempList = list()
    counter = 0
    for i in range(len(index)):
        st = startTimes[index[i]]
        et = endTimes[index[i]] + timedelta(days=1)

        waveInd = np.where((time_all < et) & (time_all >= st))
        ntrInd = np.where((tNTR < et) & (tNTR >= st))

        if len(waveInd[0]) > 0:
            newTime = startTimes[index[i]]

            tempHydroLength = subLength[i]

            while tempHydroLength > 1:
                # randLength = random.randint(1, 2)
                etNew = newTime + timedelta(days=1)
                # if etNew >= et:
                #     etNew = newTime + timedelta(days=1)
                waveInd = np.where((time_all < etNew) & (time_all >= newTime))
                fetchInd = np.where((np.asarray(timeArrayIce) == newTime))
                ntrInd = np.where((tNTR < etNew) & (tNTR >= newTime))

                deltaDays = et-etNew
                c = c + 1
                counter = counter + 1
                # if counter > 15:
                #     print('we have split 15 times')

                tempDict = dict()
                tempDict['time'] = time_all[waveInd[0]]
                tempDict['numDays'] = subLength[i]
                tempDict['hs'] = wh[waveInd[0]]
                tempDict['tp'] = tp[waveInd[0]]
                tempDict['dm'] = dm[waveInd[0]]
                tempDict['v10'] = v10[waveInd[0]]
                tempDict['u10'] = u10[waveInd[0]]
                tempDict['sst'] = sst[waveInd[0]]
                tempDict['ssr'] = ssr[waveInd[0]]
                tempDict['t2m'] = t2m[waveInd[0]]
                tempDict['ntr'] = ntr[ntrInd[0]]

                tempDict['fetch'] = areaBelow[fetchInd[0][0]]
                tempDict['cop'] = np.asarray([np.nanmin(wh[waveInd[0]]), np.nanmax(wh[waveInd[0]]),
                                              np.nanmin(tp[waveInd[0]]), np.nanmax(tp[waveInd[0]]),
                                              np.nanmean(dm[waveInd[0]]),np.nanmean(u10[waveInd[0]]),
                                              np.nanmean(v10[waveInd[0]]),np.nanmean(sst[waveInd[0]]),
                                              np.nanmean(ssr[waveInd[0]]),np.nanmean(t2m[waveInd[0]]),
                                              areaBelow[fetchInd[0][0]], np.nanmean(ntr[ntrInd[0]]),np.nanmin(t2m[waveInd[0]]),np.nanmax(t2m[waveInd[0]])])

                tempDict['hsMin'] = np.nanmin(wh[waveInd[0]])
                tempDict['hsMax'] = np.nanmax(wh[waveInd[0]])
                tempDict['tpMin'] = np.nanmin(tp[waveInd[0]])
                tempDict['tpMax'] = np.nanmax(tp[waveInd[0]])
                tempDict['dmMean'] = np.nanmean(dm[waveInd[0]])
                tempDict['u10Mean'] = np.nanmean(u10[waveInd[0]])
                tempDict['u10Max'] = np.nanmax(u10[waveInd[0]])
                tempDict['u10Min'] = np.nanmin(u10[waveInd[0]])
                tempDict['v10Max'] = np.nanmax(v10[waveInd[0]])
                tempDict['v10Mean'] = np.nanmean(v10[waveInd[0]])
                tempDict['v10Min'] = np.nanmin(v10[waveInd[0]])
                tempDict['sstMean'] = np.nanmean(sst[waveInd[0]])
                tempDict['ssrMean'] = np.nanmean(ssr[waveInd[0]])
                tempDict['ssrMin'] = np.nanmin(ssr[waveInd[0]])
                tempDict['ssrMax'] = np.nanmax(ssr[waveInd[0]])
                tempDict['t2mMean'] = np.nanmean(t2m[waveInd[0]])
                tempDict['t2mMin'] = np.nanmin(t2m[waveInd[0]])
                tempDict['t2mMax'] = np.nanmax(t2m[waveInd[0]])
                tempDict['ntrMean'] = np.nanmean(ntr[ntrInd[0]])
                tempDict['ntrMin'] = np.nanmin(ntr[ntrInd[0]])
                tempDict['ntrMax'] = np.nanmax(ntr[ntrInd[0]])
                tempList.append(tempDict)
                tempHydroLength = tempHydroLength-1
                newTime = etNew
            else:
                waveInd = np.where((time_all < et) & (time_all >= newTime))
                fetchInd = np.where((np.asarray(timeArrayIce) == newTime))
                ntrInd = np.where((tNTR < et) & (tNTR >= newTime))

                c = c + 1
                counter = counter+1
                tempDict = dict()
                tempDict['time'] = time_all[waveInd[0]]
                tempDict['numDays'] = subLength[i]
                tempDict['hs'] = wh[waveInd[0]]
                tempDict['tp'] = tp[waveInd[0]]
                tempDict['dm'] = dm[waveInd[0]]
                tempDict['v10'] = v10[waveInd[0]]
                tempDict['u10'] = u10[waveInd[0]]
                tempDict['sst'] = sst[waveInd[0]]
                tempDict['ssr'] = ssr[waveInd[0]]
                tempDict['t2m'] = t2m[waveInd[0]]
                tempDict['ntr'] = ntr[ntrInd[0]]

                tempDict['fetch'] = areaBelow[fetchInd[0][0]]
                tempDict['cop'] = np.asarray([np.nanmin(wh[waveInd[0]]), np.nanmax(wh[waveInd[0]]),
                                              np.nanmin(tp[waveInd[0]]), np.nanmax(tp[waveInd[0]]),
                                              np.nanmean(dm[waveInd[0]]),np.nanmean(u10[waveInd[0]]),
                                              np.nanmean(v10[waveInd[0]]),np.nanmean(sst[waveInd[0]]),
                                              np.nanmean(ssr[waveInd[0]]),np.nanmean(t2m[waveInd[0]]),
                                              areaBelow[fetchInd[0][0]], np.nanmean(ntr[ntrInd[0]]),np.nanmin(t2m[waveInd[0]]),np.nanmax(t2m[waveInd[0]])])
                tempDict['hsMin'] = np.nanmin(wh[waveInd[0]])
                tempDict['hsMax'] = np.nanmax(wh[waveInd[0]])
                tempDict['tpMin'] = np.nanmin(tp[waveInd[0]])
                tempDict['tpMax'] = np.nanmax(tp[waveInd[0]])
                tempDict['dmMean'] = np.nanmean(dm[waveInd[0]])
                tempDict['u10Mean'] = np.nanmean(u10[waveInd[0]])
                tempDict['u10Max'] = np.nanmax(u10[waveInd[0]])
                tempDict['u10Min'] = np.nanmin(u10[waveInd[0]])
                tempDict['v10Max'] = np.nanmax(v10[waveInd[0]])
                tempDict['v10Mean'] = np.nanmean(v10[waveInd[0]])
                tempDict['v10Min'] = np.nanmin(v10[waveInd[0]])
                tempDict['sstMean'] = np.nanmean(sst[waveInd[0]])
                tempDict['ssrMean'] = np.nanmean(ssr[waveInd[0]])
                tempDict['ssrMin'] = np.nanmin(ssr[waveInd[0]])
                tempDict['ssrMax'] = np.nanmax(ssr[waveInd[0]])
                tempDict['t2mMean'] = np.nanmean(t2m[waveInd[0]])
                tempDict['t2mMin'] = np.nanmin(t2m[waveInd[0]])
                tempDict['t2mMax'] = np.nanmax(t2m[waveInd[0]])
                tempDict['ntrMean'] = np.nanmean(ntr[ntrInd[0]])
                tempDict['ntrMin'] = np.nanmin(ntr[ntrInd[0]])
                tempDict['ntrMax'] = np.nanmax(ntr[ntrInd[0]])
                tempList.append(tempDict)
    print('we have split {} times in bmu {}'.format(counter,p))
    hydros.append(tempList)




myFmt = mdates.DateFormatter('%d')

copulaData = list()
copulaDataOnlyWaves = list()
for i in range(len(np.unique(bmuGroup))):
    tempHydros = hydros[i]
    dataCop = []
    dataCopOnlyWaves = []
    for kk in range(len(tempHydros)):
        dataCop.append(list([tempHydros[kk]['hsMax'], tempHydros[kk]['hsMin'], tempHydros[kk]['tpMax'],
                             tempHydros[kk]['tpMin'], tempHydros[kk]['dmMean'], tempHydros[kk]['u10Max'],
                             tempHydros[kk]['u10Min'], tempHydros[kk]['v10Max'], tempHydros[kk]['v10Min'],
                              tempHydros[kk]['ssrMean'], tempHydros[kk]['t2mMean'],
                             tempHydros[kk]['fetch'], tempHydros[kk]['ntrMean'], tempHydros[kk]['sstMean'],tempHydros[kk]['t2mMax'],tempHydros[kk]['t2mMin'],len(tempHydros[kk]['time']), kk]))

        if np.isnan(tempHydros[kk]['hsMax'])==1:
            print('no waves here')

        else:
            dataCopOnlyWaves.append(list([tempHydros[kk]['hsMax'],tempHydros[kk]['hsMin'],tempHydros[kk]['tpMax'],
                             tempHydros[kk]['tpMin'], tempHydros[kk]['dmMean'], tempHydros[kk]['u10Max'],
                             tempHydros[kk]['u10Min'], tempHydros[kk]['v10Max'], tempHydros[kk]['v10Min'],
                             tempHydros[kk]['ssrMean'], tempHydros[kk]['t2mMean'],
                             tempHydros[kk]['fetch'], tempHydros[kk]['ntrMean'], tempHydros[kk]['sstMean'],tempHydros[kk]['t2mMax'],tempHydros[kk]['t2mMin'],len(tempHydros[kk]['time']),kk]))

    copulaData.append(dataCop)
    copulaDataOnlyWaves.append(dataCopOnlyWaves)


bmuDataNormalized = list()
bmuDataMin = list()
bmuDataMax = list()
bmuDataStd = list()
for i in range(len(np.unique(bmuGroup))):
    temporCopula = np.asarray(copulaDataOnlyWaves[i])
    if len(temporCopula) == 0:
        bmuDataNormalized.append(np.vstack((0, 0)).T)
        bmuDataMin.append([0, 0])
        bmuDataMax.append([0, 0])
        bmuDataStd.append([0, 0])
    else:
        dataHs = np.array([sub[0] for sub in copulaDataOnlyWaves[i]])
        data = temporCopula[~np.isnan(dataHs)]
        data2 = data[~np.isnan(data[:,0])]
        if len(data2) == 0:
            print('woah, no waves here bub')
        #     data2 = data
        #     data2[:,5] = 0
            bmuDataNormalized.append(np.vstack((0, 0)).T)
            bmuDataMin.append([0, 0])
            bmuDataMax.append([0, 0])
            bmuDataStd.append([0, 0])
        else:
            maxDm = np.nanmax(data2[:,4])
            minDm = np.nanmin(data2[:,4])
            stdDm = np.nanstd(data2[:,4])
            dmNorm = (data2[:,4] - minDm) / (maxDm-minDm)
            maxSs = np.nanmax(data2[:,12])
            minSs = np.nanmin(data2[:,12])
            stdSs = np.nanstd(data2[:,12])
            ssNorm = (data2[:,12] - minSs) / (maxSs-minSs)
            bmuDataNormalized.append(np.vstack((dmNorm,ssNorm)).T)
            bmuDataMin.append([minDm,minSs])
            bmuDataMax.append([maxDm,maxSs])
            bmuDataStd.append([stdDm,stdSs])



normalizedHydros = list()
for i in range(len(np.unique(bmuGroup))):
    tempHydros = hydros[i]
    tempList = list()
    for mm in range(len(tempHydros)):
        if np.isnan(tempHydros[mm]['hsMin']):
            print('no waves')
        else:
            tempDict = dict()
            tempDict['hsNorm'] = (tempHydros[mm]['hs'] - tempHydros[mm]['hsMin']) / (tempHydros[mm]['hsMax']- tempHydros[mm]['hsMin'])
            tempDict['tpNorm'] = (tempHydros[mm]['tp'] - tempHydros[mm]['tpMin']) / (tempHydros[mm]['tpMax']- tempHydros[mm]['tpMin'])
            tempDict['timeNorm'] = np.arange(0,1,1/len(tempHydros[mm]['time']))[0:len(tempDict['hsNorm'])]
            tempDict['dmNorm'] = (tempHydros[mm]['dm']) - tempHydros[mm]['dmMean']
            tempDict['uNorm'] = (tempHydros[mm]['u10'] - tempHydros[mm]['u10Min']) / (
                        tempHydros[mm]['u10Max'] - tempHydros[mm]['u10Min'])
            tempDict['vNorm'] = (tempHydros[mm]['v10'] - tempHydros[mm]['v10Min']) / (
                        tempHydros[mm]['v10Max'] - tempHydros[mm]['v10Min'])
            # tempDict['ntrNorm'] = (tempHydros[mm]['ntr'] - tempHydros[mm]['ntrMin']) / (tempHydros[mm]['ntrMax']- tempHydros[mm]['ntrMin'])
            tempDict['ntrNorm'] = (tempHydros[mm]['ntr'] - tempHydros[mm]['ntrMean'])
            tempDict['t2mNorm'] = (tempHydros[mm]['t2m'] - tempHydros[mm]['t2mMin']) / (tempHydros[mm]['t2mMax']- tempHydros[mm]['t2mMin'])

            tempList.append(tempDict)
    normalizedHydros.append(tempList)




import pickle

normHydrosPickle = 'normalizedWaveHydrographs.pickle'
outputHydrosNorm = {}
outputHydrosNorm['normalizedHydros'] = normalizedHydros
outputHydrosNorm['bmuDataMin'] = bmuDataMin
outputHydrosNorm['bmuDataMax'] = bmuDataMax
outputHydrosNorm['bmuDataStd'] = bmuDataStd
outputHydrosNorm['bmuDataNormalized'] = bmuDataNormalized

with open(normHydrosPickle,'wb') as f:
    pickle.dump(outputHydrosNorm, f)

hydrosPickle = 'waveHydrographs.pickle'

outputHydros = {}
outputHydros['hydros'] = hydros
with open(hydrosPickle,'wb') as f:
    pickle.dump(outputHydros, f)

copPickle = 'hydrographCopulaData.pickle'

outputCopula = {}
outputCopula['copulaData'] = copulaData
outputCopula['copulaDataNoNaNs'] = copulaDataOnlyWaves
with open(copPickle,'wb') as f:
    pickle.dump(outputCopula, f)


historicalPickle = 'historicalData.pickle'

outputHistorical = {}
outputHistorical['grouped'] = grouped
outputHistorical['groupLength'] = groupLength
outputHistorical['bmuGroup'] = bmuGroup
outputHistorical['timeGroup'] = timeGroup

with open(historicalPickle,'wb') as f:
    pickle.dump(outputHistorical, f)













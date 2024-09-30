import scipy.io as sio
import os
import datetime as dt
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from scipy.io.matlab.mio5_params import mat_struct
from mpl_toolkits.basemap import Basemap


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

# define some constants
epoch = datetime(1970, 1, 1)
matlab_to_epoch_days = 719529  # days from 1-1-0000 to 1-1-1970
matlab_to_epoch_seconds = matlab_to_epoch_days * 24 * 60 * 60

def matlab_to_datetime(matlab_date_num_seconds):
    # get number of seconds from epoch
    from_epoch = matlab_date_num_seconds - matlab_to_epoch_seconds

    # convert to python datetime
    return epoch + dt.timedelta(seconds=from_epoch)




data = ReadMatfile('noaaAlaskaTides.mat')

tideTimeEmPurdoe = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['purdoeDailyData']['time']])
tideTimeEmNome = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['nomeDailyData']['time']])
tideTimeEmRed = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['redDailyData']['time']])

mslEmNome = data['nomeDailyData']['B']['msl'][1]*data['nomeDailyData']['time']+data['nomeDailyData']['B']['msl'][0]
mslEmPurdoe = data['purdoeDailyData']['B']['msl'][1]*data['purdoeDailyData']['time']+data['purdoeDailyData']['B']['msl'][0]
mslEmRed = data['redDailyData']['B']['msl'][1]*data['redDailyData']['time']+data['redDailyData']['B']['msl'][0]

tideEmPurdoe = data['purdoe']['tideoutwithPred']
tideEmNome = data['now']['tideoutwithPred']
tideEmRed = data['red']['tideoutwithPred']


tidePredNome = data['nomeDailyData']['tide']
tideWlNome = data['nomeDailyData']['wl']
tideTimeNome = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['nomeDailyData']['time']])
mslPredNome = data['nomeDailyData']['B']['msl'][1]*data['nomeDailyData']['time']+data['nomeDailyData']['B']['msl'][0]
startTime = 0
endTime = -(1488+24*365)
timeDataNome = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['nomeDailyData']['time']])[startTime:endTime]
wlNome = data['nomeDailyData']['wl'][startTime:endTime]
ssNome = data['nomeDailyData']['ss'][startTime:endTime]
seasonalNome = data['nomeDailyData']['seasonal'][startTime:endTime]
mmslaNome = data['nomeDailyData']['mmsla'][startTime:endTime]
dslaNome = data['nomeDailyData']['dsla'][startTime:endTime]
mslNome = data['nomeDailyData']['msl'][startTime:endTime]


tidePredRed = data['redDailyData']['tide']
tideWlRed = data['redDailyData']['wl']
tideTimeRed = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['redDailyData']['time']])
mslPredRed = data['redDailyData']['B']['msl'][1]*data['redDailyData']['time']+data['redDailyData']['B']['msl'][0]
startTime = 5804
endTime = -(831+24*365)
timeDataRed = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['redDailyData']['time']])[startTime:endTime]
wlRed = data['redDailyData']['wl'][startTime:endTime]
ssRed = data['redDailyData']['ss'][startTime:endTime]
seasonalRed = data['redDailyData']['seasonal'][startTime:endTime]
mmslaRed = data['redDailyData']['mmsla'][startTime:endTime]
dslaRed = data['redDailyData']['dsla'][startTime:endTime]
mslRed = data['redDailyData']['msl'][startTime:endTime]



tidePredPurdoe = data['purdoeDailyData']['tide']
tideWlPurdoe = data['purdoeDailyData']['wl']
tideTimePurdoe = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['purdoeDailyData']['time']])
mslPredPurdoe = data['purdoeDailyData']['B']['msl'][1]*data['purdoeDailyData']['time']+data['purdoeDailyData']['B']['msl'][0]
startTime = 0
endTime = -(365*24+1489)
timeDataPurdoe = np.array([matlab_to_datetime(tempTime * 24 * 60 * 60) for tempTime in data['purdoeDailyData']['time']])[startTime:endTime]
wlPurdoe = data['purdoeDailyData']['wl'][startTime:endTime]
ssPurdoe = data['purdoeDailyData']['ss'][startTime:endTime]
seasonalPurdoe = data['purdoeDailyData']['seasonal'][startTime:endTime]
mmslaPurdoe = data['purdoeDailyData']['mmsla'][startTime:endTime]
dslaPurdoe = data['purdoeDailyData']['dsla'][startTime:endTime]
mslPurdoe = data['purdoeDailyData']['msl'][startTime:endTime]




import pickle
with open(r"dwts.pickle", "rb") as input_file:
    dwtInput = pickle.load(input_file)

SLPpcs= dwtInput['PCs']
SLPtime = dwtInput['SLPtime']
SLPdatetime = [datetime(int(r[0]),int(r[1]),int(r[2])) for r in SLPtime]



plt.figure(figsize=(12,7))
p1 = plt.subplot2grid((1,1),(0,0))
p1.plot(timeDataPurdoe,wlPurdoe)
p1.plot(timeDataPurdoe,ssPurdoe)
p1.plot(timeDataPurdoe,seasonalPurdoe)
p1.plot(timeDataPurdoe,mmslaPurdoe)

tCPurdoe = list(timeDataPurdoe)
import pandas as pd
dataPurdoe = np.array([wlPurdoe,ssPurdoe,seasonalPurdoe,mmslaPurdoe,dslaPurdoe,mslPurdoe])
ogdfPurdoe = pd.DataFrame(data=dataPurdoe.T, index=tCPurdoe, columns=["wl", "ss", "seasonal","mmsla","dsla","msl"])
ogdfPurdoe = ogdfPurdoe.loc[SLPdatetime[0]:SLPdatetime[-1]]

dailySSPurdoe = ogdfPurdoe.resample("d")["ss"].max()
dailyMMSLAPurdoe = ogdfPurdoe.resample("d")["mmsla"].mean()
dailySEASONALPurdoe = ogdfPurdoe.resample("d")["seasonal"].mean()
dailyMSLPurdoe = ogdfPurdoe.resample("d")["msl"].mean()
dailyDSLAPurdoe = ogdfPurdoe.resample("d")["dsla"].mean()
badIndPurdoe = np.where(np.isnan(dailyDSLAPurdoe.values))
dailyDSLAPurdoe.values[badIndPurdoe] = 0
badIndssPurdoe = np.where(np.isnan(dailySSPurdoe.values))
dailySSPurdoe.values[badIndssPurdoe] = 0
dailyNTRPurdoe = dailySSPurdoe+dailySEASONALPurdoe+dailyMMSLAPurdoe




plt.figure(figsize=(12,7))
p1 = plt.subplot2grid((1,1),(0,0))
p1.plot(timeDataRed,wlRed)
p1.plot(timeDataRed,ssRed)
p1.plot(timeDataRed,seasonalRed)
p1.plot(timeDataRed,mmslaRed)

tCRed = list(timeDataRed)
import pandas as pd
dataRed = np.array([wlRed,ssRed,seasonalRed,mmslaRed,dslaRed,mslRed])
ogdfRed = pd.DataFrame(data=dataRed.T, index=tCRed, columns=["wl", "ss", "seasonal","mmsla","dsla","msl"])
ogdfRed= ogdfRed.loc[SLPdatetime[0]:SLPdatetime[-1]]
dailySSRed = ogdfRed.resample("d")["ss"].max()
dailyMMSLARed = ogdfRed.resample("d")["mmsla"].mean()
dailySEASONALRed = ogdfRed.resample("d")["seasonal"].mean()
dailyMSLRed = ogdfRed.resample("d")["msl"].mean()
dailyDSLARed = ogdfRed.resample("d")["dsla"].mean()
badIndRed = np.where(np.isnan(dailyDSLARed.values))
dailyDSLARed.values[badIndRed] = 0
badIndssRed = np.where(np.isnan(dailySSRed.values))
dailySSRed.values[badIndssRed] = 0
dailyNTRRed = dailySSRed+dailySEASONALRed+dailyMMSLARed





plt.figure(figsize=(12,7))
p1 = plt.subplot2grid((1,1),(0,0))
p1.plot(timeDataNome,wlNome)
p1.plot(timeDataNome,ssNome)
p1.plot(timeDataNome,seasonalNome)
p1.plot(timeDataNome,mmslaNome)

tCNome = list(timeDataNome)
import pandas as pd
dataNome= np.array([wlNome,ssNome,seasonalNome,mmslaNome,dslaNome,mslNome])
ogdfNome = pd.DataFrame(data=dataNome.T, index=tCNome, columns=["wl", "ss", "seasonal","mmsla","dsla","msl"])
ogdfNome= ogdfNome.loc[SLPdatetime[0]:SLPdatetime[-1]]
dailySSNome = ogdfNome.resample("d")["ss"].max()
dailyMMSLANome = ogdfNome.resample("d")["mmsla"].mean()
dailySEASONALNome = ogdfNome.resample("d")["seasonal"].mean()
dailyMSLNome = ogdfNome.resample("d")["msl"].mean()
dailyDSLANome = ogdfNome.resample("d")["dsla"].mean()
badIndNome = np.where(np.isnan(dailyDSLANome.values))
dailyDSLANome.values[badIndNome] = 0
badIndssNome = np.where(np.isnan(dailySSNome.values))
dailySSNome.values[badIndssNome] = 0
dailyNTRNome = dailySSNome+dailySEASONALNome+dailyMMSLANome
dailyTimeNome = ogdfNome.resample("d")



ogdfSLPs = pd.DataFrame(data=SLPpcs, index=SLPdatetime)

dailySLPtrimmedNome = ogdfSLPs.loc[timeDataNome[0]:timeDataNome[-1],:]
dailySLPtrimmedPurdoe = ogdfSLPs.loc[timeDataPurdoe[0]:timeDataPurdoe[-1],:]
dailySLPtrimmedRed = ogdfSLPs.loc[timeDataRed[0]:timeDataRed[-1],:]


data_folder = "/users/dylananderson/Documents/data/ERSSTv5/"
import xarray as xr
from dateutil.relativedelta import relativedelta

start = 1979
years = np.arange(1979, 2024)
months = np.arange(1, 13)
ogTime = []
for ii in years:
    for hh in months:
        if hh < 10:
            date = str(ii) + "0" + str(hh)
        else:
            date = str(ii) + str(hh)

        file = "ersst.v5." + date + ".nc"
        # print(file)
        if ii == 1940 and hh < 1:
            print("skipping {}/{}".format(ii, hh))
        else:
            if ii == 1979 and hh == 1:
                with xr.open_dataset(os.path.join(data_folder, file)) as ds:
                    temp = ds
                    SSTvalues = ds['sst']
                    ogTime.append(datetime(ii, hh, 1))
            elif ii == (2024 - 1) and hh > 10:
                print("skipping {}/{}".format(ii, hh))
            else:
                with xr.open_dataset(os.path.join(data_folder, file)) as ds:
                    SSTvalues = xr.concat([SSTvalues, ds['sst']], dim="time")
                    ogTime.append(datetime(ii, hh, 1))

dt = datetime(1979, 1, 1)
end = datetime((2023 - 1), 10, 1)
step = relativedelta(months=1)
sstTime = []
while dt < end:
    sstTime.append(dt)
    dt += step

data = SSTvalues.squeeze("lev")

# parse data to xr.Dataset
xds_predictor = xr.Dataset(
    {
        'SST': (('longitude', 'latitude', 'time'), data.data.T),
    },
    coords={
        'longitude': SSTvalues.lon.values,
        'latitude': SSTvalues.lat.values,
        'time': ogTime,
    }
)

var_name = "SST"
y1 = 1979
y2 = 2023
m1 = 1
m2 = 10
subset = xds_predictor.sel(longitude=slice(180, 225), latitude=slice(62, 75))

Xs = subset.longitude.values
Ys = subset.latitude.values
[XR, YR] = np.meshgrid(Xs, Ys)

plt.figure()
p1 = plt.subplot2grid((1, 1), (0, 0))
# spatialField = np.fliplr(subset["SST"][:,:,10])#np.reshape(var_anom_mean.values,(33,36))
spatialField = subset["SST"][:,:,10]#np.reshape(var_anom_mean.values,(33,36))
m = Basemap(projection='merc', llcrnrlat=62, urcrnrlat=75, llcrnrlon=175, urcrnrlon=225, lat_ts=65, resolution='l')
m.drawcoastlines()
cx, cy = m(XR, YR)
CS = m.contourf(cx, cy, spatialField.T),# np.arange(0,0.023,.003), cmap=cm.RdBu_r, shading='gouraud')





from sklearn.decomposition import PCA

nlon, nlat, ntime = np.shape(subset["SST"].values)

collapsed = np.reshape(subset["SST"].values, (nlon * nlat, ntime))




index = ~np.isnan(collapsed[:, 0])
badIndex = np.isnan(collapsed[:, 0])
ocean = [i for i, x in enumerate(index) if x]
land = [i for i, x in enumerate(badIndex) if x]
realDataAnoms = collapsed[index, :]

var_anom_mean = np.nanmean(realDataAnoms.T, axis=0)
var_anom_std = np.nanstd(realDataAnoms.T, axis=0)
timeSeries_mean = np.nanmean(realDataAnoms, axis=0)

nk_m = np.kron(np.ones(((ntime), 1)), var_anom_mean)
nk_s = np.kron(np.ones(((ntime), 1)), var_anom_std)
var_anom_demean = (realDataAnoms.T - nk_m) / nk_s
ipca = PCA()
PCs = ipca.fit_transform(var_anom_demean)

EOFs = ipca.components_
variance = ipca.explained_variance_
nPercent = variance / np.sum(variance)
APEV = np.cumsum(variance) / np.sum(variance) * 100.0
nterm = np.where(APEV <= 0.999 * 100)[0][-1]






data = PCs[:,0:6]
sstdf = pd.DataFrame(data=data, index=ogTime)
dailySST = sstdf.resample('d').interpolate()

dailySSTtrimmedNome = dailySST.loc[ogdfNome.index[0]:ogdfNome.index[-1],:]
dailySSTtrimmedPurdoe = dailySST.loc[ogdfPurdoe.index[0]:ogdfPurdoe.index[-1],:]
dailySSTtrimmedRed = dailySST.loc[ogdfRed.index[0]:ogdfRed.index[-1],:]




dailyVarsRed = np.hstack((dailySSTtrimmedRed,dailySLPtrimmedRed))
dailyVarsPurdoe = np.hstack((dailySSTtrimmedPurdoe,dailySLPtrimmedPurdoe))
dailyVarsNome = np.hstack((dailySSTtrimmedNome,dailySLPtrimmedNome))


nanPurdoe = np.where(np.isnan(dailyNTRPurdoe.values))
nanRed = np.where(np.isnan(dailyNTRRed.values))
nanNome = np.where(np.isnan(dailyNTRNome.values))

dailyVarsRed = np.delete(dailyVarsRed,nanRed[0],axis=0)
dailyVarsPurdoe = np.delete(dailyVarsPurdoe,nanPurdoe[0],axis=0)
dailyVarsNome = np.delete(dailyVarsNome,nanNome[0],axis=0)

ntrPurdoe = np.delete(dailyNTRPurdoe.values,nanPurdoe[0])
ntrNome = np.delete(dailyNTRNome.values,nanNome[0])
ntrRed = np.delete(dailyNTRRed.values,nanRed[0])

predSLP = ogdfSLPs.loc[dailySST.index[0]:dailySST.index[-1],:]
predSST = dailySST.loc[ogdfSLPs.index[0]:ogdfSLPs.index[-1],:]
predVars = np.hstack((predSST,predSLP))
predTime = predSST.index


from sklearn.linear_model import LinearRegression



allPCsToTry = np.arange(0,200)

bestPCPurdoe = []
improvingScorePurdoe = []
for qq in range(50):

    allScoresPurdoe = []
    for yy in range(len(allPCsToTry)):
        if qq == 0:
            tryPCs = yy
            xPurdoe = dailyVarsPurdoe[:, yy].reshape((-1,1))
            yPurdoe = np.array(ntrPurdoe)


        else:
            tryPCs = np.hstack([np.asarray(bestPCPurdoe).flatten(),np.asarray(allPCsToTry[yy]).flatten()])

            xPurdoe = dailyVarsPurdoe[:, tryPCs]#.reshape((-1,1))
            yPurdoe = np.array(ntrPurdoe)


        modelPurdoe = LinearRegression().fit(xPurdoe, yPurdoe)

        r_sqPurdoe = modelPurdoe.score(xPurdoe, yPurdoe)

        allScoresPurdoe.append(r_sqPurdoe)


    bestAdditionIndex = np.argmax(np.array(allScoresPurdoe))
    bestAddition = allPCsToTry[bestAdditionIndex]
    bestPCPurdoe.append(bestAddition)
    removeIndex = np.where(allPCsToTry == bestAddition)
    allPCsToTry = np.delete(allPCsToTry,removeIndex)
    improvingScorePurdoe.append(np.max(np.array(allScoresPurdoe)))
    print('Iter {}: Adding PC#{}, cumulative score:{}'.format(qq,bestAddition,np.max(np.array(allScoresPurdoe))))

allPCsToTry = np.arange(0, 200)

bestPCNome = []
improvingScoreNome = []
for qq in range(50):

    allScoresNome = []

    for yy in range(len(allPCsToTry)):
        if qq == 0:
            tryPCs = yy
            xNome = dailyVarsNome[:, yy].reshape((-1, 1))
            yNome = np.array(ntrNome)

        else:
            tryPCs = np.hstack([np.asarray(bestPCNome).flatten(), np.asarray(allPCsToTry[yy]).flatten()])
            xNome = dailyVarsNome[:, tryPCs]
            yNome = np.array(ntrNome)

        modelNome = LinearRegression().fit(xNome, yNome)
        r_sqNome = modelNome.score(xNome, yNome)
        allScoresNome.append(r_sqNome)

    bestAdditionIndex = np.argmax(np.array(allScoresNome))
    bestAddition = allPCsToTry[bestAdditionIndex]
    bestPCNome.append(bestAddition)
    removeIndex = np.where(allPCsToTry == bestAddition)
    allPCsToTry = np.delete(allPCsToTry, removeIndex)
    improvingScoreNome.append(np.max(np.array(allScoresNome)))
    print('Iter {}: Adding PC#{}, cumulative score:{}'.format(qq, bestAddition, np.max(np.array(allScoresNome))))

allPCsToTry = np.arange(0, 200)

bestPCRed = []
improvingScoreRed = []
for qq in range(50):

    allScoresRed = []
    for yy in range(len(allPCsToTry)):
        if qq == 0:
            tryPCs = yy
            xRed = dailyVarsRed[:, yy].reshape((-1, 1))
            yRed = np.array(ntrRed[0:-1])


        else:
            tryPCs = np.hstack([np.asarray(bestPCRed).flatten(), np.asarray(allPCsToTry[yy]).flatten()])
            xRed = dailyVarsRed[:, tryPCs]  # .reshape((-1,1))
            yRed = np.array(ntrRed[0:-1])

        modelRed = LinearRegression().fit(xRed, yRed)

        r_sqRed = modelRed.score(xRed, yRed)
        allScoresRed.append(r_sqRed)

    bestAdditionIndex = np.argmax(np.array(allScoresRed))
    bestAddition = allPCsToTry[bestAdditionIndex]
    bestPCRed.append(bestAddition)
    removeIndex = np.where(allPCsToTry == bestAddition)
    allPCsToTry = np.delete(allPCsToTry, removeIndex)
    improvingScoreRed.append(np.max(np.array(allScoresRed)))
    print('Iter {}: Adding PC#{}, cumulative score:{}'.format(qq, bestAddition, np.max(np.array(allScoresRed))))



plt.figure()
ax1 = plt.subplot2grid((2,3),(0,1))
ax1.plot(yRed,modelRed.predict(xRed),'.')
ax1.plot([-1.5,2],[-1.5,2],'k--')
ax1.plot([-1.5,2],[0,0],'k--')
ax1.plot([0,0],[-1.5,2],'k--')
ax1.set_xlabel('Red Tide Gauge (m, MSL)')
ax1.set_ylabel('Predicted Non-tidal Residual (m, MSL)')

ax1b = plt.subplot2grid((2,3),(1,1))
yRedStd = np.std(yRed)
yRedMean = np.mean(yRed)
yModRed = modelRed.predict(xRed)
yModRedStd = np.std(yModRed)
yModRedMean = np.mean(yModRed)

ax1b.plot(yRed,((yModRed-yModRedMean)/yModRedStd)*yRedStd+yRedMean,'.')
ax1b.plot([-1.5,2],[-1.5,2],'k--')
ax1b.plot([-1.5,2],[0,0],'k--')
ax1b.plot([0,0],[-1.5,2],'k--')
ax1b.set_xlabel('Red Tide Gauge (m, MSL)')
ax1b.set_ylabel('Predicted Non-tidal Residual (m, MSL)')




ax2 = plt.subplot2grid((2,3),(0,2))
ax2.plot(yPurdoe,modelPurdoe.predict(xPurdoe),'.')
plt.plot([-1,1.5],[-1,1.5],'k--')
ax2.plot([-1,1.5],[0,0],'k--')
ax2.plot([0,0],[-1,1.5],'k--')
ax2.set_xlabel('Purdoe Tide Gauge (m, MSL)')
ax2.set_ylabel('Predicted Non-tidal Residual (m, MSL)')


ax2b = plt.subplot2grid((2,3),(1,2))
yPurdoeStd = np.std(yPurdoe)
yPurdoeMean = np.mean(yPurdoe)
yModPurdoe = modelPurdoe.predict(xPurdoe)
yModPurdoeStd = np.std(yModPurdoe)
yModPurdoeMean = np.mean(yModPurdoe)

ax2b.plot(yPurdoe,((yModPurdoe-yModPurdoeMean)/yModPurdoeStd)*yPurdoeStd+yPurdoeMean,'.')
ax2b.plot([-1,1.5],[-1,1.5],'k--')
ax2b.plot([-1,1.5],[0,0],'k--')
ax2b.plot([0,0],[-1,1.5],'k--')
ax2b.set_xlabel('Purdoe Tide Gauge (m, MSL)')
ax2b.set_ylabel('Predicted Non-tidal Residual (m, MSL)')



ax3 = plt.subplot2grid((2,3),(0,0))
ax3.plot(yNome,modelNome.predict(xNome),'.')
plt.plot([-2.5,2.5],[-2.5,2.5],'k--')
ax3.plot([-2.5,2.5],[0,0],'k--')
ax3.plot([0,0],[-2.5,2.5],'k--')
ax3.set_xlabel('Nome Tide Gauge (m, MSL)')
ax3.set_ylabel('Predicted Non-tidal Residual (m, MSL)')


ax3b = plt.subplot2grid((2,3),(1,0))
yNomeStd = np.std(yNome)
yNomeMean = np.mean(yNome)
yModNome = modelNome.predict(xNome)
yModNomeStd = np.std(yModNome)
yModNomeMean = np.mean(yModNome)

ax3b.plot(yNome,((yModNome-yModNomeMean)/yModNomeStd)*yNomeStd+yNomeMean,'.')
ax3b.plot([-2.5,2.5],[-2.5,2.5],'k--')
ax3b.plot([-2.5,2.5],[0,0],'k--')
ax3b.plot([0,0],[-2.5,2.5],'k--')
ax3b.set_xlabel('Nome Tide Gauge (m, MSL)')
ax3b.set_ylabel('Predicted Non-tidal Residual (m, MSL)')




plt.figure()
ax10 = plt.subplot2grid((3,1),(0,0))
ax10.plot(predTime,((modelNome.predict(predVars[:, bestPCNome])-yModNomeMean)/yModNomeStd)*yNomeStd+yNomeMean)
ax10.set_xlabel('time')
ax10.set_ylabel('Nome NTR Prediction')

ax11 = plt.subplot2grid((3,1),(1,0))
ax11.plot(predTime,((modelRed.predict(predVars[:, bestPCRed])-yModRedMean)/yModRedStd)*yRedStd+yRedMean)
ax11.set_xlabel('time')
ax11.set_ylabel('Red NTR Prediction')

ax12 = plt.subplot2grid((3,1),(2,0))
ax12.plot(predTime,((modelPurdoe.predict(predVars[:, bestPCPurdoe])-yModPurdoeMean)/yModPurdoeStd)*yPurdoeStd+yPurdoeMean)
ax12.set_xlabel('time')
ax12.set_ylabel('Purdoe NTR Prediction')




plt.figure()
ax20 = plt.subplot2grid((5,1),(0,0))
ax20.plot(tideTimeNome,tideWlNome)
ax20.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax20.set_ylim([-2.5,2.5])
ax20.set_ylabel('Observed Nome TG (m)')

ax21 = plt.subplot2grid((5,1),(1,0))
ax21.plot(tideTimeEmNome,mslEmNome)
ax21.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax21.set_ylim([-0.2,0.2])
ax21.set_ylabel('MSL (m)')

ax22 = plt.subplot2grid((5,1),(2,0))
ax22.plot(tideTimeEmNome,tideEmNome)
ax22.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax22.set_ylim([-0.35,0.35])
ax22.set_ylabel('Tide (m)')

ax23 = plt.subplot2grid((5,1),(3,0))
ax23.plot(predTime,((modelNome.predict(predVars[:, bestPCNome])-yModNomeMean)/yModNomeStd)*yNomeStd+yNomeMean)
ax23.set_xlabel('time')
ax23.set_ylabel('NTR (m)')
ax23.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax23.set_ylim([-2.,2.75])


ax24 = plt.subplot2grid((5,1),(4,0))

predNome = ((modelNome.predict(predVars[:, bestPCNome])-yModNomeMean)/yModNomeStd)*yNomeStd+yNomeMean
nomeInd = np.where((predTime>datetime(1979,1,1)) & (predTime<datetime(2024,1,1)))
nomeInd2 = np.where((tideTimeEmNome>datetime(1979,1,1)) & (tideTimeEmNome<datetime(2024,1,1)))

simDailyDeltaPredTime = [(tt - predTime[0]).total_seconds() / (3600 * 24) for tt in predTime]
simHourlyDeltaTime = [(tt - tideTimeEmNome[0]).total_seconds() / (3600 * 24) for tt in tideTimeEmNome[nomeInd2]]

interpWLNome = np.interp(simHourlyDeltaTime, simDailyDeltaPredTime, predNome)

nomeExtendedWL = interpWLNome + tideEmNome[nomeInd2] + mslEmNome[nomeInd2]
ax24.plot(tideTimeEmNome[nomeInd2],nomeExtendedWL)
ax24.set_xlabel('time')
ax24.set_ylabel('Extended WL (m)')
ax24.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax24.set_ylim([-2.,3])








plt.figure()
ax20 = plt.subplot2grid((5,1),(0,0))
ax20.plot(tideTimePurdoe,tideWlPurdoe)
ax20.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax20.set_ylim([-1.5,1.5])
ax20.set_ylabel('Observed Purdoe TG (m)')

ax21 = plt.subplot2grid((5,1),(1,0))
ax21.plot(tideTimeEmPurdoe,mslEmPurdoe)
ax21.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax21.set_ylim([-0.2,0.2])
ax21.set_ylabel('MSL (m)')

ax22 = plt.subplot2grid((5,1),(2,0))
ax22.plot(tideTimeEmPurdoe,tideEmPurdoe)
ax22.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax22.set_ylim([-0.25,0.25])
ax22.set_ylabel('Tide (m)')

ax23 = plt.subplot2grid((5,1),(3,0))
ax23.plot(predTime,modelPurdoe.predict(predVars[:, bestPCPurdoe]))

ax23.set_xlabel('time')
ax23.set_ylabel('NTR (m)')
ax23.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax23.set_ylim([-1.,1.75])


ax24 = plt.subplot2grid((5,1),(4,0))

predPurdoe = ((modelPurdoe.predict(predVars[:, bestPCPurdoe])-yModPurdoeMean)/yModPurdoeStd)*yPurdoeStd+yPurdoeMean
PurdoeInd = np.where((predTime>datetime(1979,1,1)) & (predTime<datetime(2024,1,1)))
PurdoeInd2 = np.where((tideTimeEmPurdoe>datetime(1979,1,1)) & (tideTimeEmPurdoe<datetime(2024,1,1)))
simDailyDeltaPredTime = [(tt - predTime[0]).total_seconds() / (3600 * 24) for tt in predTime]
simHourlyDeltaTime = [(tt - tideTimeEmPurdoe[0]).total_seconds() / (3600 * 24) for tt in tideTimeEmPurdoe[PurdoeInd2]]
interpWLPurdoe = np.interp(simHourlyDeltaTime, simDailyDeltaPredTime, predPurdoe)
PurdoeExtendedWL = interpWLPurdoe + tideEmPurdoe[PurdoeInd2] + mslEmPurdoe[PurdoeInd2]
ax24.plot(tideTimeEmPurdoe[PurdoeInd2],PurdoeExtendedWL)
ax24.set_xlabel('time')
ax24.set_ylabel('Extended WL (m)')
ax24.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax24.set_ylim([-1.,2])







plt.figure()
ax20 = plt.subplot2grid((5,1),(0,0))
ax20.plot(tideTimeRed,tideWlRed)
ax20.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax20.set_ylim([-1.75,1.75])
ax20.set_ylabel('Observed Red TG (m)')

ax21 = plt.subplot2grid((5,1),(1,0))
ax21.plot(tideTimeEmRed,mslEmRed)
ax21.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax21.set_ylim([-0.3,0.2])
ax21.set_ylabel('MSL (m)')

ax22 = plt.subplot2grid((5,1),(2,0))
ax22.plot(tideTimeEmRed,tideEmRed)
ax22.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax22.set_ylim([-0.2,0.2])
ax22.set_ylabel('Tide (m)')

ax23 = plt.subplot2grid((5,1),(3,0))
ax23.plot(predTime,modelRed.predict(predVars[:, bestPCRed]))

ax23.set_xlabel('time')
ax23.set_ylabel('NTR (m)')
ax23.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax23.set_ylim([-0.5,1])


ax24 = plt.subplot2grid((5,1),(4,0))

predRed = ((modelRed.predict(predVars[:, bestPCRed])-yModRedMean)/yModRedStd)*yRedStd+yRedMean
RedInd = np.where((predTime>datetime(1979,1,1)) & (predTime<datetime(2024,1,1)))
RedInd2 = np.where((tideTimeEmRed>datetime(1979,1,1)) & (tideTimeEmRed<datetime(2024,1,1)))
simDailyDeltaPredTime = [(tt - predTime[0]).total_seconds() / (3600 * 24) for tt in predTime]
simHourlyDeltaTime = [(tt - tideTimeEmRed[0]).total_seconds() / (3600 * 24) for tt in tideTimeEmRed[RedInd2]]
interpWLRed = np.interp(simHourlyDeltaTime, simDailyDeltaPredTime, predRed)
RedExtendedWL = interpWLRed + tideEmRed[RedInd2] + mslEmRed[RedInd2]
ax24.plot(tideTimeEmRed[RedInd2],RedExtendedWL)
ax24.set_xlabel('time')
ax24.set_ylabel('Extended WL (m)')
ax24.set_xlim([datetime(1979,1,1),datetime(2023,5,1)])
ax24.set_ylim([-1.5,1.5])





## Ok so we need to fill the gaps....
histTimeRed = tideTimeEmRed[RedInd2]
hoursWithNoWLRed = [x for x in histTimeRed if x not in timeDataRed]
ind_dictRed = dict((k,i) for i,k in enumerate(histTimeRed))
interRed = set(hoursWithNoWLRed).intersection(histTimeRed)
indicesRed = [ ind_dictRed[x] for x in interRed ]
indicesRed.sort()

hoursWithWLRed = [x for x in histTimeRed if x in timeDataRed]
ind_dictRed = dict((k,i) for i,k in enumerate(histTimeRed))
inter2Red = set(hoursWithWLRed).intersection(histTimeRed)
indices2Red = [ ind_dictRed[x] for x in inter2Red]
indices2Red.sort()

gapFilledRedNTR = np.nan * np.ones((len(histTimeRed),))
gapFilledRedNTR[indicesRed] = interpWLRed[indicesRed]
gapFilledRedNTR[indices2Red] = ssRed
badSS = np.where(np.isnan(gapFilledRedNTR))
gapFilledRedNTR[badSS] = interpWLRed[badSS]


histTimeNome = tideTimeEmNome[nomeInd2]
hoursWithNoWLNome = [x for x in histTimeNome if x not in timeDataNome]
ind_dictNome = dict((k,i) for i,k in enumerate(histTimeNome))
interNome = set(hoursWithNoWLNome).intersection(histTimeNome)
indicesNome = [ ind_dictNome[x] for x in interNome ]
indicesNome.sort()

hoursWithWLNome = [x for x in histTimeNome if x in timeDataNome]
ind_dictNome = dict((k,i) for i,k in enumerate(histTimeNome))
inter2Nome = set(hoursWithWLNome).intersection(histTimeNome)
indices2Nome = [ ind_dictNome[x] for x in inter2Nome]
indices2Nome.sort()

gapFilledNomeNTR = np.nan * np.ones((len(histTimeNome),))
gapFilledNomeNTR[indicesNome] = interpWLNome[indicesNome]
gapFilledNomeNTR[indices2Nome] = ssNome
badSS = np.where(np.isnan(gapFilledNomeNTR))
gapFilledNomeNTR[badSS] = interpWLNome[badSS]


histTimePurdoe = tideTimeEmPurdoe[PurdoeInd2]
hoursWithNoWLPurdoe = [x for x in histTimeNome if x not in timeDataPurdoe]
ind_dictPurdoe = dict((k,i) for i,k in enumerate(histTimePurdoe))
interPurdoe = set(hoursWithNoWLPurdoe).intersection(histTimePurdoe)
indicesPurdoe = [ ind_dictPurdoe[x] for x in interPurdoe]
indicesPurdoe.sort()

hoursWithWLPurdoe= [x for x in histTimePurdoe if x in timeDataPurdoe]
ind_dictPurdoe = dict((k,i) for i,k in enumerate(histTimePurdoe))
inter2Purdoe = set(hoursWithWLPurdoe).intersection(histTimePurdoe)
indices2Purdoe = [ ind_dictPurdoe[x] for x in inter2Purdoe]
indices2Purdoe.sort()


gapFilledPurdoeNTR = np.nan * np.ones((len(histTimePurdoe),))
gapFilledPurdoeNTR[indicesPurdoe] = interpWLPurdoe[indicesPurdoe]
gapFilledPurdoeNTR[indices2Purdoe] = ssPurdoe
badSS = np.where(np.isnan(gapFilledPurdoeNTR))
gapFilledPurdoeNTR[badSS] = interpWLPurdoe[badSS]


# %distance weight (roughly) all tide locations
shish_wl = gapFilledNomeNTR *(170/370) + gapFilledRedNTR*(200/370)
pthope_wl = gapFilledPurdoeNTR *(160/920) + gapFilledRedNTR*(760/920)
wainright_wl = gapFilledPurdoeNTR *(175/620) + gapFilledRedNTR*(445/620)

wevok_wl = gapFilledPurdoeNTR*(166/856) + gapFilledRedNTR*(690/856)
wales_wl = gapFilledNomeNTR *(269/450) + gapFilledRedNTR*(181/450)
kivalina_wl = gapFilledRedNTR
ptlay_wl = gapFilledPurdoeNTR *(556/716) + gapFilledRedNTR*(160/716)

historicalPickle = 'historicalNTR.pickle'
outputHistorical = {}
outputHistorical['time'] = histTimePurdoe
outputHistorical['ntr1'] = pthope_wl

with open(historicalPickle,'wb') as f:
    pickle.dump(outputHistorical, f)


import numpy as np
from datetime import datetime, date, timedelta
from scipy.spatial import distance
import pickle
import calendar
import pandas


def toTimestamp(d):
    return calendar.timegm(d.timetuple())

with open(r"dwt49FutureSimulations500.pickle", "rb") as input_file:
    simsInput = pickle.load(input_file)

evbmus_sim = simsInput['evbmus_sim']
dates_sim = simsInput['dates_sim']

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def min_rolling(a, window, axis=1):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    rolling = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    return np.nanmin(rolling, axis=axis)

def max_rolling(a, window, axis=1):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    rolling = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    return np.nanmax(rolling, axis=axis)


directoryDict = ['/volumes/macDrive/arcticSims/pointHope/']

gevCopulaDict = ['gevCopulaSims100000pointHope.pickle']

normHyDict = ['normalizedWaveHydrographsPointHope.pickle']

hyCopDict = ['hydrographCopulaDataPointHope.pickle']

iceSimDict = ['ice.pickle']



for site in range(len(gevCopulaDict)):


    with open(gevCopulaDict[site], "rb") as input_file:
        gevCopulaSimsInput = pickle.load(input_file)
    gevCopulaSims = gevCopulaSimsInput['gevCopulaSims']

    with open(normHyDict[site], "rb") as input_file:
        normalizedWaveHydrographs = pickle.load(input_file)
    normalizedHydros = normalizedWaveHydrographs['normalizedHydros']
    bmuDataMin = normalizedWaveHydrographs['bmuDataMin']
    bmuDataMax = normalizedWaveHydrographs['bmuDataMax']
    bmuDataStd = normalizedWaveHydrographs['bmuDataStd']
    bmuDataNormalized = normalizedWaveHydrographs['bmuDataNormalized']


    with open(hyCopDict[site], "rb") as input_file:
       hydrographCopulaData = pickle.load(input_file)
    copulaData = hydrographCopulaData['copulaData']
    copulaDataNoNaNs = hydrographCopulaData['copulaDataNoNaNs']

    with open(iceSimDict[site], "rb") as input_file:
       iceSimulations = pickle.load(input_file)
    # iceSims = iceSimulations['evbmus_simHist'][212:,:]
    # iceDates = iceSimulations['dates_simHist'][212:]
    iceSims = iceSimulations['evbmus_sim'][16070:,:]
    iceDates = iceSimulations['dates_sim'][16070:]
    # iceConcHist = iceSimulations['iceConcHist']
    iceConcSims = iceSimulations['iceConcSims']
    # iceOnOff = iceSimulations['iceOnOff']
    iceOnOffSims = iceSimulations['iceOnOffSims']





    st = datetime(2023, 6, 1, 0, 0, 0)
    end = datetime(2075, 5, 31, 23, 0, 0)
    step = timedelta(days=1)
    dailyTime = []
    while st < end:
        dailyTime.append(st)#.strftime('%Y-%m-%d'))
        st += step


    st = datetime(2023, 6, 1, 0, 0, 0)
    end = datetime(2075, 6, 1, 0, 0, 0)
    step = timedelta(hours=1)
    hourlyTime = []
    while st < end:
        hourlyTime.append(st)#.strftime('%Y-%m-%d'))
        st += step

    hourlyMonths = month = np.array([tt.month for tt in hourlyTime])
    hourlyHours = np.array([tt.hour for tt in hourlyTime])

    deltaT = [(tt - hourlyTime[0]).total_seconds() / (3600*24) for tt in hourlyTime]


    def closest_node(node, nodes):
        closest_index = distance.cdist([node], nodes).argmin()
        return nodes[closest_index], closest_index


    simulationsHs = list()
    simulationsTp = list()
    simulationsDm = list()
    simulationsNTR = list()
    simulationsTime = list()
    simulationsU10 = list()
    simulationsV10 = list()
    simulationsSSR = list()
    simulationsT2M = list()
    import random


    diffSims = 100
    sim_num = 10

    simCounter = 0
    for simNum in range(diffSims):

        print('saving to ' + directoryDict[site])
        print('filling in simulation #{}'.format(simCounter))
        print('which is a variation of large scale climate #{}'.format(simNum))

        for nextSim in range(sim_num):

            simHs = []
            simTp = []
            simDm = []
            simNTR = []
            simV10 = []
            simU10 = []
            simT2M = []
            simSSR = []
            simTime = []
            for i in range(len(evbmus_sim[:,int(np.floor(np.divide(simCounter,2)))])):
                if np.remainder(i,2000) == 0:
                    print('done with {} hydrographs'.format(i))
                tempBmu = int(evbmus_sim[i,int(np.floor(np.divide(simCounter,2)))]-1)
                # randStorm = random.randint(0, 9999)
                amount = 100
                randStorm = [random.choice(range(70000)) for _ in range(amount)]
                iceAreaBelowTemp = iceOnOffSims[simNum][i+16070]

                stormDetails = gevCopulaSims[tempBmu][randStorm]

                tempSimAreaVar = stormDetails[:,11]-iceAreaBelowTemp
                tempSimAreaFinder = np.where(np.abs(tempSimAreaVar) == np.min(np.abs(tempSimAreaVar)))
                stormDetails = stormDetails[tempSimAreaFinder[0][0],:]

                if stormDetails[0] > 6:
                    print('oh boy, we''ve picked a {}m storm wave in BMU #{}'.format(stormDetails[0],tempBmu))
                durSim = 1#simBmuLengthChopped[simNum][i]

                # simDmNorm = (stormDetails[4] - np.asarray(bmuDataMin)[tempBmu,0]) / (np.asarray(bmuDataMax)[tempBmu,0]-np.asarray(bmuDataMin)[tempBmu,0])
                # simSsNorm = (stormDetails[5] - np.asarray(bmuDataMin)[tempBmu,1]) / (np.asarray(bmuDataMax)[tempBmu,1]-np.asarray(bmuDataMin)[tempBmu,1])
                # test, closeIndex = closest_node([simDmNorm,simSsNorm],np.asarray(bmuDataNormalized)[tempBmu])

                # simHsNorm = (stormDetails[0] - np.asarray(bmuDataMin)[tempBmu,0]) / (np.asarray(bmuDataMax)[tempBmu,0]-np.asarray(bmuDataMin)[tempBmu,0])
                # simTpNorm = (stormDetails[1] - np.asarray(bmuDataMin)[tempBmu,1]) / (np.asarray(bmuDataMax)[tempBmu,1]-np.asarray(bmuDataMin)[tempBmu,1])
                simHsNorm = (stormDetails[4] - np.asarray(bmuDataMin)[tempBmu, 0]) / (
                            np.asarray(bmuDataMax)[tempBmu, 0] - np.asarray(bmuDataMin)[tempBmu, 0])
                simTpNorm = (stormDetails[12] - np.asarray(bmuDataMin)[tempBmu, 1]) / (
                            np.asarray(bmuDataMax)[tempBmu, 1] - np.asarray(bmuDataMin)[tempBmu, 1])

                test, closeIndex = closest_node([simHsNorm,simTpNorm],np.asarray(bmuDataNormalized)[tempBmu])

                actualIndex = closeIndex#int(np.asarray(copulaDataNoNaNs[tempBmu])[closeIndex,6])

                tempHs = ((normalizedHydros[tempBmu][actualIndex]['hsNorm']) * (stormDetails[0]-stormDetails[1]) + stormDetails[1])#.filled()
                tempTp = ((normalizedHydros[tempBmu][actualIndex]['tpNorm']) * (stormDetails[2]-stormDetails[3]) + stormDetails[3])#.filled()
                tempDm = ((normalizedHydros[tempBmu][actualIndex]['dmNorm']) + stormDetails[4])

                tempU10 = ((normalizedHydros[tempBmu][actualIndex]['uNorm']) * (stormDetails[5]-stormDetails[6]) + stormDetails[6])#.filled()
                tempV10 = ((normalizedHydros[tempBmu][actualIndex]['vNorm']) * (stormDetails[7]-stormDetails[8]) + stormDetails[8])#.filled()
                tempSSR = stormDetails[9]
                tempT2M = stormDetails[10]
                # tempT2M = ((normalizedHydros[tempBmu][actualIndex]['t2mNorm']) * (stormDetails[14]-stormDetails[15]) + stormDetails[15])#.filled()

                tempNTR = (normalizedHydros[tempBmu][actualIndex]['ntrNorm']) + stormDetails[12]

                # tempSs = ((normalizedHydros[tempBmu][actualIndex]['ssNorm']) + stormDetails[5])
                if len(normalizedHydros[tempBmu][actualIndex]['hsNorm']) < len(normalizedHydros[tempBmu][actualIndex]['timeNorm']):
                    print('Time is shorter than Hs in bmu {}, index {}'.format(tempBmu,actualIndex))
                if stormDetails[1] < 0:
                    print('woah, we''re less than 0 over here')
                    asdfg
                # if len(tempSs) < len(normalizedHydros[tempBmu][actualIndex]['timeNorm']):
                #     print('Ss is shorter than Time in bmu {}, index {}'.format(tempBmu,actualIndex))
                #     tempLength = len(normalizedHydros[tempBmu][actualIndex]['timeNorm'])
                #     tempSs = np.zeros((len(normalizedHydros[tempBmu][actualIndex]['timeNorm']),))
                #     tempSs[0:len((normalizedHydros[tempBmu][actualIndex]['ssNorm']) + stormDetails[5])] = ((normalizedHydros[tempBmu][actualIndex]['ssNorm']) + stormDetails[5])
                # if len(tempSs) > len(normalizedHydros[tempBmu][actualIndex]['timeNorm']):
                #     print('Now Ss is longer than Time in bmu {}, index {}'.format(tempBmu,actualIndex))
                #     print('{} vs. {}'.format(len(tempSs),len(normalizedHydros[tempBmu][actualIndex]['timeNorm'])))
                #     tempSs = tempSs[0:-1]

                if iceConcSims[i+16070,simNum] < 3:
                    simHs.append(tempHs)
                    simTp.append(tempTp)
                    simDm.append(tempDm)
                    simV10.append(tempV10)
                    simU10.append(tempU10)
                    simSSR.append(tempSSR)
                    simT2M.append(tempT2M)
                    simNTR.append(tempNTR)

                else:
                    simHs.append(tempHs*0)
                    simTp.append(tempTp*0)
                    simDm.append(tempDm*0)
                    simV10.append(tempV10)
                    simU10.append(tempU10)
                    simSSR.append(tempSSR)
                    simT2M.append(tempT2M)
                    simNTR.append(tempNTR)

                    # simSs.append(tempSs)
                #simTime.append(normalizedHydros[tempBmu][actualIndex]['timeNorm']*durSim)
                #dt = np.diff(normalizedHydros[tempBmu][actualIndex]['timeNorm']*durSim)
                simTime.append(np.hstack((np.diff(normalizedHydros[tempBmu][actualIndex]['timeNorm']*durSim), np.diff(normalizedHydros[tempBmu][actualIndex]['timeNorm']*durSim)[-1])))


            cumulativeHours = np.cumsum(np.hstack(simTime))
            newDailyTime = [datetime(2023, 6, 1) + timedelta(days=ii) for ii in cumulativeHours]
            # newDailyTime = [datetime(2022, 6, 1) + timedelta(days=ii) for ii in cumulativeHours]
            simDeltaT = [(tt - newDailyTime[0]).total_seconds() / (3600 * 24) for tt in newDailyTime]
            simDailyDeltaT = [(tt - dailyTime[0]).total_seconds() / (3600 * 24) for tt in dailyTime[0:-1]]

            # simulationsTime.append(newDailyTime)
            # rng = newDailyTime
            #

            # simData = np.array(np.vstack((np.hstack(simHs).T,np.hstack(simTp).T,np.hstack(simDm).T)))
            # simData = np.array(np.vstack((np.hstack(simHs).T,np.hstack(simTp).T,np.hstack(simDm).T,np.hstack(simU10).T,np.hstack(simV10).T,np.hstack(simT2M).T)))
            simData = np.array(np.vstack((np.hstack(simHs).T,np.hstack(simTp).T,np.hstack(simDm).T,np.hstack(simU10).T,np.hstack(simV10).T)))

            # simData = np.array(np.vstack((np.hstack(simHs).T,np.hstack(simTp).T,np.hstack(simDm).T,np.hstack(simSs).T)))
            # simData = np.array((np.ma.asarray(np.hstack(simHs)),np.ma.asarray(np.hstack(simTp)),np.ma.asarray(np.hstack(simDm)),np.ma.asarray(np.hstack(simSs))))
            # simData = np.array([np.hstack(simHs).filled(),np.hstack(simTp).filled(),np.hstack(simDm).filled(),np.hstack(simSs)])

            # ogdf = pandas.DataFrame(data=simData.T,index=newDailyTime,columns=["hs","tp","dm","u10","v10"])

            print('interpolating')
            simDeltaTWaves = np.copy(simDeltaT)
            simHs = np.hstack(simHs)
            simTp = np.hstack(simTp)
            simDm = np.hstack(simDm)
            whereNan = np.where((np.isnan(simHs)))
            simHs = np.delete(simHs, whereNan)
            simTp = np.delete(simTp, whereNan)
            simDm = np.delete(simDm, whereNan)
            simDeltaTWaves = np.delete(simDeltaTWaves, whereNan)


            interpHs = np.interp(deltaT, simDeltaTWaves, simHs)
            interpTp = np.interp(deltaT, simDeltaTWaves, simTp)
            interpDm = np.interp(deltaT, simDeltaTWaves, simDm)
            interpU10 = np.interp(deltaT,simDeltaT,np.hstack(simU10))
            interpV10 = np.interp(deltaT,simDeltaT,np.hstack(simV10))
            interpSSR = np.interp(deltaT,simDailyDeltaT,np.hstack(simSSR))
            # interpT2M = np.interp(deltaT, simDeltaT, np.hstack(simT2M))
            # interpT2M = np.interp(deltaT,simDailyDeltaT,np.hstack(simT2M))
            interpNTR = np.interp(deltaT, simDeltaT, np.hstack(simNTR))


            beginOfDayIndex = np.where(hourlyHours==0)
            beginOfDayIndexPlus1 = np.where(hourlyHours==1)
            beginOfDayIndexPlus2 = np.where(hourlyHours==2)
            beginOfDayIndexPlus3 = np.where(hourlyHours==3)

            endOfDayIndex = np.where(hourlyHours==23)
            endOfDayIndexMinus1 = np.where(hourlyHours==22)
            endOfDayIndexMinus2 = np.where(hourlyHours==21)
            endOfDayIndexMinus3 = np.where(hourlyHours==20)

            endOfDayHsMinus3 = interpHs[endOfDayIndexMinus3]
            endOfDayHsMinus2 = interpHs[endOfDayIndexMinus2]
            endOfDayHsMinus1 = interpHs[endOfDayIndexMinus1]
            endOfDayHs = interpHs[endOfDayIndex]
            beginOfDayHs = interpHs[beginOfDayIndex]
            beginOfDayHsPlus1 = interpHs[beginOfDayIndexPlus1]
            beginOfDayHsPlus2 = interpHs[beginOfDayIndexPlus2]
            beginOfDayHsPlus3 = interpHs[beginOfDayIndexPlus3]

            newEndOfDayHsMinus1 = np.divide(endOfDayHsMinus3[0:-1]+2*endOfDayHsMinus2[0:-1]+2*endOfDayHsMinus1[0:-1]+2*endOfDayHs[0:-1]+beginOfDayHs[1:],8)
            newEndOfDayHs = np.divide(beginOfDayHs[1:]+beginOfDayHsPlus1[1:]+endOfDayHsMinus2[0:-1]+endOfDayHsMinus1[0:-1]+endOfDayHs[0:-1],5)
            newBeginOfDayHs = np.divide(beginOfDayHs[1:]+beginOfDayHsPlus1[1:]+beginOfDayHsPlus2[1:]+endOfDayHs[0:-1]+endOfDayHsMinus1[0:-1],5)
            newBeginOfDayHsPlus1 = np.divide(2*beginOfDayHs[1:]+2*beginOfDayHsPlus1[1:]+2*beginOfDayHsPlus2[1:]+beginOfDayHsPlus3[1:]+endOfDayHs[0:-1],8)



            interpHs[endOfDayIndexMinus1[0][0:-1]] = newEndOfDayHsMinus1
            interpHs[endOfDayIndex[0][0:-1]] = newEndOfDayHs
            interpHs[beginOfDayIndex[0][1:]] = newBeginOfDayHs
            interpHs[beginOfDayIndexPlus1[0][1:]] = newBeginOfDayHsPlus1

            interpHsMax = np.hstack([0,0,max_rolling(interpHs,5),0,0])
            # interpTpMax = np.hstack([0,0,max_rolling(interpTp,5),0,0])
            # interpDmMax = np.hstack([0,0,max_rolling(interpDm,5),0,0])

            interpHsMin = np.hstack([0,0,min_rolling(interpHs,5),0,0])
            # interpTpMin = np.hstack([0,0,min_rolling(interpTp,5),0,0])
            # interpDmMin = np.hstack([0,0,min_rolling(interpDm,5),0,0])

            earlyYear = np.where(hourlyMonths < 7)
            lateYear = np.where(hourlyMonths >= 9)
            Hs = interpHs
            Hs[earlyYear] = interpHsMin[earlyYear]
            Hs[lateYear] = interpHsMax[lateYear]



            endOfDayTpMinus3 = interpTp[endOfDayIndexMinus3]
            endOfDayTpMinus2 = interpTp[endOfDayIndexMinus2]
            endOfDayTpMinus1 = interpTp[endOfDayIndexMinus1]
            endOfDayTp = interpTp[endOfDayIndex]
            beginOfDayTp = interpTp[beginOfDayIndex]
            beginOfDayTpPlus1 = interpTp[beginOfDayIndexPlus1]
            beginOfDayTpPlus2 = interpTp[beginOfDayIndexPlus2]
            beginOfDayTpPlus3 = interpTp[beginOfDayIndexPlus3]

            newEndOfDayTpMinus1 = np.divide(endOfDayTpMinus3[0:-1]+2*endOfDayTpMinus2[0:-1]+2*endOfDayTpMinus1[0:-1]+2*endOfDayTp[0:-1]+beginOfDayTp[1:],8)
            newEndOfDayTp = np.divide(beginOfDayTp[1:]+beginOfDayTpPlus1[1:]+endOfDayTpMinus2[0:-1]+endOfDayTpMinus1[0:-1]+endOfDayTp[0:-1],5)
            newBeginOfDayTp = np.divide(beginOfDayTp[1:]+beginOfDayTpPlus1[1:]+beginOfDayTpPlus2[1:]+endOfDayTp[0:-1]+endOfDayTpMinus1[0:-1],5)
            newBeginOfDayTpPlus1 = np.divide(2*beginOfDayTp[1:]+2*beginOfDayTpPlus1[1:]+2*beginOfDayTpPlus2[1:]+beginOfDayTpPlus3[1:]+endOfDayTp[0:-1],8)

            interpTp[endOfDayIndexMinus1[0][0:-1]] = newEndOfDayTpMinus1
            interpTp[endOfDayIndex[0][0:-1]] = newEndOfDayTp
            interpTp[beginOfDayIndex[0][1:]] = newBeginOfDayTp
            interpTp[beginOfDayIndexPlus1[0][1:]] = newBeginOfDayTpPlus1

            # # plt.plot(hourlyTime,interpHs)
            # plt.figure()
            # plt.plot(hourlyTime,interpNTR)

            endOfDayNTRMinus3 = interpNTR[endOfDayIndexMinus3]
            endOfDayNTRMinus2 = interpNTR[endOfDayIndexMinus2]
            endOfDayNTRMinus1 = interpNTR[endOfDayIndexMinus1]
            endOfDayNTR = interpNTR[endOfDayIndex]
            beginOfDayNTR = interpNTR[beginOfDayIndex]
            beginOfDayNTRPlus1 = interpNTR[beginOfDayIndexPlus1]
            beginOfDayNTRPlus2 = interpNTR[beginOfDayIndexPlus2]
            beginOfDayNTRPlus3 = interpNTR[beginOfDayIndexPlus3]

            newEndOfDayNTRMinus1 = np.divide(endOfDayNTRMinus3[0:-1]+2*endOfDayNTRMinus2[0:-1]+2*endOfDayNTRMinus1[0:-1]+2*endOfDayNTR[0:-1]+beginOfDayNTR[1:],8)
            newEndOfDayNTR = np.divide(beginOfDayNTR[1:]+beginOfDayNTRPlus1[1:]+endOfDayNTRMinus2[0:-1]+endOfDayNTRMinus1[0:-1]+endOfDayNTR[0:-1],5)
            newBeginOfDayNTR = np.divide(beginOfDayNTR[1:]+beginOfDayNTRPlus1[1:]+beginOfDayNTRPlus2[1:]+endOfDayNTR[0:-1]+endOfDayNTRMinus1[0:-1],5)
            newBeginOfDayNTRPlus1 = np.divide(2*beginOfDayNTR[1:]+2*beginOfDayNTRPlus1[1:]+2*beginOfDayNTRPlus2[1:]+beginOfDayNTRPlus3[1:]+endOfDayNTR[0:-1],8)

            interpNTR[endOfDayIndexMinus1[0][0:-1]] = newEndOfDayNTRMinus1
            interpNTR[endOfDayIndex[0][0:-1]] = newEndOfDayNTR
            interpNTR[beginOfDayIndex[0][1:]] = newBeginOfDayNTR
            interpNTR[beginOfDayIndexPlus1[0][1:]] = newBeginOfDayNTRPlus1

            # plt.plot(hourlyTime,interpNTR)

            dailyMonths = month = np.array([tt.month for tt in dailyTime])

            earlyDYear = np.where(dailyMonths[0:-1] < 5)
            lateDYear = np.where(dailyMonths[0:-1] >= 10)
            midDYear = np.where(dailyMonths[0:-1] == 8)

            T2M = np.hstack(simT2M)
            T2Mmaxrolling = np.hstack([np.nanmean(T2M),max_rolling(T2M,3),np.nanmean(T2M)])
            T2Mminrolling = np.hstack([np.nanmean(T2M),np.nanmean(T2M),min_rolling(T2M,5),np.nanmean(T2M),np.nanmean(T2M)])
            T2Mrolling = np.hstack([np.nanmean(T2M),np.nanmean(T2M),moving_average(T2M,5),np.nanmean(T2M),np.nanmean(T2M)])

            # interpMinT2M = np.interp(deltaT,simDailyDeltaT,T2Mminrolling)
            # interpMaxT2M = np.interp(deltaT,simDailyDeltaT,T2Mmaxrolling)

            newT2M = T2Mrolling
            newT2M[earlyDYear[0]] = T2Mminrolling[earlyDYear[0]]
            newT2M[lateDYear[0]] = T2Mminrolling[lateDYear[0]]
            newT2M[midDYear[0]] = T2Mmaxrolling[midDYear[0]]

            interpT2M = np.interp(deltaT,simDailyDeltaT,newT2M)

            # plt.figure()
            # p1 = plt.subplot2grid((2,1),(0,0))
            # # p1.plot(newDailyTime,np.hstack(simT2M))
            # p1.plot(simDailyDeltaT,np.hstack(simT2M))
            # p1.plot(simDailyDeltaT,T2Mmoving)
            #
            # p2 = plt.subplot2grid((2,1),(1,0))
            # p2.plot(hourlyTime,interpT2M)

            # interpSs = np.interp(deltaT,simDeltaT,np.hstack(simSs))

            simDataInterp = np.array([Hs,interpTp,interpDm,interpU10,interpV10,interpSSR,interpT2M,interpNTR])

            # df = pandas.DataFrame(data=simDataInterp.T,index=hourlyTime,columns=["hs","tp","dm"])
            df = pandas.DataFrame(data=simDataInterp.T,index=hourlyTime,columns=["hs","tp","dm","u10","v10","ssr","t2m","ntr"])
            # resampled = df.resample('H')
            # interped = resampled.interpolate()
            # simulationData = interped.values
            # testTime = interped.index  # to_pydatetime()
            # testTime2 = testTime.to_pydatetime()
            # simsPickle = ('/Users/dylananderson/Documents/data/wainwright/futureSimulation{}.pickle'.format(simCounter))
            simsPickle = (directoryDict[site] + 'futureSimulation{}.pickle'.format(simCounter))
            # simsPickle = ('/Users/dylananderson/Documents/data/pointLay/futureSimulation{}.pickle'.format(simNum))
            # simsPickle = ('/Users/dylananderson/Documents/data/pointHope/futureSimulation{}.pickle'.format(simCounter))
            # simsPickle = ('/volumes/macDrive/historicalSims2/simulation{}.pickle'.format(simNum))

            outputSims= {}
            outputSims['simulationData'] = simDataInterp.T
            # outputSims['hourlyTime'] = hourlyTime
            outputSims['df'] = df
            outputSims['time'] = hourlyTime
            # outputSims['simHs'] = np.hstack(simHs)
            # outputSims['simTp'] = np.hstack(simTp)
            # outputSims['simDm'] = np.hstack(simDm)
            # outputSims['simU10'] = np.hstack(simU10)
            # outputSims['simV10'] = np.hstack(simV10)
            # outputSims['simSSR'] = np.hstack(simSSR)
            # outputSims['simT2M'] = np.hstack(simT2M)
            # outputSims['simNTR'] = np.hstack(simNTR)

            # outputSims['simSs'] = np.hstack(simSs)

            with open(simsPickle, 'wb') as f:
                pickle.dump(outputSims, f)
            simCounter = simCounter+1


    #
    # # ts = pandas.Series(np.hstack(simHs), index=newDailyTime)
    # # resampled = ts.resample('H')
    # # interp = resampled.interpolate()
    #
    # testTime = interped.index  # to_pydatetime()
    # testTime2 = testTime.to_pydatetime()
    # testData = interped.values



# simsPickle = '/home/dylananderson/projects/atlanticClimate/Sims/allSimulations.pickle'
# outputSims= {}
# outputSims['simulationsHs'] = simulationsHs
# outputSims['simulationsTime'] = simulationsTime
# outputSims['simulationsTp'] = simulationsTp
# outputSims['simulationsDm'] = simulationsDm
# outputSims['simulationsSs'] = simulationsSs
#
# with open(simsPickle, 'wb') as f:
#     pickle.dump(outputSims, f)
#
#
# plt.figure()
# ax1 = plt.subplot2grid((3,1),(0,0),rowspan=1,colspan=1)
# hs = simDataInterp[0,:]
# where0 = np.where((hs == 0))
# hs[where0] = np.nan
# ax1.plot(hourlyTime,hs)
# ax2 = plt.subplot2grid((3,1),(1,0),rowspan=1,colspan=1)
# tp = simDataInterp[1,:]
# where0 = np.where((tp < 0.5))
# tp[where0] = np.nan
# ax2.plot(hourlyTime,tp)
# ax3 = plt.subplot2grid((3,1),(2,0),rowspan=1,colspan=1)
# dm = simDataInterp[7,:]
# # where0 = np.where((dm == 0))
# # dm[where0] = np.nan
# # where360 = np.where((dm > 360))
# # dm[where360] = dm[where360]-360
# # whereNeg = np.where((dm < 0))
# # dm[whereNeg] = dm[whereNeg]+360
# ax3.plot(hourlyTime,dm)
#
# ### TODO: Need to assess the statistics of these hypothetical scenarios... Yearly max Hs? Wave Energy?
#
# ### TODO: Which requires interpolating the time series to hourly values...
#
# # for qq in len(simulationsTime):
#

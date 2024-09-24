from functions.metOcean import getMetOcean

# Adjust the start and end times to download more or less data for the emulator
startTime = [2023, 8, 1]
endTime = [2023, 9, 30]

# point hope
lonLeft = -166.8
lonRight = -166.6
latBot = 68.4
latTop = 68.6

metOcean = getMetOcean(shoreNormal=0,lonLeft=lonLeft, lonRight=lonRight, latBot=latBot, latTop=latTop, startTime=startTime, endTime=endTime)
metOcean.getERA5WavesAndWindsAndTemps(printToScreen=True)

import pickle
outdict = {}
outdict['metOcean'] = metOcean
outdict['endTime'] = endTime
outdict['startTime'] = startTime
outdict['lonLeft'] = lonLeft
outdict['lonRight'] = lonRight
outdict['latBot'] = latBot
outdict['latTop'] = latTop

with open('waves.pickle', 'wb') as f:
    pickle.dump(outdict, f)



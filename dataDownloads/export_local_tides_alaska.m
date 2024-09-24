addpath(genpath('/Users/dylananderson/Documents/projects/arcticClimate/'))
addpath(genpath('/Users/dylananderson/Documents/data/chsWest/'))

%download tides
nome.tideout = download_noaa_tides('9468756', 'MSL', 1999, 2023);
red.tideout = download_noaa_tides('9491094', 'MSL', 2003, 2023);
pru.tideout = download_noaa_tides('9497645', 'MSL', 1990, 2023);

%interolate tides
times = datenum(1999,1,1):1/24:datenum(2023,1,1);
nome_wl = interp1(nome.tideout.mtime, nome.tideout.wl, times);
red_wl = interp1(red.tideout.mtime, red.tideout.wl, times);
pru_wl = interp1(pru.tideout.mtime, pru.tideout.wl, times);

%download predicted tides
nome.tideoutwithPred = download_noaa_tides_withPred('9468756', 'MSL', 1999, 2024);
red.tideoutwithPred = download_noaa_tides_withPred('9491094', 'MSL', 2003, 2024);
pru.tideoutwithPred = download_noaa_tides_withPred('9497645', 'MSL', 1990, 2024);




%%

dat = pru.tideoutwithPred.wl;
time = pru.tideoutwithPred.wltime;
time_tide = pru.tideoutwithPred.predtime(1:10:end);
tide = pru.tideoutwithPred.pred(1:10:end);
lin = 1;
[MSL,~,MMSLA,MMSLA_hourly,DSLA,climatology,...
    climatologyDaily,month_time,parameters] = splitNTR_v2(dat,time,lin,savepath,fname);


wl = dat;
dat = DSLA;

% Use the tide to fill in spaces in the time series for spectral analysis
noaaTime = time_tide;


% Make sure tidal time matches up with wl time
 if length(time)>length(noaaTime)
     noaaTide = nan(length(time),1);
     tideData = dat;
     [~,ia,ib] = intersect(noaaTime,time);
     noaaTide(ib) = tide(ia);
     time2 = time;
     
 elseif length(time)<length(noaaTime)
     noaaTide = tide;
     tideData = nan(length(time_tide),1);
     [~,ia,ib] = intersect(noaaTime,time);
     tideData(ia) = dat(ib);    
     time2 = time;
     time = time_tide;
     
 else
     noaaTide = tide;
     tideData = dat;
 end
 
 
% find where all the NaNs are in time
removeNan = find(isnan(tideData)>0);
predictInd = find(isnan(noaaTide)<1);
 
[~,ia,ib] = intersect(removeNan,predictInd);

% Add in the noaa predicted tide.
tideDataNoGap = tideData;
tideDataNoGap(removeNan(ia)) = noaaTide(predictInd(ib))-nanmean(noaaTide);

tideData = tideDataNoGap;
moo = find(isnan(tideData)<2);
%No gaps longer exist - remove tides
valu = 3; %previously was 2...
bloc = 365.25*24*valu;  
dt = 1;
overlap = 0.50;
 

%mkdir(fname)
%cd(fname)
outdir = [pwd,'/'];
% Choose the center of each of the peaks in energy, represented in cycles
% per day
bandCenter = [0.98,1.95,2.9,3.9,4.86,5.83,6.80,7.80,8.78,9.7,10.8]; % SF data
bl = 0.25;

bandLow = bandCenter - bl;
bandHigh = bandCenter + bl;
tideBands = [bandLow; bandHigh]';
num =1;
fig_val = 0; % 1 = save figures, 0 = no figures
findNTRspectral(tideData(moo),time(moo),bloc,overlap,dt,outdir,tideBands,num,fig_val);
close all

cd("26298hr/")
load('SpectralData1.mat');

%time = combinedWaveTide(:,1);
% Find where the data was nan and get rid
TP_SS = nan(length(time),1);
[~,ia,ib] = intersect(nonTideTime,time);
TP_SS(ib) = nonTideData(ia);


ind = find(isnan(wl)>0);
TP_SS(ind) = NaN;

% Here, we remove the seasonal cycle from tide predictions (SA and SSA
% constituents). We do this as to not double add seasonal cycles in.


f = 1/365.25;
cnt = 0;
rmse = nan(1,4);
[~,ia,ib] = intersect(time_tide,time2);
ind = ib;
Y = tide(ia);

for ii = 1:4
   cnt = cnt+1;
  
    if cnt == 1
        X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia));
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 2
         X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia));
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 3
         X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) ];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia));
        
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 4
        X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) sin(4*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia))+B(5)*sin(4*pi*f*time_tide(ia));

        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    end
  
end

[~,indy] = min(rmse);


if indy == 1
   X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia))];
   [B] = regress(Y,X);
   New = B(1)+B(2)*cos(2*pi*f*time_tide(ia));

elseif indy == 2
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia))];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia));

elseif indy == 3
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) ];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia));
        

elseif indy == 4
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) sin(4*pi*f*time_tide(ia))];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia))+B(5)*sin(4*pi*f*time_tide(ia));
       
end
f = 1/365.25;
X = [ones(length(tide),1) cos(2*pi*f*time_tide) sin(2*pi*f*time_tide) cos(4*pi*f*time_tide) sin(4*pi*f*time_tide)];
Y = tide;
[B] = regress(Y,X);
New = B(1)+B(2)*cos(2*pi*f*time_tide)+B(3)*sin(2*pi*f*time_tide)+B(4)*cos(4*pi*f*time_tide)+B(5)*sin(4*pi*f*time_tide);
tide = tide-New;

[~,ia,ib] = intersect(time,time_tide);
purdoeDailyData.tide = nan(length(time),1);
purdoeDailyData.wl = nan(length(time),1);
purdoeDailyData.seasonal = nan(length(time),1);
purdoeDailyData.msl = nan(length(time),1);
purdoeDailyData.mmsla = nan(length(time),1);
purdoeDailyData.dsla = nan(length(time),1);
purdoeDailyData.tide(ia) = tide(ib);

[~,ia,ib] = intersect(time,time2);

purdoeDailyData.time = time;
%nomeDailyData.swh = combinedWaveTide(:,3);
%nomeDailyData.tp = combinedWaveTide(:,2);
%nomeDailyData.mwd = combinedWaveTide(:,4);
purdoeDailyData.wl(ia) = wl(ib);
purdoeDailyData.seasonal(ia) = climatologyDaily(ib);
dailynomeDailyDataData.msl(ia) = MSL(ib);
purdoeDailyData.mmsla(ia) = MMSLA_hourly(ib);
purdoeDailyData.dsla(ia) = DSLA(ib);
purdoeDailyData.ss = TP_SS;
purdoeDailyData.seasonal_month = climatology;
purdoeDailyData.mmsla_month = MMSLA;
purdoeDailyData.mmsla_month_time = month_time;
purdoeDailyData.B.msl = parameters.msl;
purdoeDailyData.B.season = parameters.seasons;

purdoeDailyData.monthDateVec = datevec(month_time);
purdoeDailyData.hourlyDateVec = datevec(time);


%%
dat = nome.tideoutwithPred.wl;
time = nome.tideoutwithPred.wltime;
time_tide = nome.tideoutwithPred.predtime(1:10:end);
tide = nome.tideoutwithPred.pred(1:10:end);
lin = 1;
[MSL,~,MMSLA,MMSLA_hourly,DSLA,climatology,...
    climatologyDaily,month_time,parameters] = splitNTR_v2(dat,time,lin,savepath,fname);


wl = dat;
dat = DSLA;

% Use the tide to fill in spaces in the time series for spectral analysis
noaaTime = time_tide;


% Make sure tidal time matches up with wl time
 if length(time)>length(noaaTime)
     noaaTide = nan(length(time),1);
     tideData = dat;
     [~,ia,ib] = intersect(noaaTime,time);
     noaaTide(ib) = tide(ia);
     time2 = time;
     
 elseif length(time)<length(noaaTime)
     noaaTide = tide;
     tideData = nan(length(time_tide),1);
     [~,ia,ib] = intersect(noaaTime,time);
     tideData(ia) = dat(ib);    
     time2 = time;
     time = time_tide;
     
 else
     noaaTide = tide;
     tideData = dat;
 end
 
 
% find where all the NaNs are in time
removeNan = find(isnan(tideData)>0);
predictInd = find(isnan(noaaTide)<1);
 
[~,ia,ib] = intersect(removeNan,predictInd);

% Add in the noaa predicted tide.
tideDataNoGap = tideData;
tideDataNoGap(removeNan(ia)) = noaaTide(predictInd(ib))-nanmean(noaaTide);

tideData = tideDataNoGap;
moo = find(isnan(tideData)<2);
%No gaps longer exist - remove tides
valu = 3; %previously was 2...
bloc = 365.25*24*valu;  
dt = 1;
overlap = 0.50;
 

%mkdir(fname)
%cd(fname)
outdir = [pwd,'/'];
% Choose the center of each of the peaks in energy, represented in cycles
% per day
bandCenter = [0.98,1.95,2.9,3.9,4.86,5.83,6.80,7.80,8.78,9.7,10.8]; % SF data
bl = 0.25;

bandLow = bandCenter - bl;
bandHigh = bandCenter + bl;
tideBands = [bandLow; bandHigh]';
num =1;
fig_val = 0; % 1 = save figures, 0 = no figures
findNTRspectral(tideData(moo),time(moo),bloc,overlap,dt,outdir,tideBands,num,fig_val);
close all

cd("26298hr/")
load('SpectralData1.mat');

%time = combinedWaveTide(:,1);
% Find where the data was nan and get rid
TP_SS = nan(length(time),1);
[~,ia,ib] = intersect(nonTideTime,time);
TP_SS(ib) = nonTideData(ia);


ind = find(isnan(wl)>0);
TP_SS(ind) = NaN;

% Here, we remove the seasonal cycle from tide predictions (SA and SSA
% constituents). We do this as to not double add seasonal cycles in.


f = 1/365.25;
cnt = 0;
rmse = nan(1,4);
[~,ia,ib] = intersect(time_tide,time2);
ind = ib;
Y = tide(ia);

for ii = 1:4
   cnt = cnt+1;
  
    if cnt == 1
        X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia));
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 2
         X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia));
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 3
         X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) ];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia));
        
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 4
        X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) sin(4*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia))+B(5)*sin(4*pi*f*time_tide(ia));

        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    end
  
end

[~,indy] = min(rmse);


if indy == 1
   X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia))];
   [B] = regress(Y,X);
   New = B(1)+B(2)*cos(2*pi*f*time_tide(ia));

elseif indy == 2
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia))];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia));

elseif indy == 3
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) ];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia));
        

elseif indy == 4
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) sin(4*pi*f*time_tide(ia))];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia))+B(5)*sin(4*pi*f*time_tide(ia));
       
end
f = 1/365.25;
X = [ones(length(tide),1) cos(2*pi*f*time_tide) sin(2*pi*f*time_tide) cos(4*pi*f*time_tide) sin(4*pi*f*time_tide)];
Y = tide;
[B] = regress(Y,X);
New = B(1)+B(2)*cos(2*pi*f*time_tide)+B(3)*sin(2*pi*f*time_tide)+B(4)*cos(4*pi*f*time_tide)+B(5)*sin(4*pi*f*time_tide);
tide = tide-New;

[~,ia,ib] = intersect(time,time_tide);
nomeDailyData.tide = nan(length(time),1);
nomeDailyData.wl = nan(length(time),1);
nomeDailyData.seasonal = nan(length(time),1);
nomeDailyData.msl = nan(length(time),1);
nomeDailyData.mmsla = nan(length(time),1);
nomeDailyData.dsla = nan(length(time),1);
nomeDailyData.tide(ia) = tide(ib);

[~,ia,ib] = intersect(time,time2);

nomeDailyData.time = time;
%nomeDailyData.swh = combinedWaveTide(:,3);
%nomeDailyData.tp = combinedWaveTide(:,2);
%nomeDailyData.mwd = combinedWaveTide(:,4);
nomeDailyData.wl(ia) = wl(ib);
nomeDailyData.seasonal(ia) = climatologyDaily(ib);
nomeDailyData.msl(ia) = MSL(ib);
nomeDailyData.mmsla(ia) = MMSLA_hourly(ib);
nomeDailyData.dsla(ia) = DSLA(ib);
nomeDailyData.ss = TP_SS;
nomeDailyData.seasonal_month = climatology;
nomeDailyData.mmsla_month = MMSLA;
nomeDailyData.mmsla_month_time = month_time;
nomeDailyData.B.msl = parameters.msl;
nomeDailyData.B.season = parameters.seasons;

nomeDailyData.monthDateVec = datevec(month_time);
nomeDailyData.hourlyDateVec = datevec(time);



%%
dat = red.tideoutwithPred.wl;
time = red.tideoutwithPred.wltime;
time_tide = red.tideoutwithPred.predtime(1:10:end);
tide = red.tideoutwithPred.pred(1:10:end);
lin = 1;
[MSL,~,MMSLA,MMSLA_hourly,DSLA,climatology,...
    climatologyDaily,month_time,parameters] = splitNTR_v2(dat,time,lin,savepath,fname);


wl = dat;
dat = DSLA;

% Use the tide to fill in spaces in the time series for spectral analysis
noaaTime = time_tide;


% Make sure tidal time matches up with wl time
 if length(time)>length(noaaTime)
     noaaTide = nan(length(time),1);
     tideData = dat;
     [~,ia,ib] = intersect(noaaTime,time);
     noaaTide(ib) = tide(ia);
     time2 = time;
     
 elseif length(time)<length(noaaTime)
     noaaTide = tide;
     tideData = nan(length(time_tide),1);
     [~,ia,ib] = intersect(noaaTime,time);
     tideData(ia) = dat(ib);    
     time2 = time;
     time = time_tide;
     
 else
     noaaTide = tide;
     tideData = dat;
 end
 
 
% find where all the NaNs are in time
removeNan = find(isnan(tideData)>0);
predictInd = find(isnan(noaaTide)<1);
 
[~,ia,ib] = intersect(removeNan,predictInd);

% Add in the noaa predicted tide.
tideDataNoGap = tideData;
tideDataNoGap(removeNan(ia)) = noaaTide(predictInd(ib))-nanmean(noaaTide);

tideData = tideDataNoGap;
moo = find(isnan(tideData)<2);
%No gaps longer exist - remove tides
valu = 3; %previously was 2...
bloc = 365.25*24*valu;  
dt = 1;
overlap = 0.50;
 

%mkdir(fname)
%cd(fname)
outdir = [pwd,'/'];
% Choose the center of each of the peaks in energy, represented in cycles
% per day
bandCenter = [0.98,1.95,2.9,3.9,4.86,5.83,6.80,7.80,8.78,9.7,10.8]; % SF data
bl = 0.25;

bandLow = bandCenter - bl;
bandHigh = bandCenter + bl;
tideBands = [bandLow; bandHigh]';
num =1;
fig_val = 0; % 1 = save figures, 0 = no figures
findNTRspectral(tideData(moo),time(moo),bloc,overlap,dt,outdir,tideBands,num,fig_val);
close all

cd("26298hr/")
load('SpectralData1.mat');

%time = combinedWaveTide(:,1);
% Find where the data was nan and get rid
TP_SS = nan(length(time),1);
[~,ia,ib] = intersect(nonTideTime,time);
TP_SS(ib) = nonTideData(ia);


ind = find(isnan(wl)>0);
TP_SS(ind) = NaN;

% Here, we remove the seasonal cycle from tide predictions (SA and SSA
% constituents). We do this as to not double add seasonal cycles in.


f = 1/365.25;
cnt = 0;
rmse = nan(1,4);
[~,ia,ib] = intersect(time_tide,time2);
ind = ib;
Y = tide(ia);

for ii = 1:4
   cnt = cnt+1;
  
    if cnt == 1
        X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia));
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 2
         X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia));
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 3
         X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) ];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia));
        
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    elseif cnt == 4
        X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) sin(4*pi*f*time_tide(ia))];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia))+B(5)*sin(4*pi*f*time_tide(ia));

        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New).^2)/length(ind));
    end
  
end

[~,indy] = min(rmse);


if indy == 1
   X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia))];
   [B] = regress(Y,X);
   New = B(1)+B(2)*cos(2*pi*f*time_tide(ia));

elseif indy == 2
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia))];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia));

elseif indy == 3
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) ];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia));
        

elseif indy == 4
    X = [ones(length(tide(ia)),1) cos(2*pi*f*time_tide(ia)) sin(2*pi*f*time_tide(ia)) cos(4*pi*f*time_tide(ia)) sin(4*pi*f*time_tide(ia))];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time_tide(ia))+B(3)*sin(2*pi*f*time_tide(ia))+B(4)*cos(4*pi*f*time_tide(ia))+B(5)*sin(4*pi*f*time_tide(ia));
       
end
f = 1/365.25;
X = [ones(length(tide),1) cos(2*pi*f*time_tide) sin(2*pi*f*time_tide) cos(4*pi*f*time_tide) sin(4*pi*f*time_tide)];
Y = tide;
[B] = regress(Y,X);
New = B(1)+B(2)*cos(2*pi*f*time_tide)+B(3)*sin(2*pi*f*time_tide)+B(4)*cos(4*pi*f*time_tide)+B(5)*sin(4*pi*f*time_tide);
tide = tide-New;

[~,ia,ib] = intersect(time,time_tide);
redDailyData.tide = nan(length(time),1);
redDailyData.wl = nan(length(time),1);
redDailyData.seasonal = nan(length(time),1);
redDailyData.msl = nan(length(time),1);
redDailyData.mmsla = nan(length(time),1);
redDailyData.dsla = nan(length(time),1);
redDailyData.tide(ia) = tide(ib);

[~,ia,ib] = intersect(time,time2);

redDailyData.time = time;
%nomeDailyData.swh = combinedWaveTide(:,3);
%nomeDailyData.tp = combinedWaveTide(:,2);
%nomeDailyData.mwd = combinedWaveTide(:,4);
redDailyData.wl(ia) = wl(ib);
redDailyData.seasonal(ia) = climatologyDaily(ib);
redDailyData.msl(ia) = MSL(ib);
redDailyData.mmsla(ia) = MMSLA_hourly(ib);
redDailyData.dsla(ia) = DSLA(ib);
redDailyData.ss = TP_SS;
redDailyData.seasonal_month = climatology;
redDailyData.mmsla_month = MMSLA;
redDailyData.mmsla_month_time = month_time;
redDailyData.B.msl = parameters.msl;
redDailyData.B.season = parameters.seasons;

redDailyData.monthDateVec = datevec(month_time);
redDailyData.hourlyDateVec = datevec(time);


%%


%distance weight (roughly) all tide locations
shish_wl = nome_wl *(170/370) + red_wl*(200/370);
pthope_wl = pru_wl *(160/920) + red_wl*(760/920);
wainright_wl = pru_wl *(175/620) + red_wl*(445/620);
shish_wl2 = nome_wl *(865/1065) + pru_wl*(200/1065);
pthope_wl2 = pru_wl *(430/1190) + nome_wl*(760/1190);
wainright_wl2 = pru_wl *(720/1165) + nome_wl*(445/1165);


ibad = find(isnan(shish_wl) == 1);
shish_wl(ibad) = shish_wl2(ibad);
ibad = find(isnan(pthope_wl) == 1);
pthope_wl(ibad) = pthope_wl2(ibad);
ibad = find(isnan(wainright_wl) == 1);
wainright_wl(ibad) = wainright_wl2(ibad);


%plot tides
figure, 
hold on
plot(times, shish_wl, 'b') 
plot(times, pthope_wl, 'r') 
plot(times, wainright_wl, 'c') 
datetick('x')
legend('Shishmaref', 'Pt Hope', 'Wainright')

save('noaaAlaskaTides.mat')
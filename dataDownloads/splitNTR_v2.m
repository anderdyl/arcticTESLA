function [MSL,MMSL,MMSLA,MMSLA_hourly,DSLA,climatology,climatologyDaily,month_time,parameters] = splitNTR_v2(dat,time,lin,pth,fname)
%[MONTHLYMEAN,MMSLA,SS] = SPLITNTR_v2(DAT,TIME) splits the high and low
% frequency components of the ntr. These components are outputted as:

% MMSL: The monthly mean sea level
% 
% MMSLA: The monthly mean sea level anomaly, computed from the monthly mean
% subtracted from the long term monthly average.
% 
% DSLA: What's left over after the MMSLA and seasonality are removed
%
% climatology: Long term seasonal signal
% 
% climatologyDaily: Long Term seasonal signal at the daily scale

% Created by K.Serafin 12/18/13, Oregon State University
% Updated by K. Serafin 3/31/15, Oregon State University 
%           This update changes MMSLA from a monthly signal to an hourly
%           signal where the data is found using a running average.
% Updated by K. Serafin 9/25/2015, Oregon State University
%           This update changes the way we calculate MSL - from linear
%           regression fit to the hourly time series to a linear regression
%           fit to the annual averages.
% Updated by K.Serafin 6/17/2016, Oregon State University
%           This update changes the linear regression for MSL to occur over
%           the "event" year, e.g., May 2012 - April 2013 = 2012. Also
%           added figures for each of the parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove the mean from the dataset
if exist('fname')
else
    fname = ' ';
end

%% 1.0 Remove mean from dataset

wl = dat;
dat = dat - nanmean(dat);

%% 2. Remove seasons from SWL
[year,month] = datevec(time);

% Calculate long term seasonality 
cnt1 = 0;
climatology = NaN(12,1);
climatologyDaily = NaN(length(year),1);

for ii = 1:12 % For January through February
    cnt1=cnt1+1;
    climatology(cnt1,1) = ii; % 1st column is the month
     ind_good_data = find(month==ii);
    if ~isempty(ind_good_data)>0
        climatology(cnt1,2) = nanmean(dat(month==ii));
        climatology(cnt1,3) = datenum(2000,ii,15,0,0,0);
       climatologyDaily(ind_good_data) =  nanmean(dat(month==ii));
    else
         climatology(cnt1,2)=NaN;
         climatology(cnt1,3) = datenum(2000,ii,15,0,0,0);

         climatologyDaily(ind_good_data) =  NaN;
    end
end

% Find the best fit for the seasonal cycle using annual and semi-annual
% sinusoids.

figure
plot([1:12],[climatology(5:12,2); climatology(1:4,2)],'k-s','markerfacecolor','k')
hold on
set(gca,'fontweight','demi','fontsize',12,'linewidth',2)
grid on
title(['Seasonality ', fname])
ylabel('elevation (m)')
xlim([1 12])
set(gca,'xtick',[1:1:12],'xticklabel',{'M';'J';'J';'A';'S';'O';'N';'D';'J';'F';'M';'A'})
hline(0,'r--')

dateR = datevec(time);
ind = find(dateR(:,3)==15); % To get mid-month results

% Find the best fit to different sinusoids. Stop when we get the lowest
% Root Mean Squared Errors (RMSE)
f = 1/365.25;
cnt = 0;
Y = dat;
rmse = nan(1,4);

for ii = 1:4
   cnt = cnt+1;
  
    if cnt == 1
        X = [ones(length(dat),1) cos(2*pi*f*time)];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time);
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New(ind)).^2)/length(ind));
    elseif cnt == 2
         X = [ones(length(dat),1) cos(2*pi*f*time) sin(2*pi*f*time)];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time)+B(3)*sin(2*pi*f*time);
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New(ind)).^2)/length(ind));
    elseif cnt == 3
         X = [ones(length(dat),1) cos(2*pi*f*time) sin(2*pi*f*time) cos(4*pi*f*time) ];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time)+B(3)*sin(2*pi*f*time)+B(4)*cos(4*pi*f*time);
        
        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New(ind)).^2)/length(ind));
    elseif cnt == 4
        X = [ones(length(dat),1) cos(2*pi*f*time) sin(2*pi*f*time) cos(4*pi*f*time) sin(4*pi*f*time)];
        [B] = regress(Y,X);
        New = B(1)+B(2)*cos(2*pi*f*time)+B(3)*sin(2*pi*f*time)+B(4)*cos(4*pi*f*time)+B(5)*sin(4*pi*f*time);

        rmse(cnt) = sqrt(sum((climatologyDaily(ind)-New(ind)).^2)/length(ind));
    end
  
end

[~,indy] = min(rmse);


if indy == 1
   X = [ones(length(dat),1) cos(2*pi*f*time)];
   [B] = regress(Y,X);
   New = B(1)+B(2)*cos(2*pi*f*time);

elseif indy == 2
    X = [ones(length(dat),1) cos(2*pi*f*time) sin(2*pi*f*time)];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time)+B(3)*sin(2*pi*f*time);

elseif indy == 3
    X = [ones(length(dat),1) cos(2*pi*f*time) sin(2*pi*f*time) cos(4*pi*f*time) ];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time)+B(3)*sin(2*pi*f*time)+B(4)*cos(4*pi*f*time);
        

elseif indy == 4
    X = [ones(length(dat),1) cos(2*pi*f*time) sin(2*pi*f*time) cos(4*pi*f*time) sin(4*pi*f*time)];
    [B] = regress(Y,X);
    New = B(1)+B(2)*cos(2*pi*f*time)+B(3)*sin(2*pi*f*time)+B(4)*cos(4*pi*f*time)+B(5)*sin(4*pi*f*time);
       
end
parameters.seasons = B;

climatologyDaily = New;
ind2 = dateR(:,1)==2020&dateR(:,3)==15&dateR(:,4)==1;
climatology(:,2) = climatologyDaily(ind2);

hold on
plot([1:12],[climatology(5:12,2); climatology(1:4,2)],'r-s','markerfacecolor','k')
hh=legend('Monthly Average','Best-Fit Sinusoids');
set(hh,'location','best')
print(['wave_climatology.png'],'-dpng','-r300')

% print([pth fname,'_climatology.png'],'-dpng','-r300')


%% 3. Compute MMSLA using original dataset - seasonality
dat = wl-climatologyDaily;

% First find the Monthly MSLA one point per month
[year,month] = datevec(time);
cnt = 0;
for ii = year(1):year(end)

    for kk = 1:12
        cnt = 1+cnt;
        ind = find(year==ii&month==kk);
        ind2 = find(isnan(dat(ind))<1);
        if length(ind2)<(28*24*.66) % Has to have at least 66% of the data to count
            MMSLA(cnt) = NaN;
            month_time(cnt) = datenum(ii,kk,15,0,0,0);
        else
            MMSLA(cnt) = nanmean(dat(ind));
            month_time(cnt) = datenum(ii,kk,15,0,0,0);
        end
    end
end

%% 3.0 Compute MSL trend

dat = MMSLA';

if lin>0
    
    % Now fit a trend to it
    ind = find(isnan(dat)<1);
    Y = dat(ind);
    X = [ones(length(dat(ind)),1) month_time(ind)'];
    [B] = regress(Y,X);
    parameters.msl = B;
    MSL = B(1)+B(2)*time;    
    
    figure
    plot(month_time,dat,'k-s','markerfacecolor','k')
    hold on
    plot(time,MSL,'g','linewidth',2)
    datetick
    set(gca,'fontweight','demi','fontsize',12,'linewidth',2)
    grid on
    title(['Long-term Trend ', fname])
    ylabel('elevation (m)')
    xlabel('time')
    print(['wave_longterm.png'],'-dpng','-r300')
%       print([pth fname,'_longterm.png'],'-dpng','-r300')


else
    MSL = nanmean(dat);
end

%% 4.0 Now re-caluclate MMSLA after trend and seasons are removed
% First find the Monthly MSLA one point per month
dat = wl-climatologyDaily-MSL;
[year,month] = datevec(time);
cnt = 0;
for ii = year(1):year(end)

    for kk = 1:12
        cnt = 1+cnt;
        ind = find(year==ii&month==kk);
        ind2 = find(isnan(dat(ind))<1);
        if length(ind2)<(28*24*.66) % Has to have at least 66% of the data to count
            MMSLA(cnt) = NaN;
            month_time(cnt) = datenum(ii,kk,15,0,0,0);
        else
            MMSLA(cnt) = nanmean(dat(ind));
            month_time(cnt) = datenum(ii,kk,15,0,0,0);
        end
    end
end
dateR = datevec(time);

% Also want this monthly anomaly in an hourly time series
newTime = datevec(month_time);
MMSLA_hourly = nan(length(dat),1);
for kk = year(1):year(end)
    for ii = 1:12
        ind = find(newTime(:,1)==kk&newTime(:,2)==ii);
        if isnan(MMSLA(ind))>0
        else
            ind2 = find(dateR(:,1)==kk&dateR(:,2)==ii);
            MMSLA_hourly(ind2) = MMSLA(ind(1));
        end
    end
    
end

figure
plot(month_time,MMSLA,'k','linewidth',2)
hold on
%[y]=moving(MMSLA,12);
[y] = movmean(MMSLA,12);
plot(month_time,y,'m')
set(gca,'fontweight','demi','fontsize',12,'linewidth',2)
grid on
title(['Monthly Mean Anomaly ', fname])
ylabel('elevation (m)')
datetick
%hline(0,'r--')
h=legend('MMA','12-m running avg');
set(h,'location','best')
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[0 0 8 4])
% print([pth fname,'_MMA.png'],'-dpng','-r300')
print(['wave_MMA.png'],'-dpng','-r300')


%% Now find the MMSL by adding climatology to MMSLA
MMSL = MMSLA_hourly+climatologyDaily;


%% DSLA (climatology and MSL already removed)

DSLA = wl-climatologyDaily-MSL-MMSLA_hourly;
%close all

end


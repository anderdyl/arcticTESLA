function [nonTideTime,nonTideData,nonTideDatafilt] = findNTRspectral(tideData,serialTime,bloc,overlap,dt,outdir,tideBands,num,fig_val)
% [NONTIDETIME, NONTIDEDATA, NONTIDEDATAFILT] = FINDNTRSPECTRAL(TIDEDATA,
% SERIALTIME,BLOC,OVERLAP,DT,OUTDIR,TIDEBANDS,NUM,FIG_VAL) uses the 
% Bromirski et al., 2003 approach (with slight modifications) of spectral
% filtering to remove tide bands from hourly storm surge/non-tidal 
% residual time series.
%
%   INPUTS:
%       tideData: input water level time series
%
%       serialTime: input time vector; matlab time
%
%       bloc: blocking time period, in hours. If seasonal and MMSLA are
%           removed, this matters less, otherwise you will filter out 
%           the low frequency signals by blocking.
%
%       overlap: amt of overlap of chunks of time from block
%
%       dt: time step of input time series (e.g., 1 hr)
%
%       outdir: location to save data and figures
%
%       num: number in name of simulation (so you can do multiple in a loop
%           to compare statistics)
%
%       fig_val: 1 = save, 0 or nothing means do not.
%
%
%   OUTPUTS:
%       nonTideTime: time for non-tidal residual time series
%
%       nonTideData: non-tidal residual time series
%
%       nonTideDatafilt: non-tidal residual time series
%
%
% Created by K. Serafin Oregon State University 03/2012
%           last edited by K.Serafin 8/1/2016 (to include more
%           documentation)
%
% Bromirski, P.D., Flick, R.E. and Cayan, D.R., 2003. Storminess 
% variability along the California coast: 1858�2000. Journal of Climate,
% 16(6), pp.982-993.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Does the figure value exist?
if ~exist('fig_val'), fig_val = 0; end

% Get the variance of the original timeseries
vari.original = nanvar(tideData,1);

% Demean and detrend timeseries and recompute variance 
origDemean = tideData-nanmean(tideData);
stats_orig = regstats(origDemean,serialTime,'linear',{'beta','yhat','r'});

tideData = stats_orig.r;
vari.origDeTrend = var(stats_orig.r,1);

% Display average of time series (should be close to zero if mean removed)
nanmean(stats_orig.r)
numCalc = num2str(bloc);

% name directory by the blocking time in hours and go to that directory
mkdir([numCalc,'hr'])
outpath = [outdir, numCalc,'hr/'];

% Find how many blocks there are in the data set
blocOverlap = overlap*bloc;
blocInd = 1:blocOverlap:length(tideData);


% Block data off, tides and serialTime
% Remove the last data block for now...
%tideBlock = cell(length(blocInd)-1,1);
%timeBlock = cell(length(blocInd)-1,1);
for ii = 1:length(blocInd)-1
    if ii < length(blocInd)-1
       tideBlock{ii} = tideData(blocInd(ii):blocInd(ii+2)-1);
       timeBlock{ii} = serialTime(blocInd(ii):blocInd(ii+2)-1); 
    elseif ii == length(blocInd)-1
           if blocInd(ii)+bloc < length(tideData)
               tideBlock{ii} = tideData(blocInd(ii):blocInd(ii)+bloc-1);
               timeBlock{ii} = serialTime(blocInd(ii):blocInd(ii)+bloc-1); 
           else
             
           end
     end
               
end

% Length of any of the tideBlocks will suffice (except maybe last!)
N = length(tideBlock{1});

ts = dt*(0:length(tideBlock{1})-1)';

%% De-Mean, De-Trend & Window All Data

% Find the mean for each cell
tideMean = cellfun(@nanmean,tideBlock);

% Demean the data cell format
tideMeanRM = cellfun(@(a,b) a-b,tideBlock,num2cell(tideMean),'uni',0);

%tideDetrend = cell(length(timeBlock),1);
%STATS = cell(length(timeBlock),1);

% Detrend the data
for ii = 1:length(timeBlock)
    [STATS{ii}] = regstats(tideMeanRM{ii},timeBlock{ii},'linear',{'beta','yhat','r'});
    tideDetrend{ii} = STATS{ii}.r;    % The residuals
end


% Window the data using a tukey window this is to keep the majority of the
% data. We don't care about it being perfect in the frequency domain (where
% other windowing procedures might work better).
window = tukeywin(N,.5);

dataRm = cell(length(timeBlock),1);
yt = cell(length(timeBlock),1);

for ii = 1:length(timeBlock)
    dataRm{ii} = tideDetrend{ii}.*(1 - window);    % This is the data removed by windowing. Must be replaced later.
    yt{ii} = tideDetrend{ii}.*window;           % This is the rest. This what we FFT.
end

%% FFT
% Some general parameters for plotting:
if rem(N,2) == 1
    jj = -(N-1)/2:(N-1)/2;
else
    jj = -N/2:(N/2)-1; 
end
T = N*dt; % # Discrete Samples
Ny = 1/(2*dt); % Nyquest Freq
FF = jj/T;  %Fourier Freq
f = (FF(1):1/(N*dt):FF(end)); 

% Allocate all parameters
Yf = cell(length(timeBlock),1);
Yfc = cell(length(timeBlock),1);
Sj = cell(length(timeBlock),1);
Sjp = cell(length(timeBlock),1);

% Now the data is ready to be tranformed into the frequency domain
for ii = 1:length(timeBlock)
    Yf{ii} = fft(yt{ii},N);
    Yf{ii} = fftshift(Yf{ii})*(1/N);
    Yfc{ii} = conj(Yf{ii}); % Complex conjugate
    
    % Calculate the one sided PSD
    Sj{ii} = (N*dt).*Yfc{ii}.*Yf{ii}; % PSD
    ind_pos = find(f>0);
    Sj{ii}=Sj{ii}(ind_pos);
    cnt1 = 0;
    f_pos = f(ind_pos);
    Sjp{ii} = ones(length(f_pos),1);
    for kk = 1:length(f_pos)
        cnt1=cnt1+1;
    if f_pos(kk)==0 || f_pos(kk)==Ny
        Sjp{ii}(cnt1,:) = Sj{ii}(kk);
    else
        Sjp{ii}(cnt1,:) = 2*Sj{ii}(kk);
    end
    end
end

%% Now begin the process of removing the tidal signal from the data to
% leave only the NTR
% You will be overwriting the fft so save again 
Yforig = Yf;

YftideRM = cell(length(timeBlock),1);
YftideRMc = cell(length(timeBlock),1);
SjtideRM = cell(length(timeBlock),1);
SjptideRM = cell(length(timeBlock),1);


for ii = 1:length(timeBlock)

    for kk = 1:length(tideBands)
        
        % Find the # of spectral estimates surrounding tide band
        SE1 = find(f<tideBands(kk,1)/24);
        SE2 = find(f>tideBands(kk,2)/24);
        % Find the # of spectral estimates in this tideBand
        numSEband = find(f>=(tideBands(kk,1)/24)&f<=(tideBands(kk,2)/24));
        numSEnegband = find(f>=(-tideBands(kk,2)/24)&f<=(-tideBands(kk,1)/24));
        numEsts = round(.10*length(numSEband));
        
        % Now just deal with the positive side and apply this to the
        % negative side
        numSEa = SE1(end-floor(numEsts*.15):end-1):SE1(end-1);
        numSE2 = SE2(1):SE2(floor(numEsts*.85));
        numSEb = [numSEa,numSE2];

        % find a trend across tide band
        [stats] = regstats(Yf{ii}(numSEband),f(numSEband)','linear',{'beta','yhat','r'});

        % Compute the range of variability for the 10 spectral estimates
        allVals = Yf{ii}(numSEb);
  
        varRangea = range(allVals);
         
       
        % randomized variance estimates for real and imaginary
        rR = -1 + (1-(-1)).*rand(length(numSEband),1);
        rI = -1 + (1-(-1)).*rand(length(numSEband),1);
       
        % multiply the variance estimate of each tide band by the randomized
        % estimates
        newB = rR.*real(varRangea); % This is the same as splitting them up
                               %individually and taking magnitude as well
         
       newI = rI.*imag(varRangea);
        newR = complex(newB,newI);
        % Add randomized variance estimates to trend derived estimates.
        aVar = newR+stats.yhat;
        bothVar = aVar;
 
       % Uncomment if you want to plot the frequency bands
        
%         if fig_val ==1
%            % Plot frequency bands
%                 figure(2)
%                 plot(f(numSE)*24,abs(Yf{ii}(numSE)),'r.-','LineWidth',2)
%                 hold on
%                 plot(f(numSE),imag(Yf{ii0}(numSE)),'b.-','LineWidth',2)
%                 plot(f(numSE),bothVar,'m.-')
%                 plot(f(numSE),imag(bothVar),'c.-')
%                 plot(f(numSE10),Yf{ii}(numSE10),'k*')
%                 plot(f(numSE10),imag(Yf{ii}(numSE10)),'g*')
%                 h = legend('real spectra','imaginary spectra','real estimates','imaginary estimates'...
%                     ,'real spectral estimates','imaginary spectral estimates')
%                 set(h,'location','northwest')
%                 set(h,'FontSize',8)
%                 legend boxoff
%                 ylabel('Yf')
%                 xlabel('frequency')
%         else 
%         end

        % Sub this into the FFT at both the neg an pos freq
        YftideRM{ii} = Yf{ii};
        YftideRM{ii}(numSEband) = bothVar;
        YftideRM{ii}(numSEnegband) = conj(flipud(bothVar));
        %YftideRM{ii}(numSEneg) = bothVar;

        YftideRMc{ii} = conj(YftideRM{ii}); % Complex conjugate

        % Calculate the one sided PSD for the new spectrum
  
        SjtideRM{ii} = (N*dt).*YftideRMc{ii}.*YftideRM{ii}; % PSD
        SjptideRM{ii} = SjtideRM{ii}(ind_pos);
        cnt1 = 0;
        for mm = 1:length(f_pos)
            cnt1=cnt1+1;
            if f_pos(mm)==0 || f_pos(mm)==Ny
                SjptideRM{ii}(cnt1,:) = SjptideRM{ii}(mm);
            else
                SjptideRM{ii}(cnt1,:) = 2*SjptideRM{ii}(mm);
            end

        end

    
   
    Yf{ii} = YftideRM{ii};
    end  

    if fig_val == 1
    % This may make a ton of figures so I turn visiability off and just
    % save
        figure('Visible','off')
        semilogy(f_pos*24,Sjp{ii},'b')
        hold on
        semilogy(f_pos*24,SjptideRM{ii},'r')
        ylabel('Amp. Spectrum (m)')%,'fontSize',14,'color','w');
        xlabel('Frequency (cpd)')%,'fontsize',14,'color','w')
        title(['' num2str(ii) ''])
        set(gca,'fontsize',14)

        % Save without viewing
        savename = ['PSD_' num2str(ii), '.png' ''];
        fname_print=[outpath,savename];
        print('-dpng','-r400',fname_print);
    else 
    end
    
    
end

%% Recover Original TimeSeries and with Data Removed
% Shift time series back

yn = cell(length(tideBlock),1);
ynTideRM = cell(length(tideBlock),1);
dataRec = cell(length(tideBlock),1);
dataRectideRM = cell(length(tideBlock),1);
for ii = 1:length(tideBlock)
    Yf{ii} = ifftshift(Yf{ii})*N;
    Yforig{ii} = ifftshift(Yforig{ii})*N;
    yn{ii} = ifft(Yforig{ii});
    ynTideRM{ii} = ifft(Yf{ii});
    
    %% Replace Windowed Data
    dataRec{ii} = yn{ii} + yn{ii}.*(1 - window);
    dataRectideRM{ii} = ynTideRM{ii} + ynTideRM{ii}.*(1 - window);  



end

nonTideD = cell(length(tideBlock),1);
nonTideT = cell(length(tideBlock),1);
TideD = cell(length(tideBlock),1);
% Put final time series back together - getting rid of areas of overlap
for ii = 1 :length(tideBlock)
    nonTideD{ii} = dataRectideRM{ii}(round(bloc/4):bloc-round(bloc/4));
    nonTideT{ii} = timeBlock{ii}(round(bloc/4):bloc-round(bloc/4));
    TideD{ii} = dataRec{ii}(round(bloc/4):bloc-round(bloc/4));
end


% Concatonate all data together
nonTideData = cell2mat(nonTideD);
nonTideData = nonTideData(:);

nonTideTime = cell2mat(nonTideT);
nonTideTime = nonTideTime(:);

withTideData = cell2mat(TideD);
withTideData = withTideData(:);


%Triangular filter data
fwf = [1/3 1/3 1/3];
nonTideDatafilt = conv(nonTideData,fwf);
nonTideDatafilt = nonTideDatafilt(2:end-1);

% See what the variances are
vari.tide = var(withTideData,1);
vari.nontide = var(nonTideData,1);
vari.nontidefilt = var(nonTideDatafilt,1);

close all

% This figure is of the timeseries
if fig_val ==1
for ii = 1:length(tideBlock)
    
    % Plot Original TS windowing mags
    figure('Visible','off')
    subplot 211
    plot(tideDetrend{ii})
    hold on
    plot(dataRec{ii},'r')
    plot(yn{ii},'k')
    xlim([0 bloc])
    title('Original Time Series')
    h1 = legend('Original TS','IFFT w/ window correction','IFFT TS');
    set(h1,'Location','SouthWest')
    legend boxoff
    set(h1,'FontSize',6)
   
    
    subplot 212
    plot(round(bloc/4):bloc-round(bloc/4),tideDetrend{ii}(round(bloc/4):bloc-round(bloc/4)))
    hold on
    plot(round(bloc/4):bloc-round(bloc/4),dataRec{ii}(round(bloc/4):bloc-round(bloc/4)),'r')
    plot(round(bloc/4):bloc-round(bloc/4),yn{ii}(round(bloc/4):bloc-round(bloc/4)),'k')
    title('Removal of 50% of exterior')
    xlim([0 bloc])
    %ylim([-2 4])
    
    % Save without viewing
    savename = ['tideGaugeWindowing_' num2str(ii), '.png' ''];
    fname_print=[outpath,savename];
    print('-dpng','-r300',fname_print);
    
    % Plot NTR windowing mags
    figure('Visible','off')
    subplot 211
    plot(tideDetrend{ii})
    hold on
    plot(1:bloc,dataRectideRM{ii},'r')
    plot(1:bloc,ynTideRM{ii},'k')
    xlim([0 bloc])
    title('Original Time Series')
    h = legend('Original TS','IFFT NTR w/ window correction','IFFT NTR');
    set(h,'Location','SouthWest')
    legend boxoff
    set(h,'FontSize',6)
    
    subplot 212
    plot(round(bloc/4):bloc-round(bloc/4),tideDetrend{ii}(round(bloc/4):bloc-round(bloc/4)))
    hold on
    plot(round(bloc/4):bloc-round(bloc/4),dataRectideRM{ii}(round(bloc/4):bloc-round(bloc/4)),'r')
    plot(round(bloc/4):bloc-round(bloc/4),ynTideRM{ii}(round(bloc/4):bloc-round(bloc/4)),'k')
    title('Removal of 50% of exterior')
    xlim([0 bloc])
    
    savename = ['NTRwindowing_' num2str(ii), '.png' ''];
    fname_print=[outpath,savename];
    print('-dpng','-r300',fname_print);
    
end

else
end


% break the serial time into months
[~,month,~,~,~,~] = datevec(nonTideTime);
wave_month = month;

% Derive monthly means for all data
% vary criteria -> mindata = 13 (40%); mindata = 20 (~66%)

cnt1=0;

ntrClimatology = NaN(12,3);

for ii=min(wave_month):max(wave_month) % For January through February
    cnt1=cnt1+1;
    ntrClimatology(cnt1,1)=ii; % 1st column is the month
    ind_good_data=find(wave_month==ii);
    if ~isempty(ind_good_data)>0
        ntrClimatology(cnt1,2)=length(ind_good_data);
        ntrClimatology(cnt1,3)=nanmean(nonTideDatafilt(find(wave_month==ii)));
    else
         ntrClimatology(cnt1,3)=NaN;
    end
end

output_name = ['SpectralData',num2str(num),'.mat'];

save([outpath output_name], 'ntrClimatology','nonTideDatafilt',...
    'nonTideTime','Sjp','SjptideRM','f_pos',...
    'nonTideData','vari','withTideData','tideBands')

end


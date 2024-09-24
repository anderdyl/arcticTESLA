function tideout = download_noaa_tides_withPred(gauge, datum, startYear, endYear)
%function to download NOAA tides year by year
    %inputs - startYear - 4 digit year to start downloading data from in matlab format
    %         endYear - 4 digit year to end downloading data from in matlab format, it can be the present year
    %         'gauge' - buoy name as a string, e.g., '9435380' for South Beach
    %         'datum' - e.g., 'MHW' or 'NAVD', note that not all stations have all datums (many are missing NAVD)
    
    %initialize variables
    wl = [];
    time = [];
    pred = [];
    time2 = [];
    %loop through each year of data
    for yr = startYear:endYear
        disp(yr)
        %website = ['https://opendap.co-ops.nos.noaa.gov/erddap/tabledap/IOOS_Hourly_Height_Verified_Water_Level.mat?STATION_ID%2CDATUM%2CBEGIN_DATE%2CEND_DATE%2Ctime%2CWL_VALUE&STATION_ID=%22', gauge,'%22&DATUM%3E=%22', datum,'%22&BEGIN_DATE%3E=%22', num2str(yr),'0101%2000%3A00%22&END_DATE%3E=%22', num2str(yr),'1231%2023%3A59%22'];
        %        
        website = append('https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=', num2str(yr),'0101&end_date=', num2str(yr),'1231&station=', gauge,'&product=hourly_height&datum=', datum,'&time_zone=gmt&units=metric&format=csv');
        %disp(website)
        %https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=20200101&end_date=20201231&station=9449424&product=hourly_height&datum=MLLW&time_zone=gmt&units=english&format=csv
        
        out_name = 'tempwaves.csv'; %need to store temporary data to current folder
        try 
            options = weboptions('Timeout',15);
            websave(out_name, website);%,'timeout',10);
            
            fid=fopen(out_name);
    	    header1 = textscan(fid,'%s%s%s%s%s%s%s',1);
            [data2] = textscan(fid,'%16c %f%f%f%f','delimiter',','); 
            fclose(fid);
            
            %load(out_name);
            delete(out_name);
            
            year=str2num(data2{1}(:,1:4));
            month=str2num(data2{1}(:,6:7));
            day=str2num(data2{1}(:,9:10));
            hour=str2num(data2{1}(:,12:13));
            min=str2num(data2{1}(:,15:16));
            sec=zeros(size(min));
    
            datevar = datenum(year,month,day,hour,min,sec);
            
            wl = [wl;(data2{2})];
            time = [time;datevar];

%             load(out_name);
%             delete(out_name);
%             disp(yr)
%             wltemp = IOOS_Hourly_Height_Verified_Wat.WL_VALUE;
%             timetemp = double(IOOS_Hourly_Height_Verified_Wat.time/86400 + datenum(1970,1,1));
% 
%             wl = [wl; wltemp];
%             time = [time; timetemp];
        catch
            continue
        end
         
        %website2 = ['https://opendap.co-ops.nos.noaa.gov/erddap/tabledap/IOOS_Predictions_Water_Level.mat?STATION_ID%2CDATUM%2CBEGIN_DATE%2CEND_DATE%2Ctime%2CWL_VALUE&STATION_ID=%22', gauge,'%22&DATUM%3E=%22', datum,'%22&BEGIN_DATE%3E=%22', num2str(yr),'0101%2000%3A00%22&END_DATE%3E=%22', num2str(yr),'1231%2023%3A59%22'];
        website2 = append('https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=', num2str(yr),'0101&end_date=', num2str(yr),'1231&station=', gauge,'&product=predictions&datum=', datum,'&time_zone=gmt&units=metric&format=csv');
        %website2 = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=', num2str(yr),'0101&end_date=', num2str(yr),'1231&station=', gauge,'&product=predictions&datum=', datum,'&time_zone=gmt&units=metric&format=csv'];

        out_name2 = 'tempwaves2.csv'; %need to store temporary data to current folder
        websave(out_name2, website2);
        
        fid=fopen(out_name2);
    	header1 = textscan(fid,'%s%s%s',1);
        [data] = textscan(fid,'%16c %f','delimiter',','); 
        fclose(fid);
        
        %load(out_name);
        delete(out_name2);
        
        year=str2num(data{1}(:,1:4));
        month=str2num(data{1}(:,6:7));
        day=str2num(data{1}(:,9:10));
        hour=str2num(data{1}(:,12:13));
        min=str2num(data{1}(:,15:16));
        sec=zeros(size(min));

        datevar = datenum(year,month,day,hour,min,sec);
        
        pred = [pred;(data{2})];
        time2 = [time2;datevar];
        
%         
%         predtemp = IOOS_Hourly_Height_Verified_Wat.WL_VALUE;
% 
%         pred = [pred; wltemp];
%         time = [time; timetemp];
%         
        
    end
    
    %output data as structure variable
    tideout.wltime = time;
    tideout.wl = double(wl);
    tideout.predtime = time2;
    tideout.pred = double(pred);
    
end

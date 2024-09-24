function tideout = download_noaa_tides(gauge, datum, startYear, endYear)
%function to download NOAA tides year by year
    %inputs - "startYear" - 4 digit year to start downloading data from in matlab format
    %         "endYear" - 4 digit year to end downloading data from in matlab format, it can be the present year
    %         "gauge" - buoy name as a string, e.g., '9435380' for South Beach
    %         "datum" - e.g., 'MHW' or 'NAVD', note that not all stations have all datums (many are missing NAVD)
    
    %initialize variables
    wl = [];
    time = [];
    
    %loop through each year of data
    for yr = startYear:endYear
        website = ['https://opendap.co-ops.nos.noaa.gov/erddap/tabledap/IOOS_Hourly_Height_Verified_Water_Level.mat?STATION_ID%2CDATUM%2CBEGIN_DATE%2CEND_DATE%2Ctime%2CWL_VALUE&STATION_ID=%22', gauge,'%22&DATUM%3E=%22', datum,'%22&BEGIN_DATE%3E=%22', num2str(yr),'0101%2000%3A00%22&END_DATE%3E=%22', num2str(yr),'1231%2023%3A59%22'];
        out_name = 'tempwaves.mat'; %need to store temporary data to current folder
        websave(out_name, website);
        load(out_name);
        delete(out_name);

        wltemp = IOOS_Hourly_Height_Verified_Wat.WL_VALUE;
        timetemp = double(IOOS_Hourly_Height_Verified_Wat.time/86400 + datenum(1970,1,1));

        wl = [wl; wltemp];
        time = [time; timetemp];
    end
    
    %output data as structure variable
    tideout.mtime = time;
    tideout.wl = double(wl);
    
end

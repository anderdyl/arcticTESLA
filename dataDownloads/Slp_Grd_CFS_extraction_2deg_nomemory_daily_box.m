% The following function extracts CFSR slp fields on a monthly basis and
% manipulates the resulting in fields several ways to check how order and
% resolution affect the results
%


% Method 2
%   extract hourly
%   gradients hourly
%   average 1 day
%   atm predictor daily
%   downscale to 2 degrees


function [SLP_mem_daily_NPSP, GRD_mem_daily_NPSP, TIME_daily, X_in, Y_in] = Slp_Grd_CFS_extraction_2deg_nomemory_daily_box (areaR_lat, areaR_lon, year, month, pathSLP);%,boundx,boundy,ESTELA_obj)



yeart=1996;
montht=11;
% Open netCDF example file
ncid = netcdf.open(fullfile(pathSLP,['prmsl.gdas.' num2str(yeart) num2str(montht,'%02d') '.grb2.nc']),'NOWRITE');


%Latitud
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,4);
% Get variable ID of the first variable, given its name.
varid = netcdf.inqVarID(ncid,varname);
% Get the value of the first variable, given its ID.
latitudCFSR = netcdf.getVar(ncid,varid);

%Longitude
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,5);
% Get variable ID of the first variable, given its name.
varid = netcdf.inqVarID(ncid,varname);
% Get the value of the first variable, given its ID.
longitudCFSR = netcdf.getVar(ncid,varid);

netcdf.close(ncid)

pos_lon1=find(longitudCFSR==areaR_lon(1));
pos_lon2=find(longitudCFSR==areaR_lon(2));
pos_lat1=find(latitudCFSR==areaR_lat(1));
pos_lat2=find(latitudCFSR==areaR_lat(2));

%Latitude - longitude
latitud=latitudCFSR(pos_lat1:pos_lat2);
if areaR_lon(1)>areaR_lon(2)
    longitud=[longitudCFSR(pos_lon1:end); longitudCFSR(1:pos_lon2)];
else
    longitud=longitudCFSR(pos_lon1:pos_lon2);
end

[x,y]=meshgrid(longitud,latitud);
x=double(x);y=double(y);

% Pre-allocating variable names
SLP = [];
SLP_diarias = [];
GRD_diarias = [];
FECHAS_diarias = [];

SLP_B_fullres = [];
GRD_B_fullres = [];
DATES_hourly = [];
SLP_B_1degres = [];
GRD_B_1degres = [];
% Looping through the 2 months to extract hourly full resolution data
display('extracting SLP fields at full resolution')
for mm=1:1
    
%     if mm == 1
%         mon = month -1;
%         if mon == 0
%             mon = 12;
%             yr = year-1;
%         else
%             yr = year;
%         end
%     else
%         mon = month;
%         yr = year;
%     end
%     disp([num2str(yr),'-',num2str(mon)])


    yr = year;
    mon = month;
    disp([num2str(yr),'-',num2str(mon)])
    
    if (yr==2011 & mon>=4) | yr>=2012
        % Open netCDF example file
        ncid = netcdf.open(fullfile(pathSLP,['prmsl.cdas1.' num2str(yr) num2str(mon,'%02d') '.grb2.nc']),'NOWRITE');
    else
        % Open netCDF example file
        ncid = netcdf.open(fullfile(pathSLP,['prmsl.gdas.' num2str(yr) num2str(mon,'%02d') '.grb2.nc']),'NOWRITE');
    end
    
    
    %Time
    % Get name and length of first dimension
    [dimname, dimlen] = netcdf.inqDim(ncid,0);
    % Get the name of the first variable.
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,1);
    % Get variable ID of the first variable, given its name.
    varid = netcdf.inqVarID(ncid,varname);
    % Get the value of the first variable, given its ID.
    time = netcdf.getVar(ncid,varid);
    DATE4=zeros(length(time),4);
    for tt=1:length(time)
        DATE4(tt,1)=str2num(time(1:4,tt)');
        DATE4(tt,2)=str2num(time(5:6,tt)');
        DATE4(tt,3)=str2num(time(7:8,tt)');
        DATE4(tt,4)=str2num(time(9:10,tt)');
    end
    
    

    %Slp
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,6);
    varid = netcdf.inqVarID(ncid,varname);
    if areaR_lon(1)>areaR_lon(2)
        slp1 = netcdf.getVar(ncid,varid,[pos_lon1-1 pos_lat1-1 0],[(length(longitudCFSR)-pos_lon1)+1 (pos_lat2-pos_lat1+1) length(time)]);
        slp2 = netcdf.getVar(ncid,varid,[0 pos_lat1-1 0],[pos_lon2 (pos_lat2-pos_lat1+1) length(time)]);
        slp = [slp1; slp2];
        clear slp1 slp2
    else
        slp = netcdf.getVar(ncid,varid,[pos_lon1-1 pos_lat1-1 0],[(pos_lon2-pos_lon1)+1 (pos_lat2-pos_lat1)+1 length(time)]);
    end    
    netcdf.close(ncid)
    
    [m,n,ntime] = size(slp);

    slp_ = zeros(ntime,m*n);
    
    for tt=1:ntime
        slp_(tt,:) = reshape(squeeze(slp(:,:,tt))',m*n,1);
    end
    
    SLP = [SLP;slp_];
    DATES_hourly = [DATES_hourly;DATE4];    
   
end
    

    display('averaging SLP full resolution to daily')
    c = 1;
    for n = 1:24:length(DATES_hourly(:,3))
        
        DATES_daily(c,:) = DATES_hourly(n,:);
        
        SLP_daily(c,:) = nanmean(SLP(n:n+23,:));
       
        c = c+1;
    end
    
    [mg,ng] = size(x);
    SLP_grid_daily = reshape(SLP_daily,length(DATES_daily),mg,ng);
    
    
    p = projcrs(3411);
    [u,v] = projfwd(p,y,x);


    %[predu,predv] = meshgrid([-3.5e6:0.025e6:1e6],[-0.75e6:0.05e6:4.0e6]);
    %[predu,predv] = meshgrid([-3.5e6:0.05e6:1e6],[-0.75e6:0.05e6:4.0e6]);
    [predu,predv] = meshgrid([-3.5e6:0.05e6:3e6],[-3e6:0.05e6:4.0e6]);

    for uu = 1:length(DATES_daily)
        temp = squeeze(SLP_grid_daily(uu,:,:));
        F = scatteredInterpolant(u(:),v(:),temp(:));
        vq = F(predu,predv);
        SLP_gridded(uu,:) = vq(:);
        
    end
    
    %Vq = interp2(u,v,squeeze(SLP_grid_daily(1,:,:)),predu,predv);
    
    
% 
%     
% %     % Cropping to the bounding area, at 1 degree resolution
%         %   need and index identifying only lat/lons divisible by 1
%         X_ = x(:);
%         Y_ = y(:);
%         
%     c = 1;
%     for ff = 1:length(X_)
%        
%         tempor = rem(X_(ff),1);
%         if tempor == 0
%             tempor2 = rem(Y_(ff),1);
%             if tempor2 == 0
%                 ind2deg(c) = ff;
%                 c = c+1;
%             end
%         end
%     end

% 
%     X_2degres = x(ind2deg);
%     Y_2degres = y(ind2deg);
%     % For Point Hope
%     x_1 = reshape(X_2degres,41,[]);
%     y_1 = reshape(Y_2degres,41,[]); 
%        
%     % For Nags Head
%     x_1 = reshape(X_2degres,56,[]);
%     y_1 = reshape(Y_2degres,56,[]);    
%     % For Nags Head 2nd Attempt
%     x_1 = reshape(X_2degres,67,[]);
%     y_1 = reshape(Y_2degres,67,[]); 
   
   SLP_2degres = SLP_gridded;%daily(:,ind2deg);
   GRD_2degres = gradient_calculator(SLP_gridded,predv);    
   
   X_in = predu;
   Y_in = predv;    
   
   
%    GRD_daily = gradient_calculator(SLP_daily,y_1);    
% 
% 
%     % Cropping to the bounding area, full resolution
%     [in,on] = inpolygon(x_1,y_1,boundx,boundy);
%     SLP_in = SLP_2degres(:,in);
%     GRD_in = GRD_2degres(:,in);
%     X_in = x_1(in);
%     Y_in = y_1(in);
%     
% 
%     % Getting the ESTELA handle
%     load(ESTELA_obj,'C','full','Polar')
%     
%     X_estela = full.Xraw;
%     temp = find(X_estela<0);
%     X_estela(temp) = X_estela(temp)+360;
%     move = X_estela(1:360,:);
%     X_estela(1:360,:) = [];
%     X_estela(361:720,1:317) = move;
%     
%     Y_estela = full.Y;
%     move = Y_estela(1:360,:);
%     Y_estela(1:360,:) = [];
%     Y_estela(361:720,1:317) = move;
%        
% 
%     datetemp = DATES_daily(2:end,2)-DATES_daily(1:end-1,2);
%     if DATES_daily(1,2) == 12
%         sindE = find(datetemp < 0);
%         sindE = sindE+1;
%     else
%         sindE = find(datetemp>0);
%         sindE = sindE+1;
%     end
%     
%     TIME_daily = DATES_daily(sindE:end,:);
% 
% 
%     display('interpolating travel times for daily resolution')
% 
%     
%     for ff = sindE:length(DATES_daily(:,1));
%         
%         monE = DATES_daily(ff,2);
%         
%         if monE <=2 || monE == 12
%             W = C.traveldays_interp.DJF;
%         elseif monE >= 3 && monE <=5
%             W = C.traveldays_interp.MAM;
%         elseif monE >= 6 && monE <=8
%             W = C.traveldays_interp.JJA;
%         elseif monE >= 9 && monE <=11
%             W = C.traveldays_interp.SON;
%         end
%         
%         move = W(1:360,:);
%         W(1:360,:) = [];
%         W(361:720,1:317) = move;
%         temp = ceil(flipud(W')*1);
%         temp_interp = interp2(X_estela',flipud(Y_estela'),temp,x_1,y_1);
%         
%         temp_interp = temp_interp(in(:));
%         
%         xvector = x_1(in(:));
%         yvector = y_1(in(:));
%         travelv = temp_interp;
%         
%         madeup_slp = NaN.*ones(length(travelv),1);
%         madeup_grdslp = NaN.*ones(length(travelv),1);
%         
%         for hh = 1:25
%             index = find(temp_interp == hh);
%             madeup_slp(index) = SLP_in(ff-hh,index);
%             madeup_grdslp(index) = GRD_in(ff-hh,index);
%         end
%         
%         index = find(temp_interp > hh);
%             madeup_slp(index) = SLP_in(ff-hh,index);
%             madeup_grdslp(index) = GRD_in(ff-hh,index);        
%         
%         [m_regional,n_regional]=size(x_1);
%         
%         modo_slp = ones(m_regional*n_regional,1)*NaN;
%         modo_slp(in(:)) = madeup_slp;
%         modo_grdslp = ones(m_regional*n_regional,1)*NaN;
%         modo_grdslp(in(:)) = madeup_grdslp;
%         reshaped_slp = reshape(modo_slp,m_regional,n_regional);
%         reshaped_grdslp = reshape(modo_grdslp,m_regional,n_regional);
%         
%         slp_mem_daily(ff,:) = reshaped_slp(in(:));
%         grdslp_mem_daily(ff,:) = reshaped_grdslp(in(:));
%         
%     end
%       
%     
%     
%     %daily_TIME = DATES_hourly(sindE:end,:);
%     %SLP_1degdaily_NPSP  = SLP_1deg_daily(sindE:end,:);
%     %GRD_1degdaily_NPSP  = GRD_1deg_daily(sindE:end,:);
%     SLP_mem_daily_NPSP  = slp_mem_daily(sindE:end,:);
%     GRD_mem_daily_NPSP  = grdslp_mem_daily(sindE:end,:);    
    SLP_mem_daily_NPSP  = SLP_2degres(1:end,:);
    GRD_mem_daily_NPSP  = GRD_2degres(1:end,:);    
    TIME_daily = DATES_daily;
    
    
    



    
    
    





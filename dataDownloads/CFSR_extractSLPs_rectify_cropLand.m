% CFSR_extractSLPs_master_script.m
%   created by D. Anderson 11/2017

% The following is a wrapping script used to extract SLPs from CFSR within
% a bounding region. 

% Needs following sub-functions:
%   Slp_Grd_CFS_extraction_2deg_daily.m 
%   m_map toolbox
%   gradient_calculator.m

clear
% Point Hope
areaR_lat = [90; 44.0];
areaR_lon = [0; 359];  
ESTELA_obj = 'pointHopeTest.mat';

pathSLP = '/media/dylananderson/Elements/CFS/prmsl';

slp_mem = [];
grd_mem = [];
time = [];
c = 1;

for year = 1979:2023;
    if year == 1979
        for month = 2:12;
            
            [SLP_mem_daily_NPSP, GRD_mem_daily_NPSP, TIME_daily, X_in, Y_in] = Slp_Grd_CFS_extraction_2deg_nomemory_daily_box (areaR_lat, areaR_lon, year, month, pathSLP);%,boundx,boundy,ESTELA_obj);

            slp_mem = [slp_mem;SLP_mem_daily_NPSP];
            grd_mem = [grd_mem;GRD_mem_daily_NPSP];
            time = [time;TIME_daily];
        end
    elseif year == 2023
        for month = 1:5;

            [SLP_mem_daily_NPSP, GRD_mem_daily_NPSP, TIME_daily, X_in, Y_in] = Slp_Grd_CFS_extraction_2deg_nomemory_daily_box (areaR_lat, areaR_lon, year, month, pathSLP);%,boundx,boundy,ESTELA_obj);

            slp_mem = [slp_mem;SLP_mem_daily_NPSP];
            grd_mem = [grd_mem;GRD_mem_daily_NPSP];
            time = [time;TIME_daily];           
            
        end        
    else
        for month = 1:12;

            [SLP_mem_daily_NPSP, GRD_mem_daily_NPSP, TIME_daily, X_in, Y_in] = Slp_Grd_CFS_extraction_2deg_nomemory_daily_box (areaR_lat, areaR_lon, year, month, pathSLP);%,boundx,boundy,ESTELA_obj);

            slp_mem = [slp_mem;SLP_mem_daily_NPSP];
            grd_mem = [grd_mem;GRD_mem_daily_NPSP];
            time = [time;TIME_daily];
        end        
    end
end


%%


[mq,nq] = size(X_in);

load('paleta2.mat');
figure
%wt = ones(mq*nq,1)*NaN;
temp = SLP_mem_daily_NPSP(2,:)./100;
%wt(sea_sq) =  temp;
reshaped = reshape(temp,mq,nq);
iso1_N = 970:2:1045;

p = projcrs(3411);
%[u,v] = projfwd(p,y,x);
[lat,lon] = projinv(p,X_in,Y_in);

% m_proj('stereo','lat',[-65 65],'lon',[120 280]);
m_proj('stereographic','lat',90,'long',-45,'radius',40);
m_coast('patch',[0.5 0.5 0.5]);
m_grid('box','fancy','xtick',[],'ytick',[])
hold on
[cs,h] = m_contourf(lon,lat,reshaped,iso1_N,'linewidth',0.5,'linecolor','none');

m_coast('patch',[0.5 0.5 0.5]);
%%

figure
contourf(X_in,Y_in,reshaped,'linewidth',0.5,'linecolor','none')
hold on
load('coastlines')
p = projcrs(3411);
[u_coast,v_coast] = projfwd(p,coastlat,coastlon);
plot(u_coast,v_coast,'k-','linewidth',2)
%
% % Point hope
% boundx=[-2e6,-1.75e6,-1.45e6,-1.25e6,-1e6,-8e5,-4.5e5,-2e5,1.5e5,6e5,...
%     8e5,9.5e5,1e6,8e5,5.5e5,2.5e5,-3e5,-7e5,-1.1e6,-1.6e6,-1.8e6,...
%     -2e6,-1.95e6,-1.7e6,-1.85e6,-2.15e6,-2.5e6,-2.85e6,-3.2e6,...
%     -3.35e6,-3.45e6,-3.5e6,-3.35e6,-3.25e6,-3e6,-2.7e6,-2.55e6,...
%     -2.55e6,-2.3e6,-2e6,-1.9e6,-2.05e6,-2.25e6,-2.1e6];
% boundy=[-2.5e5,-3e5,-2.5e5,-3.5e5,-5e5,-6e5,-5.5e5,-7e5,-6e5,-6e5,...
%     -4.5e5,-5e4,5.5e5,1e6,1.4e6,1.7e6,1.8e6,1.9e6,1.85e6,1.7e6,1.8e6,...
%     2.15e6,2.7e6,3.1e6,3.55e6,3.7e6,3.7e6,3.65e6,3.45e6,...
%     3.1e6,2.7e6,2.05e6,1.65e6,1.4e6,1.65e6,1.65e6,1.5e6,...
%     1.3e6,1.05e6,9.5e5,6.5e5,4e5,1e5,-1.5e5];

%% All Arctic
boundx=[-2.299e6,-1.75e6,-1.45e6,-1.25e6,-1e6,-8e5,-4.5e5,-2e5,1.5e5,6e5,1.1e6,1.6e6...
    2.05e6,2.45e6,2.15e6,1.4e6,8.5e5,2.5e5,-3e5,-7e5,-1.1e6,-1.6e6,-1.65e6,...
    -1.95e6,-1.95e6,-1.7e6,-1.85e6,-2.15e6,-2.5e6,-2.85e6,-3.2e6,...
    -3.35e6,-3.5e6,-3.5e6,-3.45e6,-3.25e6,-3e6,-2.7e6,...
    -2.65e6,-2.3e6,-2.2e6,-2.25e6,-2.35e6,-2.5e6,-2.45e6,-2.3e6];
boundy=[-0.5e6,-0.85e6,-0.95e6,-1.0e6,-1.1e6,-1.11e6,-1.08e6,-1.05e6,-1.08e6,-1.12e6,-1.15e6,-0.95e6...
    -0.75e6,-5e4,1.05e6,1.45e6,1.7e6,2.05e6,2.1e6,2.2e6,2.15e6,2.0e6,2.05e6,...
    2.25e6,2.7e6,3.1e6,3.55e6,3.7e6,3.7e6,3.65e6,3.45e6,...
    3.1e6,2.7e6,2.05e6,1.65e6,1.3e6,1.5e6,1.55e6,...
    1.3e6,1.05e6,9.5e5,6.5e5,4e5,1e5,-1.5e5,-0.49e6];

plot(boundx,boundy,'r-','linewidth',1.5)

xlim([-3.5e6, 3e6])
ylim([-3e6,4e6])

%%
[in,on] = inpolygon(X_in,Y_in,boundx,boundy);
% % SLP_in = SLP_2degres(:,in);
% % GRD_in = GRD_2degres(:,in);
% % X_in = x_1(in);
% % Y_in = y_1(in);
% 
% plot(X_in(in),Y_in(in),'ro')
% plot(X_in(on),Y_in(on),'bo')

Xsea = X_in(in);
Ysea = Y_in(in);
SLPsea = slp_mem(:,in);
GRDsea = grd_mem(:,in);

%%
[mq,nq] = size(X_in);

load('paleta2.mat');
figure
%wt = ones(mq*nq,1)*NaN;
temp = SLP_mem_daily_NPSP(2,:)./100;
%wt(sea_sq) =  temp;
reshaped = reshape(temp,mq,nq);
iso1_N = 970:2:1045;

p = projcrs(3411);
%[u,v] = projfwd(p,y,x);
[lat,lon] = projinv(p,X_in,Y_in);

[boundxLat,boundyLon] = projinv(p,boundx,boundy);

% m_proj('stereo','lat',[-65 65],'lon',[120 280]);
m_proj('stereographic','lat',90,'long',-45,'radius',45);
m_coast('patch',[0.5 0.5 0.5]);
m_grid('box','fancy','xtick',[],'ytick',[])
hold on
[cs,h] = m_contourf(lon,lat,reshaped,iso1_N,'linewidth',0.5,'linecolor','none');

m_coast('patch',[0.5 0.5 0.5]);
m_plot(boundyLon,boundxLat,'ro-')

 
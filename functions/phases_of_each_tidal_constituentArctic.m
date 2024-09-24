addpath(genpath('/Users/dylananderson/Documents/projects/arcticClimate/'))
addpath(genpath('/Users/dylananderson/Documents/data/chsWest/'))
%addpath(genpath('C:/Users/anderdyl/Dropbox/downloaded_functions/'))
% addpath(genpath('/media/dylananderson/Elements/downloaded_functions/'))
%addpath(genpath('/home/server/pi/homes/danderso/SERDP/'))
%addpath(genpath('/home/server/pi/homes/danderso/Research/'))
addpath(genpath('/Users/dylananderson/Documents/projects/serdp/'))


clear

% load('frfTideFitting.mat')
load('arcticTides.mat')

tide = purdoeDailyData.tide;
time = purdoeDailyData.time;

% for fff = 1:length(timeStampHours)
% 
%     time(fff) = datenum([1970,1,1,0,0,timeStampHours(fff)*(60*60)]);
% 
% end


hourlyTime = time;%[datenum(1999,1,1,0,0,0):1/24:datenum(2019,12,31,23,0,0)];
hourlyTide = tide;%interp1q(time',predTide',hourlyTime');


ind = find(hourlyTime >= datenum(1990,1,1,0,0,0) & hourlyTime <= datenum(2024,12,31,23,0,0));
tide = hourlyTide(ind);
time = hourlyTime(ind)+4/24;







% %load('dailyData_SanDiego.mat')
% %load('LaJolla_daily.mat')
% %plot(dailyData.time,dailyData.tide)
% ind = find(dailyData.time >= datenum(1980,1,1,0,0,0) & dailyData.time <= datenum(1998,8,7,1,0,0));
% tide = dailyData.tide(ind);
% time = dailyData.time(ind);

[nameu,fu,tidecon,xout]=t_tide(tide,'start time',[1998 1 1 0 0 0],'latitude',65.5);

[sorted,sortind] = sort(tidecon(:,1),'descend');
sortednames = nameu(sortind,:);
sortedtidecon = tidecon(sortind,:);
sortedfu = fu(sortind,:);
% M2 = sortedtidecon(1,:);
% S2 = sortedtidecon(3,:);
% N2 = sortedtidecon(5,:);
% K2 = sortedtidecon(7,:);
% K1 = sortedtidecon(2,:);
% O1 = sortedtidecon(4,:);
% P1 = sortedtidecon(6,:);
% Q1 = sortedtidecon(8,:);

consts = 8;%64;
%time_emulator = [datenum(1700,6,1,0,0+ffff*2,0):1/24:datenum(2800,5,31,23,0+ffff*2,0)];
%time_emulator = [datenum(2018,1,1,0,0,0):1/24:datenum(2018,5,31,23,0,0)];
time_emulator = [datenum(1979,1,1,0,0,0):1/24:datenum(2122,6,1,0,0,0)];
%time_emulator = [datenum(1980,1,1,0,0,0):1/24:datenum(2013,12,31,23,0,0)];



yout=t_predic(time_emulator,sortednames(1:consts,:),sortedfu(1:consts,:),sortedtidecon(1:consts,:),'latitude',65.5);

% 
% 
% index = find(time_emulator > datenum(1930,1,1,0,0,0) & time_emulator < datenum(2010,1,1,0,0,0));
% 
% index2 = find(dailyData.time > datenum(1930,1,1,0,0,0) & dailyData.time <= datenum(2010,1,1,0,0,0));
% figure
% p9 = plot(dailyData.time(index2), dailyData.tide(index2)-yout(index)')
% hold on
% st = dailyData.time(index2);
% subset =  dailyData.tide(index2)-yout(index)';
% c = 1;
% for qqq = 1:12*40;
%    
%     sst = st(c:c+60*24);
%     ssub = subset(c:c+60*24);
%     ind = find(ssub == max(ssub));
%     
%     p10 = plot(sst(ind),ssub(ind),'ro');
%     store_t(qqq) = sst(ind);
%     store_wl(qqq) = ssub(ind);
%     
%     c = c + 60*24;
% end
% 
% 
% 
% xmod = [store_t'];
% 
% y = store_wl';
% 
% modelfun = @(b,x) b(1) + b(2)*cos(2*pi/(18.6*365.25)*xmod(:,1) +b(3));
% 
% 
% beta = [.34 .025 0];
% 
% mdl = fitnlm(xmod,y,modelfun,beta);
% coef = table2array(mdl.Coefficients);
% 
% 
% %figure
% col = rgb('black');
% p11 = plot(store_t,mdl.Fitted,'color',col,'linewidth',2);
% datetick 
% legend([p9 p10 p11],'Residuals','2 month maximum','sine fit')
% 
% %xlim([t(1) t(end)])
% ylabel('monthly sea level (mm)')
% title('San Diego','color','k')
% 

%%
figure(1)
%plot(time,tide)
%hold on
plot(time_emulator,yout)
hold on

[y] = runmean(yout,24*183);

figure(10)
%plot(time,tide)
%hold on
plot(time_emulator,y)
hold on


[f,Sf] = FFTRaw(y(1:365.25*19*24),1); % Can also look at whole time series

figure
semilogy(f*24,Sf,'b')
xlabel('cycles per day')




jdmid=mean(time_emulator(1:2*fix((length(time_emulator)-1)/2)+1));

time2 = time_emulator-jdmid;

% % My own predictions
% for qq = 1:consts
%     tides(qq,1:length(time)) = sortedtidecon(qq,1).*cos(2*pi*24*sortedfu(qq).*(time2) + sortedtidecon(qq,3)*pi/180);
% %        tides(qq,1:length(time)) = sortedtidecon(qq,1).*cos(2*pi*24*sortedfu(qq).*(time2-9/24) + sortedtidecon(qq,3)*pi/180);  this works for 1 through 3?
% end
% 
% 
% tsum = sum(tides);
% plot(time,tsum,'k','linewidth',1.5)
% xlim([datenum(1980,1,1,0,0,0) datenum(1980,2,1,0,0,0)])
% 


  ap=sortedtidecon(1:consts,1)/2.*exp(-i*sortedtidecon(1:consts,3)*pi/180);
  am=conj(ap);
  % jdmid=mean(tim(1:2*fix((length(tim)-1)/2)+1));
tim = time_emulator;
freq = sortedfu(1:consts);
    const=t_getconsts;
  ju=zeros(size(sortedfu(1:consts)));

  % Check to make sure names and frequencies match expected values.

  for k=1:size(sortednames(1:consts,:),1),
    ju(k)=strmatch(sortednames(k,:),const.name);
  end;
  %if any(freq~=const.freq(ju)),
  %  error('Frequencies do not match names in input');
  %end;
  lat = 32.5;
  ltype='nodal';

% Get the astronical argument with or without nodal corrections.  
if ~isempty(lat) & abs(jdmid)>1,				  
  [v,u,f]=t_vuf(ltype,jdmid,ju,lat);				  
elseif abs(jdmid)>1, % a real date				  
  [v,u,f]=t_vuf(ltype,jdmid,ju);				  
else								  
   v=zeros(length(ju),1);					  
   u=v; 							  
   f=ones(length(ju),1);					  
end;								  


ap2=ap.*f.*exp(+i*2*pi*(u+v));
am2=am.*f.*exp(-i*2*pi*(u+v));

tim=tim-jdmid;

[n,m]=size(tim);
tim=tim(:)';
ntim=length(tim);

nsub=10000; % longer than one year hourly.
for j1=1:nsub:ntim
  j2=min(j1 + nsub - 1,ntim);
  yout(j1:j2)=sum(exp( i*2*pi*freq*tim(j1:j2)*24).*ap2(:,ones(1,j2-j1+1)),1)+ ...
              sum(exp(-i*2*pi*freq*tim(j1:j2)*24).*am2(:,ones(1,j2-j1+1)),1);
end;
yout=reshape(yout,n,m);

plot(time_emulator,yout)
  
for ff = 1:consts;
    
    newy(ff,1:length(time_emulator)) = exp(i*2*pi*freq(ff).*time2*24).*ap2(ff) + exp(-i*2*pi*freq(ff).*time2*24).*am2(ff);
    
end

tsum2 = sum(newy);
%plot(time_emulator,tsum2,'m','linewidth',1)



%%

yall = newy;%[M2;S2;N2;K2;K1;O1;P1;Q1];

for qqq = 1:consts

    y = yall(qqq,1:200);
    x = time_emulator(1:200);
    
    
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of �y�
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset
%modelfun = @(b,x)  b(1).*(cos(2*pi*24.*x.*b(2) + b(3)*pi/180)) + b(4);    % Function to fit
fit = @(b,x)  sortedtidecon(qqq,1).*(cos(2*pi*24.*x.*sortedfu(qqq) + b(1)*pi/180));% + b(4);    % Function to fit

fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
 s2 = fminsearch(fcn, [sortedtidecon(qqq,3)*pi/180]);%;  ym])     % Minimise Least-Squares
%beta0 = [sortedtidecon(qqq,1) sortedfu(qqq) 0 0];
%s2 = nlinfit(x,y(qqq,:),modelfun,beta0);

phases(qqq) = s2;
fitted =  sortedtidecon(qqq,1).*(cos(2*pi*24.*x.*sortedfu(qqq) + s2(1)*pi/180));
% figure(qqq)
% plot(x,y,'b', x, fitted, 'r')
clear yu yl yr yz zx per ym fit fcn s2
% grid
end

%%

% My own predictions
for qq = 1:consts
    tides(qq,1:length(time_emulator)) = sortedtidecon(qq,1).*cos(2*pi*24*sortedfu(qq).*(time_emulator) + phases(qq)*pi/180);
%        tides(qq,1:length(time)) = sortedtidecon(qq,1).*cos(2*pi*24*sortedfu(qq).*(time2-9/24) + sortedtidecon(qq,3)*pi/180);  this works for 1 through 3?
end

figure(1)
tsum = sum(tides);
plot(time_emulator,tsum,'k','linewidth',1.5)
xlim([datenum(1980,1,1,0,0,0) datenum(1980,2,1,0,0,0)])


tideOut = tsum;
indexEmulate = find(time_emulator >= datenum(2022,3,1,0,0,0));

synTime = time_emulator(indexEmulate);
synTide = tsum(indexEmulate);

%{

%%
clearvars -except newy phases sortedfu sortedtidecon time_emulator ffff
close all
% M2 = newy(1,:);
% S2 = newy(3,:);
% N2 = newy(5,:);
% K2 = newy(7,:);
% K1 = newy(2,:);
% O1 = newy(4,:);
% P1 = newy(6,:);
% Q1 = newy(8,:);
%load('dailyData_SanDiego.mat')
%load('LaJolla_daily.mat')
time = time_emulator;%dailyData.time;
M2phase = phases(1);
S2phase = phases(3);
N2phase = phases(2);
%K2phase = phases(7);
K1phase = phases(4);
%O1phase = phases(4);
%P1phase = phases(6);
%Q1phase = phases(8);

M2freq = 24*sortedfu(1);
S2freq = 24*sortedfu(3);
N2freq = 24*sortedfu(2);
%K2freq = 24*sortedfu(7);
K1freq = 24*sortedfu(4);
%O1freq = 24*sortedfu(4);
%P1freq = 24*sortedfu(6);
%Q1freq = 24*sortedfu(8);

M2amp = sortedtidecon(1,1);
S2amp = sortedtidecon(3,1);
N2amp = sortedtidecon(2,1);
%K2amp = sortedtidecon(7,1);
K1amp = sortedtidecon(4,1);
%O1amp = sortedtidecon(4,1);
%P1amp = sortedtidecon(8,1);
%Q1amp = sortedtidecon(8,1);

M2 = M2amp.*cos(2*pi*M2freq.*time_emulator + M2phase*pi/180);
S2 = S2amp.*cos(2*pi*S2freq.*time_emulator + S2phase*pi/180);
N2 = N2amp.*cos(2*pi*N2freq.*time_emulator + N2phase*pi/180);
%K2 = K2amp.*cos(2*pi*K2freq.*time_emulator + K2phase*pi/180);
K1 = K1amp.*cos(2*pi*K1freq.*time_emulator + K1phase*pi/180);
%O1 = O1amp.*cos(2*pi*O1freq.*time_emulator + O1phase*pi/180);
%P1 = P1amp.*cos(2*pi*P1freq.*time_emulator + P1phase*pi/180);
%Q1 = Q1amp.*cos(2*pi*Q1freq.*time_emulator + Q1phase*pi/180);

all = [M2;S2;N2;K1];%;K2;K1;O1;P1;Q1];
fourtides = sum(all);
allamp = [M2amp;S2amp;N2amp;K1amp;];%K2amp;K1amp;O1amp;P1amp;Q1amp];

figure
%plot(time_emulator,dailyData.tide)
%hold on
plot(time_emulator,fourtides)

% %%
% allshort = all(:,1:end);
% shorttime = time_emulator(1:end);
% 
% for qqq = 1:8
%     
% %     figure(qqq)
% 
% dphase = allshort(qqq,2:end)-allshort(qqq,1:end-1);
% Totphi = acos(allshort(qqq,:)/allamp(qqq))*180/pi;
% dphidt = Totphi(2:end)-Totphi(1:end-1);
% 
% 
% 
% 
% 
% if qqq == 1 || qqq == 2 || qqq == 3 || qqq == 4
% sind = find(abs(dphidt) < 27 & Totphi(2:end) > 90);
% pind = find(dphase >= 0 & abs(dphidt) > 27);
% elseif qqq == 5 
% sind = find(abs(dphidt) < 10 & Totphi(2:end) > 90);
% pind = find(dphase >= 0 & abs(dphidt) > 10);
% elseif qqq == 6 
% sind = find(abs(dphidt) < 10 & Totphi(2:end) > 90);
% pind = find(dphase >= 0 & abs(dphidt) > 10);
% elseif qqq == 7 
% sind = find(abs(dphidt) < 10 & Totphi(2:end) > 90);
% pind = find(dphase >= 0 & abs(dphidt) > 10);
% elseif qqq == 8
% sind = find(abs(dphidt) < 10 & Totphi(2:end) > 90);
% pind = find(dphase >= 0 & abs(dphidt) > 10);
% end
% 
% index = find(dphidt <= 0);
% Totphi(pind+1) = 360-Totphi(pind+1);%(180-Totphi(index)) + 180;%
% Totphi(sind+1) = 360-Totphi(sind+1);
% % %plot(shorttime(2:end),dphidt,'-o')
% % subplot(211)
% % plot(shorttime,allshort(qqq,:))
% % hold on
% % %plot(shorttime(pind+1),allshort(qqq,pind+1),'ro')
% % subplot(212)
% %  plot(shorttime,Totphi,'-o')
% %  hold on
% % plot(shorttime(sind+1),Totphi(sind+1),'ro')
% % plot(shorttime(index),Totphi(index),'k.')
% degrees(qqq,1:length(Totphi)) = Totphi;
% 
% end
% 
% 
% %save(['tide_emulation_',num2str(ffff),'.mat'],'degrees','time_emulator')
% %close all
% 
% 
% 


%}
save('tide_emulation_Purdoe.mat','time_emulator','tideOut','synTime','synTide')

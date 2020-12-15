%% JRA55 forcing data
clear all;
close all;
clc;
%%
maindir='G:\JRA55\';

history_file_path=[maindir];

history_file_name=['JRA55_03hr_forcing_2005.nc'];
history_file=[history_file_path, history_file_name];

% diary([history_file_path,'history_file_info.txt'])% 输出重定向
% ncdisp(history_file);%查看nc格式数据文件 描述信息
% diary off %关闭输出重定向

% LON =ncread(history_file,'LON');% ice area  (aggregate)（无量纲）
% 
% LAT =ncread(history_file,'LAT');% ice area  (aggregate)（无量纲）
% 
% time   =ncread(history_file,'time');% ice area  (aggregate)（无量纲）

airtmp =ncread(history_file,'airtmp');% ice area  (aggregate)（无量纲）

spchmd =ncread(history_file,'spchmd');% ice area  (aggregate)（无量纲）

wndewd =ncread(history_file,'wndewd');% ice area  (aggregate)（无量纲）

wndnwd =ncread(history_file,'wndnwd');% ice area  (aggregate)（无量纲）

glbrad =ncread(history_file,'glbrad');% ice area  (aggregate)（无量纲）

dlwsfc =ncread(history_file,'dlwsfc');% ice area  (aggregate)（无量纲）

ttlpcp =ncread(history_file,'ttlpcp');% ice area  (aggregate)（无量纲）

% save -v7.3 JRA55_03hr_forcing_2005.mat LON LAT time airtmp spchmd wndewd wndnwd ttlpcp glbrad  dlwsfc
 save -v7.3 JRA55_03hr_forcing_2005_lonlattime.mat LON LAT time 
%%
% load H:\程序\JRA55_03hr_forcing_2005.mat;
load  JRA55_03hr_forcing_2005_lonlattime.mat;
%%
LONLON=LON(:);
lon_goal =220;
[sLON,index_lon]= sort(abs(LONLON - lon_goal) );
LATLAT=LAT(:);
lat_goal =70 ;
[sLAT,index_lat]= sort(abs(LATLAT - lat_goal) );

% enhance roubusty and compatible
% use more variables to replace consts
for n=1:1000
    
    
    % n =2 ;
    
    lon_find_alone= LON(index_lon(1:n));
    
    
    % n =2 ;
    
    lat_find_alone= LAT(index_lat(1:n));
    
    % n=1000;
    [RESULT,ia,ib] = intersect(index_lat(1:n), index_lon(1:n));
    if ~isempty(RESULT)
        
        break;
    end
    disp(['n=',num2str(n)])
end

lon_find =LON (RESULT(1))
lat_find =LAT (RESULT(1))

index_lon_find=index_lon(ib)
index_lat_find=index_lat(ia)
[i, j] =ind2sub([320,384],RESULT(1))
% LON (i,j )
% LAT (i,j )

%% airtmp
var = airtmp(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('air temperature','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('air temperature','fontsize',18)
%%
t_1hr =var_1hr;

%% spchmd
var = spchmd(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('Specific Humidity','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('Specific Humidity','fontsize',18)
%%
q_1hr =var_1hr;

%% uwind
var = wndewd(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('u wind','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('u wind','fontsize',18)
%%
uwind_1hr =var_1hr;

%% vwind
var = wndnwd(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('vwind','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('vwind','fontsize',18)
%%
vwind_1hr =var_1hr;

%% longwave
var = dlwsfc(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('longwave','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('longwave','fontsize',18)
%%
longwave_1hr =var_1hr;

%% shortwave
var = glbrad(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('shortwave','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('shortwave','fontsize',18)
%%
shortwave_1hr =var_1hr;


%% precipitation
var = ttlpcp(i,j,:);
var_1p= squeeze( var );
%%
figure
plot(var_1p,'linewidth',2)
xlabel('time(03hr spacing)','fontsize',18)
title('precipitation','fontsize',18)
set(gca,'fontsize',16)
%% 插值
x=0:3:8760-3; x=x';
xq=0:1:8760-1; xq = xq';
var_1hr = interp1(x, var_1p, xq);

%% 补NAN
var_1hr(8759) = var_1hr(8758);
var_1hr(8760) = var_1hr(8758);
%%
figure
plot(var_1hr,'linewidth',2)
xlabel('time(01hr spacing)','fontsize',18)
title('precipitation','fontsize',18)
%%
precipitation_1hr =var_1hr;

%%
% data= 0;
% data=cat(2, shortwave_1hr, longwave_1hr, uwind_1hr, vwind_1hr, t_1hr, q_1hr, precipitation_1hr);
% save jra55_2005_220_70_01hr.mat data ;
% data= roundn(data, -8)














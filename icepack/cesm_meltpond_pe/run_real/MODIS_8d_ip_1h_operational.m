%% MODIS_8d_ip_1h_operational.m
% 从MODIS的8d文件生成1小时文件
% 20200406，luyang,v1
%%
clc;
close all;
clear;
%% 导入原始数据
year = 2005; lon = 90; lat = 78 ;
char_data_year = num2str(year);
char_lon = num2str(lon);
char_lat = num2str(lat);
       

filename=['./','modis_',char_data_year,'_',char_lon,'_',char_lat,'_08day.txt'];

fid =fopen(filename);

formatSpec = '%s';
N = 2;
headline = textscan(fid, formatSpec, N );

formatSpec = '%d %f';
A = textscan(fid, formatSpec );

fclose(fid);


% day_obs = A{1,1};
% mpf_obs = zeros(1,365) ;

for i = 1:16
%     dayindex = day_obs(i);
%     mpf_obs(dayindex) = A{1,2}(i,1);
    mpf_obs(i) = A{1,2}(i,1);
end
mpf_obs(mpf_obs == 103) = nan;



%% 插值到1小时
hour = 3096:1:5976;
hour = hour';

mpf_obs( mpf_obs==103 ) = nan; 

x = 3096:192:5976;
xq =  3096:1:5976;
mpf_ip_1h = interp1( x ,mpf_obs, xq ,'spline') ;



mpf_ip_1h = mpf_ip_1h';
mpf_ip_1h( mpf_ip_1h < 0) =0;


data =cat (2, hour, mpf_ip_1h);


%% 保存为txt文件


filename=['./','modis_',char_data_year,'_',char_lon,'_',char_lat,'_01hour_ip.txt'];

fid =fopen(filename,'wt');
formatSpec = '%9s\t';
fprintf(fid,formatSpec,'hourofyear');
formatSpec = '%3s\n';
fprintf(fid,formatSpec,'mpf');
[mm,nn]=size(data);
for ii=1:1:mm
    for jj=1:1:nn
        if jj==1
            formatSpec = '%3d\t';
            fprintf(fid,formatSpec,data(ii,jj));
        else
            formatSpec = '%10.5f\n';
            fprintf(fid,formatSpec,data(ii,jj));
        end
    end
end
fclose(fid);

%% 数据说明
% data_file =
%
% I:\数据\2000-2011网格化产品\MODIS__MeltPondFraction__UHAM-CliSAP-ICDC__v02__12.5km__2000257.nc
%
% Source:
%            I:\数据\2000-2011网格化产品\MODIS__MeltPondFraction__UHAM-CliSAP-ICDC__v02__12.5km__2000257.nc
% Format:
%            classic
% Global Attributes:
%            Conventions               = 'CF-1.6'
%            Source                    = 'Satellite visible spectroradiometry'
%            institution               = 'Institute of Oceanography, University of Hamburg'
%            title                     = '8-day composite MODIS Arctic sea ice melt pond cover fraction [in percent] derived from 500 m resolution artifical neural network output'
%            time_coverage_start       = '2000'
%            time_coverage_end         = '2011'
%            product_version           = '2.0'
%            keywords                  = 'Earth Science > Cryosphere > Sea Ice > Melt Ponds\n Earth Science > Oceans > Sea Ice > Sea Ice Concentration\n Earth Science > Climate Indicators > Cryospheric Indicators > Sea Ice Concentration,Geographic Region > Northern Hemisphere,Vertical Location > Sea Surface, ICDC, University of Hamburg, Germany'
%            cdm_data_type             = 'grid'
%            date_created              = '2015-Jul-10'
%            creator_name              = 'UHAM-ICDC'
%            creator_url               = 'http://icdc.zmaw.de'
%            creator_email             = 'stefan.kern@uni-hamburg.de'
%            geographic_region         = 'Northern hemisphere'
%            geospatial_lat_min        = '60.0N'
%            geospatial_lat_max        = '90.0N'
%            geospatial_lon_min        = '0.0E'
%            geospatial_lon_max        = '360.0E'
%            geospatial_lon_resolution = '12.5km'
%            geospatial_lat_resolution = '12.5km'
%            license                   = 'Free and open access'
%            platform                  = 'EOS-TERRA'
%            sensor                    = 'MODIS'
%            Projection                = 'Polar-Stereographic with tangential plane at 70degN (NSIDC)'
%            Comment                   = 'This is version 2.0 of the melt pond fraction and open water fraction data set from University of Hamburg for periods May 9 to September 13 every year, main changes to version 1 are a) a bias reduction of the melt pond fraction of 8% and b) a bias reduction of the open water fraction by 3%; see M鋕ynen, M., S. Kern, A. R鰏el, and L. T. Pedersen, On the estimation of melt pond fraction on the Arctic sea ice with ENVISAT WSM images, Trans. Geosci. Rem. Sens., 52(11), 7366-7379, doi:10.1109/TGRS.2014.2311476, 2014, and Kern, S., M. Zygmuntowska, K. Khvorostovsky, G. Spreen, N. Ivanova, and A. Beitsch, ESA CCI Sea Ice ECV Project Report D4.1 Product Validation & Intercomparison Report (PVIR), SICCI-PVIR, v1.1,25-02-2015; Melt pond fraction and its standard deviation are set to 101 (the fill value) for grid cells with less than 63 or 563 good values (63 and 563 corresponds to 10% and 90% of the theoretical number of 500 m sub grid cells), Melt pond fraction and its standard deviation are set to 0 for ice concentrations less than 15%, Land mask from AMSRE_gsfc_n.hdf, 0 is open water, 1 is land; see Roesel, A., L. Kaleschke, G. Birnbaum, Melt ponds on Arctic sea ice determined from MODIS satellite data using an artificial neural network, The Cryosphere, 6, 431-446, doi:10.5194/tc-6-431-2012, 2012'
% Dimensions:
%            x    = 608
%            y    = 896
%            time = 1
% Variables:
%     x
%            Size:       608x1
%            Dimensions: x
%            Datatype:   single
%            Attributes:
%                        standard_name = 'projection_x_coordinate'
%                        long_name     = 'x coordinate of projection'
%                        units         = 'km'
%                        grid_spacing  = '12.5 km'
%     y
%            Size:       896x1
%            Dimensions: y
%            Datatype:   single
%            Attributes:
%                        standard_name = 'projection_y_coordinate'
%                        long_name     = 'y coordinate of projection'
%                        units         = 'km'
%                        grid_spacing  = '12.5 km'
%     lat
%            Size:       608x896
%            Dimensions: x,y
%            Datatype:   single
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude coordinate'
%                        units         = 'degrees_north'
%     lon
%            Size:       608x896
%            Dimensions: x,y
%            Datatype:   single
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude coordinate'
%                        units         = 'degrees_east'
%     time
%            Size:       1x1
%            Dimensions: time
%            Datatype:   double
%            Attributes:
%                        standard_name = 'time'
%                        units         = 'days since 2000-01-01'
%     mpf
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        long_name     = 'melt_pond_area_fraction_at_top_of_sea_ice'
%                        units         = 'percent'
%                        valid_range   = [0.00e+00 1.00e+02]
%                        missing_value = 103
%                        _FillValue    = 101
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     mpf_stdev
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        long_name     = 'melt_pond_area_fraction_at_top_of_sea_ice standard_deviation'
%                        cell_methods  = 'area: standard_deviation (interval: 500 m comment: sampled 25 x 25 grid cells)'
%                        units         = 'percent'
%                        valid_range   = [0.00e+00 1.00e+02]
%                        missing_value = 103
%                        _FillValue    = 101
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     mpf_noofdata
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        long_name     = 'number_of_valid_data_per_grid_cell'
%                        units         = '1'
%                        valid_range   = [0.00e+00 6.25e+02]
%                        missing_value = 1e+03
%                        _FillValue    = 1e+03
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     owf
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        standard_name = 'sea_area_fraction'
%                        long_name     = 'sea_area_fraction'
%                        units         = 'percent'
%                        valid_range   = [0.00e+00 1.00e+02]
%                        missing_value = 103
%                        _FillValue    = 101
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     mpf90
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        long_name     = 'melt_pond_area_fraction_at_top_of_sea_ice_assuming_clear_sky'
%                        units         = 'percent'
%                        valid_range   = [0.00e+00 1.00e+02]
%                        missing_value = 103
%                        _FillValue    = 101
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     mpf90_stdev
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        long_name     = 'melt_pond_area_fraction_at_top_of_sea_ice_assuming_clear_sky standard_deviation'
%                        cell_methods  = 'area: standard_deviation (interval: 500 m comment: sampled 25 x 25 grid cells)'
%                        units         = 'percent'
%                        valid_range   = [0.00e+00 1.00e+02]
%                        missing_value = 103
%                        _FillValue    = 101
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     owf90
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        long_name     = 'sea_area_fraction_assuming_clear_sky'
%                        units         = 'percent'
%                        valid_range   = [0.00e+00 1.00e+02]
%                        missing_value = 103
%                        _FillValue    = 101
%                        add_offset    = 0
%                        coordinates   = 'lat lon time'
%     land
%            Size:       608x896x1
%            Dimensions: x,y,time
%            Datatype:   int16
%            Attributes:
%                        standard_name = 'land_binary_mask'
%                        long_name     = 'land_mask'
%                        units         = '1'
%                        valid_range   = [0.00e+00 1.00e+00]
%                        coordinates   = 'lat lon time'
%%
clear all;
close all;
clc;

data_dir = '/home/share/modis_meltpond_areafraction_data/';
data_name_prefix = 'MODIS__MeltPondFraction__UHAM-CliSAP-ICDC__v02__12.5km__';
data_year=2005;
data_day=129;
char_data_year=num2str(data_year);
char_data_day=num2str(data_day);
data_file = [data_dir, data_name_prefix, char_data_year, char_data_day, '.nc'];
%% 确定位置
%  220，70N 点
lon_goal = 220;    lat_goal = 85 ;
char_lon=num2str(lon_goal);char_lat=num2str(lat_goal);

LON =ncread(data_file,'lon');
LAT =ncread(data_file,'lat');
LONLON=LON(:);
LATLAT=LAT(:);
[sLON,index_lon]= sort(abs(LONLON - lon_goal) );
[sLAT,index_lat]= sort(abs(LATLAT - lat_goal) );

% enhance roubusty and compatible
% use more variables to replace consts
for n=1:1000
    
    
    
    % n=1000;
    [RESULT,ia,ib] = intersect(index_lat(1:n), index_lon(1:n));
    if ~isempty(RESULT)
        
        break;
    end
    disp(['n=',num2str(n)])
end

lon_find =LON (RESULT(1));
lat_find =LAT (RESULT(1));

index_lon_find=index_lon(ib);
index_lat_find=index_lat(ia);
[i,j] =ind2sub([608,896],RESULT(1)); % i,j 是这一步最终需要的产品

% 检查
LON (i,j )
LAT (i,j )

%% 提取

mpf_thisyear = zeros( 1, 16 );
mpf_thisyear(:) = nan;

k=0;
for number = 129: 8: 249
    k=k+1;
    data_file = [data_dir, data_name_prefix, char_data_year, num2str( number ), '.nc'];
    mpf = ncread( data_file, 'mpf');
    mpf_thisyear ( k ) = mpf( i, j );
%     ncdisp( data_file );
end
mpf_thisyear(mpf_thisyear == 103) = nan;




%% 保存
day=129:8:249;
day=day';
mpf_thisyear = mpf_thisyear';
mpf_thisyear = mpf_thisyear*0.01; % 转换为小数制
mpf_thisyear( isnan(mpf_thisyear) ) = 103;
data =cat (2,day,mpf_thisyear);


filename=['/home/luy/icepack/cesm_meltpond_pe/run_real/','modis_',char_data_year,'_',char_lon,'_',char_lat,'_08day.txt'];

fid =fopen(filename,'wt');
formatSpec = '%9s\t';
fprintf(fid,formatSpec,'dayofyear');
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

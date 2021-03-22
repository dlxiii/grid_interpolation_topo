% read netcdf topographic data to make a rectangular grid
% Yulong Wang
% 2021/03/13
clear all;
clearvars; clc;

% Which system am I using?
if ismac    % On Mac
    basedir = '/Users/yulong/GitHub/';
    basedir2 = '/Volumes/Yulong/';
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
elseif isunix       % Unix?
    basedir = '/home/usr0/n70110d/';
    addpath([basedir,'github/fvcomtoolbox/']);
    addpath([basedir,'github/fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'github/fvcomtoolbox/utilities/']);
    addpath([basedir,'github/fvcomtoolbox/custom/']);
elseif ispc     % Or Windows?
    basedir = 'C:/Users/Yulong WANG/Documents/GitHub/';      
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
end

%% IO setting

ncfile.name = '/Volumes/GoogleDrive/Shared drives/Wang Yulong 2017.9/data/bathymetry/meiji_elevation.nc';
otfile =      '/Volumes/GoogleDrive/Shared drives/Wang Yulong 2017.9/data/bathymetry/meiji_elevation_50grids.csv';

% ncfile.name = '/Volumes/Yulong/GitHub/database_tokyobay/FVCOM_case/depth_0090-05+06+07.nc';
% otfile = '/Volumes/GoogleDrive/Shared drives/Wang Yulong 2017.9/data/bathymetry/reiwa_elevation_50grids.csv';

% ncfile.name = './inp/meiji_elevation.nc';
% otfile = './otp/meiji_elevation_500grids.csv';

%% READ

% ncdisp(ncfile);
% from nc information, we list the following variables
%=================================================================================
% variable   long name                         units             dimensions
%=================================================================================
% lon:	     longitude                         degrees_east     [lon]
% lat:	     latitude                          degrees_north    [lat]
% Band1:     GDAL Band Number 1                meters           [lon,lat]
%=================================================================================

vartoread = {'lon','lat'};
vartoread2d = {'Band1'};

for i = 1:length(vartoread)
    ncfile.(vartoread{i}) = ncread(ncfile.name,vartoread{i});
    ncfile.(vartoread{i}) = double(ncfile.(vartoread{i}));
    fprintf(['Read 0d var ',vartoread{i},'.\n']);
end

if exist('vartoread2d','var') == 1
    % 2d variables
    for i = 1:length(vartoread2d)
        ncfile.(vartoread2d{i}) = ncread(ncfile.name,vartoread2d{i},[1, 1],[Inf, Inf]);
        fprintf(['Read 2d var ',vartoread2d{i},'.\n']);
    end
end

% basic analysis
ncfile.lon_res = (max(ncfile.lon) - min(ncfile.lon))/(length(ncfile.lon)-1);
ncfile.lat_res = (max(ncfile.lat) - min(ncfile.lat))/(length(ncfile.lat)-1);
% cvt.R = 6378.137;
% cvt.a = sin((ncfile.lon_res/2))^2 + cos(max(ncfile.lat)*pi/180)*cos(min(ncfile.lat)*pi/180)*(sin(ncfile.lon_res/2))^2;
% cvt.c = 2 * atan2(sqrt(cvt.a), sqrt(1-cvt.a));
% cvt.d = cvt.R * cvt.c *1000;
% if assume that 1 deg = 111km, or 111000m.
fprintf(['Resolution of data is ',num2str(ncfile.lon_res*111000),'m X ',num2str(ncfile.lat_res*111000),'m.\n']);

%% EXTRACT DATA

% before interpolation, extract the interesting area.
% [grid_lon, grid_lat] = meshgrid(ncfile.lon, ncfile.lat);
% idx = grid_lon < 140.40 & grid_lon > 139.13 & grid_lat < 35.71 & grid_lat > 34.86;
%[3577,7060],[5937,9360]
lon_lat = [139.55,140.15,34.86,35.71];
lon_lat_list = [...
    lon_lat(1),lon_lat(3);...
    lon_lat(1),lon_lat(4);...
    lon_lat(2),lon_lat(3);...
    lon_lat(2),lon_lat(4);...
    ];
lon = ncfile.lon(...
    ncfile.lon > lon_lat(1) &...
    ncfile.lon < lon_lat(2));
lat = ncfile.lat(...
    ncfile.lat > lon_lat(3) &...
    ncfile.lat < lon_lat(4));
[grid_lon, grid_lat] = meshgrid(lon, lat);
grid_lon = reshape(grid_lon,[length(lon)*length(lat),1]);
grid_lat = reshape(grid_lat,[length(lon)*length(lat),1]);
[row,~] = find(ncfile.lon > lon_lat(1) & ncfile.lon < lon_lat(2));
[~,col] = find(ncfile.lat' > lon_lat(3) & ncfile.lat' < lon_lat(4));
ele = rot90(ncfile.Band1(min(row):max(row),min(col):max(col)));
ele_new = ele;
ele_new(isnan(ele_new)) = 1;
grid_ele = reshape(ele_new,[length(lon)*length(lat),1]);
% imagesc(ele_new)

%% INTERPOLATION
% interpolation
data.interp = scatteredInterpolant(...
    grid_lon,...
    grid_lat,...
    grid_ele,...
    'natural');
use_lonlat = 0;
if use_lonlat
    % use lon lat coordinate
    data.lonint = (lon_lat(1):50/111000:lon_lat(2))';
    data.latint = (lon_lat(3):50/111000:lon_lat(4))';
    data.valint = data.interp({data.lonint,data.latint})';
else
    % use x y coordinate
    % Generate coordinates
    inputConf.utmZone = {'54 N'};
    UTMzone = regexpi(inputConf.utmZone,'\ ','split');
    for s = 1:length(lon_lat_list)
        [lon_lat_list(s,3),lon_lat_list(s,4),~,~] = wgs2utm(...
            lon_lat_list(s,2),lon_lat_list(s,1),...
            str2double(char(UTMzone{1}(1))),char(UTMzone{1}(2)));
    end
    clear s UTMzone
    
    mesh = 50;
    data.xint = (max(lon_lat_list(1:2,3)):mesh:min(lon_lat_list(3:4,3)))';
    data.yint = (max(lon_lat_list(1,4),lon_lat_list(3,4)):mesh:min(lon_lat_list(2,4),lon_lat_list(4,4)))';
    [data.grid_x, data.grid_y] = meshgrid(data.xint, data.yint);
    
    for x = 1:length(data.xint)
        for y = 1:length(data.yint)
            [data.grid_lat(y,x),data.grid_lon(y,x)] = utm2deg(...
                data.grid_x(y,x),data.grid_y(y,x),inputConf.utmZone{1});
            data.valint(y,x) = data.interp({data.grid_lon(y,x),data.grid_lat(y,x)})';
        end
    end
    clear s UTMzone
end

%% EXPORT
fileID = fopen(otfile,'w');
fprintf(fileID,'%12s %12s %12s\n','lon,','lat,','elevation,');
if use_lonlat
    for i = 1:length(data.lonint)
        for j = 1:length(data.latint)
            fprintf(fileID,'%11.6f%1s %11.6f%1s %11.2f%1s \n',...
                data.lonint(i),',',data.latint(j),',',data.valint(length(data.latint)+1-j,i),',');
        end
    end
else
    for i = 1:length(data.xint)
        for j = 1:length(data.yint)
            fprintf(fileID,'%11.6f%1s %11.6f%1s %11.2f%1s \n',...
                data.grid_lon(j,i),',',data.grid_lat(j,i),',',data.valint(length(data.yint)+1-j,i),',');
        end
    end
end
fclose(fileID);
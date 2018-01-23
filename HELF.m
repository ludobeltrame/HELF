function infectiveMeta = HELF (time,temp,tempMin,tempMax,lat,ti_values,rain,params_hydro,params_fluke)
% this function runs the Hydro-Epidemiological model for Liver Fluke (HELF)
%
% INPUTS:
% - time: array of size (# days in the sim period,1) containing serial date 
%         numbers for the simulation period
% - temp, tempMin, tempMax: matrixes of size (# days in the sim period,
%                           # grid cells overlapping with catchment area) 
%                           containing daily mean, min and max temperature  
%                           [°C] for the simulation period)
% - lat: array of size (# grid cells overlapping with catchment area,1)
%        containing their latitude [°]
% - ti_values: name of tif file containing the matrix (with # of elements 
%              equal to the # of grid cells in the domain) of Topographic
%              Index (TI) values
% - rain: array of size (# days in the sim period,1) containing catchment 
%         average daily rainfall [mm] for the simulation period)
% - params_hydro: array of size (1,7) containing parameters for the hydro 
%                 model component of HELF
% - params_fluke: array of size (1,22) containing parameters for the fluke 
%                 model component of HELF
% OUTPUT:
% - infective_meta: matrix of size (# days in the sim period,# TI classes)
%                   containing the abundance of infective metacercariae on
%                   pasture for the sim period

%% calculate potential evapotranspiration

pet_gridded = calculateDailyPetHargreaves (time, lat, temp, tempMin, tempMax);
% (matrix of size (# days in the sim period,# grid cells
% overlapping with catchment area) containing daily pet [mm] over the  
% sim period)

pet_catchment_average = mean(pet_gridded,2);
% (array of size (# days in the sim period,1) containing catchment average 
% daily pet for the sim period)

%% generate Topographic Index (TI) classes

% read TI map
[topindexAtb, topindexR] = geotiffread(ti_values);

% choose number of TI classes
topIndexCount = 25; 

% generate TI classes
[topindexResult, topindexBinLookup] = makeTIclasses(topindexAtb, topIndexCount);

TI = topindexResult(:,1);   % TI value for each class
freq = topindexResult(:,2); % portion of catchment associated with each class

%% run hydrological component of HELF 

sTimeSeries = HELF_hydro_model_component(rain,pet_catchment_average,TI,freq,params_hydro);
% (matrix of size (# days in the sim period,# TI classes) containing daily
% sat def map for the sim period)

%% run liver fluke component of HELF

% recover parameter labels
paramOrder = {'S1MinTemp' 'S1PeakTemp' 'S1PeakRate' 'S1MortRate' 'S1MaxTemp' 'S1TrigSatDef' 'S1TrigMinTemp' 'S2SatDef' 'S2SatSlope' 'S2MinTemp'	'S2PeakTemp' 'S2MaxTemp' 'S3SatDef' 'S3SatSlope' 'S3MinATemp' 'S3PeakRate' 'S3MortRate'	'S4TempMin' 'S4TempPeak' 'S4TempMax' 'S4PeakRate' 'S4MinRate'};

% choose whether to keep track of the development progress within the
% life-cycle stages (i.e. cohortMode=1 vs. cohortMode=0)
cohortMode = 1;

% calculate average temperature over the catchment
temp_catchment_average = mean(temp,2);
% (array of size (# days in the sim period,1) containing catchment average 
% daily temperature for the sim period

% obtain abundance of developed eggs, infected snails, developed snail
% infections and infective metacercariae on pasture for each day and TI
% class
[developedEggs, infectedSnails, developedSnailInf, infectiveMeta] = HELF_fluke_model_component (params_fluke, paramOrder, cohortMode, time, sTimeSeries, temp_catchment_average);


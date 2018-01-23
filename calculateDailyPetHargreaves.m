function [pet] = calculateDailyPetHargreaves (timeDimention, t, tmin, tmax, latitude)
% this function calculates daily potential evapotranspiration using 
% Hargreaves equation (*)
% INPUTS:
% - timeDimention: array of size (# days in the sim period,1) containing 
%                  serial date numbers for the sim period
% - t, tmin, tmax: matrixes of size (# days in the sim period,# grid cells
%                  overlapping with catchment area) containing daily mean, 
%                  min and max temperature [°C] for the sim period
% - latitude: array of size (# grid cells overlapping with catchment area,
%             1) containing their latitude [°]
% OUTPUT:
% - pet: matrix of size (# days in the sim period,# grid cells overlapping
%        with catchment area) containing daily potential evapotranspiration
%        [mm] for the sim period

% (*) references: 
% - Allen et al. (1998) Crop evapotranspiration - Guidelines for computing 
%   crop water requirements. FAO irrigation and drainage paper 56.
% - Droogers and Allen (2002) Estimating reference evapotranspiration under 
%   inaccurate data conditions. Irrigations and Drainage Systems 16: 33-45.

%% transform lat degrees in lat rad

lat = degtorad(latitude);

%% define julian days for sim period (array of daily values)

julianTimeDim = [(1:365)'; (1:366)'; (1:365)'; (1:365)'; (1:365)'; (1:366)'; (1:365)'; (1:365)'; (1:365)'; (1:366)'; (1:365)'; (1:365)'];
% *** this is for 1999-2010, has to be changed for a different period ***

%% calculate relative distance between earth and sun (array of daily values) [-]

dr = 1 + 0.033 * cos ((2*pi*julianTimeDim)/365); 

%% calculate solar declination (array of daily values) [rad]

delta = 0.409 * sin (((2*pi*julianTimeDim)/365)-1.39); 

%% calculate sunset hour angle (matrix # of days, # of grid cells) [rad]

ws = nan(length(julianTimeDim),length(lat));

for i=1:length(lat)
    ws(:,i) = acos((-tan(lat(i)))*(tan(delta)));
end

%% calculate incoming solar energy (matrix # of days, # of grid cells) [mm/day]

Gsc = 0.0820; % solar constant [MJ / m2 min]

Ra=nan(size(ws));

for i=1:length(lat)
    Ra(:,i) = ((24*60)/pi * Gsc * dr).*(((ws(:,i)*sin(lat(i))).*sin(delta))+((cos(lat(i))*cos(delta)).*sin(ws(:,i))));
end

%% calculate pet

TD = abs(tmax - tmin);

pet = nan(size(Ra)); % mm/d

for i = 1:size(pet,1)
    for j = 1:size(pet,2)
        pet(i,j) = 0.0023 * Ra(i,j) * sqrt(TD(i,j)) * 0.4082 * (t(i)+17.8);
    end
end

end


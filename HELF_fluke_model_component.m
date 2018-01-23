function [S1, S2, S3, S4] = HELF_fluke_model_component (params, paramOrder, cohortMode, time, s_time_series, temp)
% this function simulates the liver fluke life-cycle driven by temperature 
% and soil moisture conditions
% INPUTS:
% - params          - array of size (1,22) containing parameter values 
% - paramOrder      - array of size (1,22) containing parameter names 
% - cohortMode      - scalar which can take value either 1 for cohort
%                     modelling or 0 otherwise
% - time            - array of size (# days in the sim period,1) containing serial date numbers for the sim period
% - s_time_series   - matrix of size (# days in the sim period,# TI
%                     classes) containing soil saturation maps for the sim period
% - temp            - array of size (# days in the sim period,1) containing 
%                     catchment average daily mean temperature [°C] for the 
%                     sim period
% OUTPUTS:
% - S1: matrix of size (# days in the sim period,# TI classes) containing 
%       the abundance of eggs developed on pasture                     
% - S2: matrix of size (# days in the sim period,# TI classes) containing
%       the abundance of infected snails
% - S3: matrix of size (# days in the sim period,# TI classes) containing 
%       the abundance of developed snail infections              
% - S4: matrix of size (# days in the sim period,# TI classes) containing 
%       the abundance of infective metacercariae on pasture
                  
%% set an egg scenario

egg = zeros(size(time));

% uncomment for option 1: introduction of a constant number of eggs ON THE 
% FIRST DAY OF THE SIMULATION PERIOD ONLY
% egg(1) = 100;

% uncomment for option 2: introduction of a constant number of eggs EVERY 
% DAY IN THE SIMULATION PERIOD 
egg = ones(size(egg)) * 100;

% uncomment for option 3: introduction of a constant number of eggs EVERY 
% DAY IN THE SIMULATION PERIOD, DURING GRAZING MONTHS ONLY (e.g. MAY TO NOV)
% for dayIndex=1:size(time,1)
%     v = datevec(time(dayIndex));
%     m = v(2);
%     if(m>=5 && m <12)
%         egg(dayIndex) = 100;
%     end
% end

%% recover parameters and inputs to run the model

sTimeSeries = s_time_series;

count = 0;

for index = 1:size(sTimeSeries,2) 
    count = count + 1;
    satu = sTimeSeries(:,index);  
    FlukeEggToMetaProcessor(count) = FlukeEggToMeta.Create(cohortMode, params, paramOrder); % set cohort mode and model parameters based on paramOrder defined by the user
    FlukeEggToMetaProcessor(count).InitialiseTimeSeries(time, temp, satu); % set model inputs (temperature and sat def)
end

for index = 1:(size(FlukeEggToMetaProcessor,2)) 
    FlukeEggToMetaProcessor(index).Run(egg); % run the model (given the chosen egg scenario)
end

%% extrapolate model outputs

% initialise variables to store model outputs in
% (matrixes of size (# days in the sim period,# TI classes)
S1 = zeros(size(time,1),size(FlukeEggToMetaProcessor,2));
S2 = zeros(size(time,1),size(FlukeEggToMetaProcessor,2));
S3 = zeros(size(time,1),size(FlukeEggToMetaProcessor,2));
S4 = zeros(size(time,1),size(FlukeEggToMetaProcessor,2));

for index = 1:(size(FlukeEggToMetaProcessor,2))
    S1(:,index) = FlukeEggToMetaProcessor(index).S1SumEgg;
    S2(:,index) = FlukeEggToMetaProcessor(index).S2SumLarva;
    S3(:,index) = FlukeEggToMetaProcessor(index).S3SumInfectedSnail;
    S4(:,index) = FlukeEggToMetaProcessor(index).S4SumMetacercariae;
end

end


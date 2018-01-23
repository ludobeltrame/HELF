classdef FlukeEggToMeta < handle
    
    properties
        
        S1CohortEgg             % S1Eggs
        S1AgeRate               % S1Develop
        S1Mortality             % S1Mort
        
        S2CohortLarva           % S2Larvae
        S2Suitability           % S2Suitability
        
        S3CohortInfectedSnail   % S3SnailInfections
        S3AgeRate               % S3Develop
        S3Mortality             % S3Mort
        
        S4CohortMetacercariae   % S4Metacercariae
        S4AgeRate               % S4Survival
        
        % outputs 
        S1SumEgg
        S2SumLarva
        S3SumInfectedSnail
        S4SumMetacercariae
        
        % inputs
        timeDimention
        temperatureTimeSeries
        saturationTimeSeries
        
    end
    
    properties(SetAccess = private, GetAccess = private)
        
        S1TriggerMinSaturationDeficit
        S1TriggerMinTemp              
        
        S2SuitabilitySatu
        S2SuitabilityTemp
        
        S3AgeRateSatu
        S3AgeRateTemp
        
    end
  
%%
    methods(Static)
                
        function obj = Create(cohortMode, params, paramOrder)
            
            % initialise all parameters we could possibly have to zero
            
            S1AgeRateMinTemp = 0;
            S1AgeRatePeakTemp= 0;
            S1AgeRateMaxTemp= 0; 
            S1AgeRatePeak= 0;
            S1MortalityRate= 0;
            S1AgeRateMortalityMultiple= 0;
            
            S1TriggerMinSaturationDeficit= 0;
            S1TriggerMinTemp= 0;
            
            S2MinSaturationDeficit= 0;
            S2PeakSaturationDeficit= 0;
            S2SuccessMinTemp= 0;
            S2SuccessPeakTemp= 0;
            S2SuccessMaxTemp= 0;
            
            S3SnailMinSaturationDeficit= 0;
            S3SnailPeakSaturationDeficit= 0;
            S3SnailAgeRateMinATemp= 0;
            S3SnailAgeRateMinBTemp= 0;
            S3SnailAgeRateMaxATemp= 0;
            S3SnailAgeRateMaxBTemp= 0;
            S3SnailAgeRatePeak= 0;
            S3SnailMortalityRate= 0;
            
            S4SurvivabilityTempMin= 0;
            S4SurvivabilityTempPeak= 0;
            S4SurvivabilityTempMax= 0;
            S4SurvivabilityPeakRate= 0;
            S4SurvivabilityMinRate=0;
            
            % recover model parameters based on user-specified paramOrder
            
            for paramIndex=1:length(paramOrder)
                switch(char(paramOrder(paramIndex))) 
                    case 'S1MinTemp'
                        S1AgeRateMinTemp = params(paramIndex);
                    case 'S1PeakTemp'
                        S1AgeRatePeakTemp = params(paramIndex); 
                    case 'S1MaxTemp' 
                        % if S1MaxTemp is set by the user, the eggs development  
                        % rate function is triangular.
                        % if S1MaxTemp is NOT set by the user, it remains equal 
                        % to zero, and the development rate function is linear.
                        S1AgeRateMaxTemp = params(paramIndex);
                    case 'S1PeakRate'
                        S1AgeRatePeak = params(paramIndex);
                    case 'S1MortRate'
                        S1MortalityRate = params(paramIndex);
                    case 'S1MortMultiple'
                        S1AgeRateMortalityMultiple = params(paramIndex);
                    case 'S1TrigSatDef'
                        S1TriggerMinSaturationDeficit = params(paramIndex);
                    case 'S1TrigMinTemp'
                        S1TriggerMinTemp = params(paramIndex);
                        
                    case 'S2SatDef'
                        S2MinSaturationDeficit = params(paramIndex);
                    case 'S2SatSlope'
                        S2PeakSaturationDeficit = S2MinSaturationDeficit * params(paramIndex);
                        % (watch out: this is based on the assumption that
                        % 'S2SatDef' was set first)
                    case 'S2MinTemp'
                        S2SuccessMinTemp = params(paramIndex);
                    case 'S2PeakTemp'
                        S2SuccessPeakTemp = params(paramIndex);
                    case 'S2MaxTemp'
                        S2SuccessMaxTemp = params(paramIndex);
                        
                    case 'S3SatDef'
                        S3SnailMinSaturationDeficit = params(paramIndex);
                    case 'S3SatSlope'
                        S3SnailPeakSaturationDeficit =S3SnailMinSaturationDeficit * params(paramIndex);
                        % (watch out: this is based on the assumption that
                        % 'S3SatDef' was set first)                        
                    case 'S3MinATemp'
                        % option1: hp of constant snail infection
                        % development above a min temp threshold.
                        % (watch out: S3MinATemp' has to be set before
                        % 'S3PeakRateTemp' and 'S3MaxTemp')
                        S3SnailAgeRateMinATemp = params(paramIndex);
                        S3SnailAgeRateMinBTemp = params(paramIndex);
                        S3SnailAgeRateMaxATemp = 100000;
                        S3SnailAgeRateMaxBTemp = 100000;
                    case 'S3PeakRateTemp'
                        % option2: hp of temperature dependent development 
                        % rate for snail infection.
                        % (minB has to be set equal to maxA to make a triangular function)  
                        S3SnailAgeRateMinBTemp = params(paramIndex);
                        S3SnailAgeRateMaxATemp = params(paramIndex);
                    case 'S3MaxTemp'
                        S3SnailAgeRateMaxBTemp = params(paramIndex);
                    case 'S3PeakRate'
                        S3SnailAgeRatePeak = params(paramIndex);
                    case 'S3MortRate'
                        S3SnailMortalityRate = params(paramIndex);
                        
                    case 'S4TempMin'
                        S4SurvivabilityTempMin = params(paramIndex);
                    case 'S4TempPeak'
                        S4SurvivabilityTempPeak = params(paramIndex);
                    case 'S4TempMax'
                        S4SurvivabilityTempMax = params(paramIndex);
                    case 'S4PeakRate'
                        S4SurvivabilityPeakRate = params(paramIndex);
                    case 'S4MinRate'
                        S4SurvivabilityMinRate = params(paramIndex);
                end
            end
            
            % create object with the defined parameters
            
            obj = FlukeEggToMeta(cohortMode,...
                S1AgeRateMinTemp, ...
                S1AgeRatePeakTemp, ...
                S1AgeRateMaxTemp, ... 
                S1AgeRatePeak, ...
                S1MortalityRate, ...
                S1AgeRateMortalityMultiple, ...
                ...
                S1TriggerMinSaturationDeficit, ...
                S1TriggerMinTemp, ...
                ...
                S2MinSaturationDeficit, ...
                S2PeakSaturationDeficit, ...
                S2SuccessMinTemp, ...
                S2SuccessPeakTemp, ...
                S2SuccessMaxTemp, ...
                ...
                S3SnailMinSaturationDeficit, ...
                S3SnailPeakSaturationDeficit, ...
                S3SnailAgeRateMinATemp, ...
                S3SnailAgeRateMinBTemp, ...
                S3SnailAgeRateMaxATemp, ...
                S3SnailAgeRateMaxBTemp, ...
                S3SnailAgeRatePeak, ...
                S3SnailMortalityRate, ...
                ...
                S4SurvivabilityTempMin, ...
                S4SurvivabilityTempPeak, ...
                S4SurvivabilityTempMax, ...
                S4SurvivabilityPeakRate, ...
                S4SurvivabilityMinRate);
            
        end
        
    end
    
%%
    methods
                
        function obj = FlukeEggToMeta(cohortMode,...
                S1AgeRateMinTemp, ...
                S1AgeRatePeakTemp, ...
                S1AgeRateMaxTemp, ... 
                S1AgeRatePeak, ...
                S1MortalityRate, ...
                S1AgeRateMortalityMultiple, ...
                ...
                S1TriggerMinSaturationDeficit, ...
                S1TriggerMinTemp, ...
                ...
                S2MinSaturationDeficit, ...
                S2PeakSaturationDeficit, ... 
                S2SuccessMinTemp, ...
                S2SuccessPeakTemp, ...
                S2SuccessMaxTemp, ...
                ...
                S3SnailMinSaturationDeficit, ...
                S3SnailPeakSaturationDeficit, ...
                S3SnailAgeRateMinATemp, ...
                S3SnailAgeRateMinBTemp, ...
                S3SnailAgeRateMaxATemp, ...
                S3SnailAgeRateMaxBTemp, ...
                S3SnailAgeRatePeak, ...
                S3SnailMortalityRate, ...
                ...
                S4SurvivabilityTempMin, ...
                S4SurvivabilityTempPeak, ...
                S4SurvivabilityTempMax, ...
                S4SurvivabilityPeakRate, ...
                S4SurvivabilityMinRate ...
                )
            
            obj.S1TriggerMinSaturationDeficit = S1TriggerMinSaturationDeficit;
            obj.S1TriggerMinTemp = S1TriggerMinTemp;
            
            % checks on model params
            
            if ~(S1AgeRateMinTemp < S1AgeRatePeakTemp)
                error('S1 temp parameters invalid');
            end
            
            if ~(S1AgeRatePeak > 0)
                error('S1AgeRatePeak invalid');
            end
            
            if ~(S2MinSaturationDeficit > S2PeakSaturationDeficit) 
                error('S2MinSaturationDeficit must be > S2PeakSaturationDeficit');
            end
            
            if ~(S3SnailMinSaturationDeficit > S3SnailPeakSaturationDeficit) 
                error('S3SnailMinSaturationDeficit must be > S3SnailPeakSaturationDeficit');
            end
            
            if ~(S2SuccessMinTemp < S2SuccessPeakTemp && ...
                    S2SuccessPeakTemp < S2SuccessMaxTemp)
                error('S2 temp parameters invalid');
            end
            
            if ~(S3SnailAgeRateMinATemp <= S3SnailAgeRateMinBTemp && ...
                    S3SnailAgeRateMinBTemp <= S3SnailAgeRateMaxATemp && ...
                    S3SnailAgeRateMaxATemp <= S3SnailAgeRateMaxBTemp)
                error('Snail temp parameters invalid');
            end
            
            if ~(S4SurvivabilityTempMin < S4SurvivabilityTempPeak && ...
                    S4SurvivabilityTempPeak < S4SurvivabilityTempMax)
                error('S4Survivability temp parameters invalid');
            end
            
            if ~(S4SurvivabilityPeakRate > 0)
                error('S4SurvivabilityPeakRate parameter invalid');
            end
            
            if ~(S4SurvivabilityMinRate > 0)
                error('S4SurvivabilityMinRate parameter invalid');
            end
            
            if ~(S4SurvivabilityPeakRate < S4SurvivabilityMinRate)
                error('S4Survivability rate parameters invalid');
            end
            
            if ~(S3SnailAgeRatePeak > 0)
                error('S3SnailAgeRatePeak parameter invalid');
            end
            
            % calculation of stage-specific development/mortality rates
            
            % EGGS
            
            % if cohortMode=1, use cohort-based version of the model. 
            % otherwise, use pooled version.
            if(cohortMode)
                obj.S1CohortEgg = stage_cycle_cohort();
            else
                obj.S1CohortEgg = stage_cycle_pooled();
            end
            
            obj.S1CohortEgg.stage_name = 'Eggs';
            obj.S1CohortEgg.InitialiseDevelopmentDistribution();
            
            % development rate
         
            if(S1AgeRateMaxTemp == 0) % (option of linear development rate)
                rateM = (S1AgeRatePeak - 0) / (S1AgeRatePeakTemp - S1AgeRateMinTemp); % slope
                rateB = -(rateM * S1AgeRateMinTemp);                                  % intercept
                obj.S1AgeRate = rateCalcLinear(rateB, rateM);
            else                      % (option of triangular development rate)
                obj.S1AgeRate = rateLinearFuzzyMembership(S1AgeRateMinTemp, S1AgeRatePeakTemp, S1AgeRatePeakTemp, S1AgeRateMaxTemp);
                obj.S1AgeRate.rateMember = S1AgeRatePeak;
                obj.S1AgeRate.rateNonMember = 0;
            end
            
            % mortality rate
            
            obj.S1Mortality = rateCalcConstant(S1MortalityRate);
            obj.S1CohortEgg.developmentMortality = S1AgeRateMortalityMultiple;
            
            % MIRACIDIA
            
            % here development rate is equal to the probability of finding 
            % a snail intermediate host (function of both sat def and temp), 
            % and mortality rate is set to 1-probability of finding one.
            
            obj.S2SuitabilitySatu = rateLinearFuzzyMembership(-100000,-100000,S2PeakSaturationDeficit,S2MinSaturationDeficit); % hp of no upper limit to saturation for snail presence
            obj.S2SuitabilityTemp = rateLinearFuzzyMembership(S2SuccessMinTemp,S2SuccessPeakTemp,S2SuccessPeakTemp,S2SuccessMaxTemp);
            obj.S2Suitability = rateCombine(obj.S2SuitabilitySatu, obj.S2SuitabilityTemp, 1);
            
            % if cohortMode=1, use cohort-based version of the model. 
            % otherwise, use pooled version.
            if(cohortMode)
                obj.S2CohortLarva = stage_cycle_cohort();
            else
                obj.S2CohortLarva = stage_cycle_pooled();
            end
            
            obj.S2CohortLarva.stage_name = 'Larva';
            obj.S2CohortLarva.development_rate = 1;
            
            % SNAIL INFECTIONS
            
            % if cohortMode=1, use cohort-based version of the model. 
            % otherwise, use pooled version.
            if(cohortMode)
                obj.S3CohortInfectedSnail = stage_cycle_cohort();
            else
                obj.S3CohortInfectedSnail = stage_cycle_pooled();
            end
            
            obj.S3CohortInfectedSnail.stage_name = 'Infected Snail';
            
            % development rate

            obj.S3AgeRateSatu = rateLinearFuzzyMembership(-100000,-100000,S3SnailPeakSaturationDeficit,S3SnailMinSaturationDeficit); % hp of no upper limit to saturation for snail presence and snail infection development 
            obj.S3AgeRateTemp = rateLinearFuzzyMembership(S3SnailAgeRateMinATemp,S3SnailAgeRateMinBTemp,S3SnailAgeRateMaxATemp,S3SnailAgeRateMaxBTemp);
            obj.S3AgeRate = rateCombine(obj.S3AgeRateSatu, obj.S3AgeRateTemp, S3SnailAgeRatePeak);
            
            % mortality rate
            
            obj.S3Mortality = rateCalcConstant(S3SnailMortalityRate);
            
            % METACERCARIAE
            
            % if cohortMode=1, use cohort-based version of the model. 
            % otherwise, use pooled version.
            if(cohortMode)
                obj.S4CohortMetacercariae = stage_cycle_cohort();
            else
                obj.S4CohortMetacercariae = stage_cycle_pooled();
            end
            
            % survival rate
            
            obj.S4CohortMetacercariae.mortality_rate = 0;
            obj.S4AgeRate = rateLinearFuzzyMembership(S4SurvivabilityTempMin,S4SurvivabilityTempPeak,S4SurvivabilityTempPeak,S4SurvivabilityTempMax);
            obj.S4AgeRate.rateMember = S4SurvivabilityPeakRate;
            obj.S4AgeRate.rateNonMember = S4SurvivabilityMinRate;
            
        end
        
        %% recover model inputs
        
        function InitialiseTimeSeries(obj, timeDimention, temperatureTimeSeries, saturationTimeSeries)
            
            obj.timeDimention = timeDimention;
            obj.temperatureTimeSeries = temperatureTimeSeries;
            obj.saturationTimeSeries = saturationTimeSeries;
            
            obj.S1AgeRate.timeSeries = temperatureTimeSeries;
            
            obj.S2SuitabilitySatu.timeSeries = saturationTimeSeries;
            obj.S2SuitabilityTemp.timeSeries = temperatureTimeSeries;
            
            obj.S3AgeRateSatu.timeSeries = saturationTimeSeries;
            obj.S3AgeRateTemp.timeSeries = temperatureTimeSeries;
            
            obj.S4AgeRate.timeSeries = temperatureTimeSeries;
            
        end
        
        %% calculate abundance of individuals maturing to next life-cycle stage
        
        function [ output_args ] = Run(obj, eggTimeSeries )
            
            % initialisations
            
            obj.S1CohortEgg.allocate(length(eggTimeSeries));
            obj.S2CohortLarva.allocate(length(eggTimeSeries));
            obj.S3CohortInfectedSnail.allocate(length(eggTimeSeries));
            obj.S4CohortMetacercariae.allocate(length(eggTimeSeries));
            
            days = size(eggTimeSeries,1);
            
            obj.S1SumEgg = zeros(days ,1);
            obj.S2SumLarva = zeros(days ,1);
            obj.S3SumInfectedSnail = zeros(days ,1);
            obj.S4SumMetacercariae = zeros(days ,1);
            
            % calculations
            
            for day = 1:days
                                
                obj.S1SumEgg(day) = sum(obj.S1CohortEgg.GetCount());
                obj.S2SumLarva(day) = sum(obj.S2CohortLarva.GetCount());
                obj.S3SumInfectedSnail(day) = sum(obj.S3CohortInfectedSnail.GetCount());
                obj.S4SumMetacercariae(day) = sum(obj.S4CohortMetacercariae.GetCount());
                
                % update development/mortality rates
                
                % EGGS
                
                obj.S1CohortEgg.development_rate = obj.S1AgeRate.getRate(day);
                obj.S1CohortEgg.mortality_rate   = obj.S1Mortality.getRate(day);
                
                % apply trigger for eggs hatching into miracidia:
                if(obj.temperatureTimeSeries(day) > obj.S1TriggerMinTemp &&...
                        obj.saturationTimeSeries(day) < obj.S1TriggerMinSaturationDeficit)
                    obj.S1CohortEgg.progressionTrigger = 1;
                else
                    obj.S1CohortEgg.progressionTrigger = 0;
                end
                
                % MIRACIDIA
                
                obj.S2CohortLarva.mortality_rate = 1 - obj.S2Suitability.getRate(day);
                
                % SNAIL INFECTIONS
                
                obj.S3CohortInfectedSnail.development_rate = obj.S3AgeRate.getRate(day);
                obj.S3CohortInfectedSnail.mortality_rate   = obj.S3Mortality.getRate(day);
                
                % METACERCARIAE
                
                obj.S4CohortMetacercariae.development_rate = obj.S4AgeRate.getRate(day);
                
                % calculate abundance of "matured" 
                
                S0egg = eggTimeSeries(day);
                
                matureTotal = S0egg * 1; % (can be changed to e.g. 0.01 to take 1% of eggs deposited on pasture per day)
                
                [matureTotal, removedMortality] = obj.S1CohortEgg.process_day(matureTotal);
                [matureTotal, removedMortality] = obj.S2CohortLarva.process_day(matureTotal);
                [matureTotal, removedMortality] = obj.S3CohortInfectedSnail.process_day(matureTotal);
                [matureTotal, removedMortality] = obj.S4CohortMetacercariae.process_day(matureTotal);

            end
            
            % output the abundance of infective metacercariae on pasture
            
            output_args = obj.S4SumMetacercariae;
            
            % remove temporary matrix allocation to reduce required memory

            obj.S1CohortEgg.free();
            obj.S2CohortLarva.free();
            obj.S3CohortInfectedSnail.free();
            obj.S4CohortMetacercariae.free();
            
        end
        
        %% plot model outputs
        
        function plotResult(obj)
            plot(obj.timeDimention, obj.S1SumEgg,'DisplayName','S1 Embrionic Egg');
            hold all;
            plot(obj.timeDimention, obj.S2SumLarva,'DisplayName','S2 Larva');
            plot(obj.timeDimention, obj.S3SumInfectedSnail,'DisplayName','S3 Infected Snail');
            plot(obj.timeDimention, obj.S4SumMetacercariae,'DisplayName','S4 Metacercariae');
            datetick('x', 'yy');
            axis([-inf,inf,-inf,inf]);
            hold off;
        end
        
        %% recover daily rates
        
        function result = getRateSeries(obj, rateProvider)
            result = zeros(size(obj.timeDimention));
            for day = 1:length(obj.timeDimention)
                result(day) = rateProvider.getRate(day);
            end
        end
        
    end
end

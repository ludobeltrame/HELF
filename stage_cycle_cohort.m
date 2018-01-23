classdef stage_cycle_cohort < handle
    
    properties
        
        stage_name

        mortality_rate;
        
        development_rate;
        
        progressionTrigger = 1; 
        % rationale behind use of the trigger: progress to next life-cycle
        % stage may be prevented until certain environmental conditions 
        % become suitable.
        
        developmentMortality = 0 
        % optional mortality mechanism (in case a trigger is used) for
        % individuals that have developed but cannot progress to next stage
        % due to triggering conditions not satisfied.
        
        maxQueueSize = 0;
        
        log_grouped_count = 0;
        
    end
        
    properties(SetAccess = private)
        
        rows
        cols
        
        development_distribution_used = 0; 
        development_distribution
        
        cohort_count
        cohort_development_age
        
    end
    
    properties(SetAccess = private, GetAccess = private)
        first_cohort_index = 1;
        last_cohort_index = 1;
    end
    
%%
    
    methods
        
        function obj = stage_cycle_cohort()
            obj.development_rate = 1; 
            obj.mortality_rate = 0;
        end
        
        function InitialiseDevelopmentDistribution(obj, paramA, paramB)
            % (even individuals exposed to the same environmental conditions
            % may develop at different times -> use of  weibull distribution 
            % to distribute development times)
            if(nargin < 3)
                paramA = 0.25;
                paramB = 1.5;
            end
            obj.development_distribution_used = 1;
            obj.development_distribution  = makedist('weibull', 'a', paramA, 'b', paramB);
        end
        
        function count = GetCount(obj)
            count = sum(obj.cohort_count);
        end
        
        function allocate(obj,maxQueueSize)
            obj.maxQueueSize = maxQueueSize;
            obj.cohort_count = zeros(maxQueueSize,1);
            obj.cohort_development_age = zeros(maxQueueSize,1);
        end
        
        function free(obj)
            obj.cohort_count = 0;
            obj.cohort_development_age = 0;
        end
        
        function [developed, removed_mortality] = process_day(obj, input)
            % this function progresses the stage under consideration by one
            % time step.
            % - input: number of individuals entering the stage, as they have developed from the previous stage
            % - developed: number of individuals in the stage that mature to the next one
            % - removed_mortality: number of individuals that die 
            
            removed_mortality = obj.increment_mortality();
            developed = obj.increment_development();
            removed_age = obj.increment_DevelopmentMortality();
            obj.input_from_previous(input); 
            removed_mortality = removed_mortality + removed_age;            
        end        
        
        function input_from_previous(obj, input)
            % this function adds individuals that have developed from the
            % previous stage.
            
            % once the count in a cohort goes below 0, stop processing that
            % group
            while(obj.first_cohort_index < obj.last_cohort_index ...
                    && obj.cohort_count(obj.first_cohort_index) < 1)
                obj.cohort_count(obj.first_cohort_index) = 0;
                obj.first_cohort_index = obj.first_cohort_index + 1;
            end
            
            if(obj.cohort_development_age(obj.last_cohort_index) == 0)
                % if development has not occurred, cohort can be grouped,
                % thus reducing processing time
                obj.log_grouped_count = obj.log_grouped_count + 1;
                obj.cohort_count(obj.last_cohort_index) = ...
                    obj.cohort_count(obj.last_cohort_index) + input;
            else 
                obj.last_cohort_index = obj.last_cohort_index + 1;
                active_count = obj.last_cohort_index - obj.first_cohort_index;
                
                if(active_count >= obj.maxQueueSize)
                    error('maxQueueSize reached');
                end
                
                if(obj.last_cohort_index >= obj.maxQueueSize)
                    % if the end of the cohort buffer size has been reached, 
                    % shift buffer to make room
                    n = -1 * (obj.first_cohort_index - 1);
                    obj.cohort_count = circshift(obj.cohort_count,n,1);
                    obj.cohort_development_age = circshift(obj.cohort_development_age,n,1);
                    obj.first_cohort_index = 1;
                    obj.last_cohort_index = active_count;
                end
                
                % add cohorts that have matured from previous stage
                obj.cohort_count(obj.last_cohort_index) = input;
                obj.cohort_development_age(obj.last_cohort_index) = 0;                
            end
        end
        
        function [developed_step] = increment_development(obj)
            % this function progresses all cohorts according to development
            % rate and returns grid counts that mature to next stage.
                        
            range = obj.first_cohort_index:obj.last_cohort_index;
            
            % update age of each cohort based on the development rate
            obj.cohort_development_age(range) = obj.cohort_development_age(range) + obj.development_rate;
            
            % calculate the total number of matured ("developed_step")
 
            if(obj.progressionTrigger == 1)

                if(obj.development_distribution_used == 1)
                    mask = cdf(obj.development_distribution, obj.cohort_development_age(range)-1);  
                else
                    mask = obj.cohort_development_age(range) >= 1;
                end
                
                % number of matured from each cohort
                matured = obj.cohort_count(range) .* mask;
                
                % number of individuals remaining in the current stage
                obj.cohort_count(range) = obj.cohort_count(range) - matured;
                
                % total number of matured individuals equal to the sum of
                % the matured cohorts
                developed_step = sum(matured);
            else
                developed_step = 0;
            end
            
        end
        
        function [removed] = increment_mortality(obj)
            % this function removes individuals from all cohorts according
            % to mortality rate and returns grid of death counts.
            
            range = obj.first_cohort_index:obj.last_cohort_index;
            
            if(obj.mortality_rate > 0) 
                
                removed_cohort = obj.cohort_count(range) .* obj.mortality_rate;
                obj.cohort_count(range) = obj.cohort_count(range) - removed_cohort;
                removed = sum(removed_cohort);
            else
                removed = 0;
            end
        end
        
        function [removed] = increment_DevelopmentMortality(obj)
            % this function removes individuals from all cohorts according
            % to developmentMortality and returns grid of death counts.
            
            % developmentMortality: additional optional mortality mechanism 
            % for cohorts that have reached developement but havent matured 
            % into the next life-cycle stage due to unfavourable triggering 
            % conditions (occurs when cohort age is above a certain threshold).
        
            range = obj.first_cohort_index:obj.last_cohort_index;
            
            if(obj.developmentMortality > 1)
                
                mask = obj.cohort_development_age(range) >= obj.developmentMortality;
                
                removed_cohort = obj.cohort_count(range) .* mask;
                obj.cohort_count(range) = obj.cohort_count(range) - removed_cohort;
                removed = sum(removed_cohort);
                
            else
                removed = 0;
            end
        end
        
    end
end

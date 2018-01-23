classdef rateCalcConstant < handle
    
    properties
        rate = 1;
    end
    
    methods
        function obj = rateCalcConstant(rate)
            obj.rate = rate;
        end
        
        function rate = getRate(obj, day)
            rate = obj.rate;
        end
    end
    
end

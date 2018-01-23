classdef rateCalcLinear < handle
   
    properties
        timeSeries
        
        a
        b
    end
    
    methods
        function obj = rateCalcLinear(a, b)
            obj.a = a; % intercept
            obj.b = b; % slope
        end
        
       function rate = getRate(obj, day)
           temp = obj.timeSeries(day);
           rate = obj.a + obj.b*temp;
           if(rate < 0)
               rate = 0;
           end
       end
    end
    
end
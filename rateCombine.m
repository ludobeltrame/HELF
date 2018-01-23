classdef rateCombine < handle
    
    properties
        
        inputA  
    
        inputB 
        
        rate   
    end
    
    methods
        function obj = rateCombine(inputA, inputB, rate)
            obj.inputA = inputA;
            obj.inputB = inputB;
            obj.rate = rate;
        end
        
        function rate = getRate(obj, day)
            valueA = obj.inputA.getRate(day);
            valueB = obj.inputB.getRate(day);
            
            rate = obj.rate * (valueA * valueB);
        end
    end
    
end

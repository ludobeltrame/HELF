classdef rateLinearFuzzyMembership < handle
    % rateLinearFuzzyMembership gives 'rateMember' when input value is 
    % within range (minB-maxA) and 'rateNonMember' when it is outside the
    % range, and uses a linear interpolation to switch between the 2 cases.
    %   
    %     minB     maxA
    %      |        |
    %       _________          _ rateMember
    %      /         \
    %     /           \
    % ___/             \_____  _ rateNonMember
    %    |             |
    %   minA          maxB
    
   properties
       timeSeries
       minA
       minB
       maxA
       maxB
       
       rateMember = 1;
       rateNonMember = 0;
    end
    
    methods
        function obj = rateLinearFuzzyMembership(minA, minB, maxA, maxB)
            obj.minA = minA;
            obj.minB = minB;
            obj.maxA = maxA;
            obj.maxB = maxB;
        end
        
        function rate = getRate(obj, day)
            value = obj.timeSeries(day);
            if(value > obj.minA && value < obj.maxB)
                if(value < obj.minB)
                    rate =  obj.rateNonMember + ((value - obj.minA) / (obj.minB - obj.minA) * (obj.rateMember - obj.rateNonMember));
                elseif(value > obj.maxA)
                    rate = obj.rateMember - ((value - obj.maxA) / (obj.maxB - obj.maxA) * (obj.rateMember - obj.rateNonMember));
                else
                    rate = obj.rateMember;
                end
            else
                rate = obj.rateNonMember;
            end
        end
    end
    
end


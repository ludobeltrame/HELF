function [sat_def, q_rout] = HELF_hydro_model_component (ppt, pet, TI, freq, param)
% this function calculates soil moisture for each TI class over a catchment,
% as well as streamflow at the catchment outlet, based concepts from the
% catchment-scale rainfall-runoff model TOPMODEL (*)
% INPUTS:
% - ppt  : array of size (# days in the sim period,1) containing rainfall 
% - pet  : array of size (# days in the sim period,1) containing potential evapotranspiration 
% - TI   : array of size (# of classes,1) containing Topographic Index values for each class
% - freq : array of size (# of classes,1) containing portion of catchment associated with each TI class
% - param: array of size (1,7) containing the model parameters (lnTe, m, td, Srz_init, Srz_max, a, b)
% OUTPUTS:
% - sat_def: matrix of size (# days in the sim period,# of classes) containing the simulated sat deficit
% - q_rout : array of size (# days in the sim period,1) containing the simulated total routed flow

% (*) references: 
% - Beven et al. (1995) TOPMODEL. In: Sing VP (Ed), Computer Models of Watershed Hydrology. Water Resources Publications, Colorado. pp. 627-668.

%% set up

% recover model parameters:
lnTe     = param(1); 
m        = param(2); 
td       = param(3); 
Srz_init = param(4); 
Srz_max  = param(5); 
a        = param(6); 
b        = param(7);

% number of time samples (i.e. # of days in the sim period):
T = length(ppt);

% number of space samples (i.e. # of TI classes considered):
Sp = length(TI);

% initialise variables:
q_tot = zeros(T,1);          % catchment average total flow
q_rout = zeros(T,1);         % routed catchment average total flow
q_over = zeros(T,1);         % catchment average overland flow
overland_flow = zeros(T,Sp); % local contribution to overland flow
q_vertical = zeros(T,1);     % catchment average vertical flow
vertical_flow = zeros(T,Sp); % local vertical flow 
actual_et = zeros(T,Sp);     % local actual evapotranspiration
q_sub = zeros(T,1);          % catchment average subsurface flow
q_inf = zeros(T,1);          % catchment average infiltration
Srz = zeros(T,Sp);           % local root zone storage deficit
Suz = zeros(T,Sp);           % local unsaturated zone storage
Saverage = zeros(T,1);       % catchment average storage deficit
S = zeros(T,Sp);             % local saturated zone deficit (calculations)
sat_def = zeros(T,Sp);       % local saturated zone deficit (model output)

% recharge:
r = max(mean(ppt)-mean(pet), 0.1);

% catchment average topographic index:
lambda = sum(TI.*freq)/sum(freq);

% initial storage deficit in the root zone:
Srz(1,:) = Srz_init;

% initial catchment average storage deficit:
sbar = -m*(log(r/exp(lnTe))+lambda); 
Saverage(1) = sbar;

%% calculations

% loop through time:

for t = 1:T
    
    % calculate catchment average subsurface flow
    q_sub(t) = exp(lnTe-lambda)*exp(-Saverage(t)/m); 
    
    % calculate catchment average infiltration 
    % (by default all rain infiltrates)
    q_inf(t) = ppt(t); 
    
    % loop through space:
    
    for i = 1:Sp
        
        % calculate and store local storage deficit
        S(t,i) = Saverage(t)+m*(lambda-TI(i));
                sat_def(t,i) = S(t,i);   
        
        if S(t,i)<0
            S(t,i)=0;
        end
        
        % rain first enters the root zone (written as deficit) 
        % (root zone storage has here the same function as the interception
        % store in other topmodel versions)
        Srz(t,i) = Srz(t,i)-q_inf(t);
            
        % only when root zone is filled water goes to the unsaturated zone
        % (written as storage)
        if Srz(t,i)<0
            Suz(t,i)=Suz(t,i)-Srz(t,i);
            Srz(t,i)=0;
        end
        
        % if Suz exceeds S, then saturation-excess overland flow occurs
        % (S determines the storage capacity of Suz)
        if Suz(t,i)>S(t,i)
            overland_flow(t,i)=Suz(t,i)-S(t,i);
            Suz(t,i)=S(t,i);
        end
          
        % calculate local vertical flow entering the saturated zone 
        % (use of a storage deficit-dependent time delay)
        if S(t,i)>0 
            vertical_flow(t,i)=Suz(t,i)/(S(t,i)*td); 
            if vertical_flow(t,i)>Suz(t,i) 
                vertical_flow(t,i)=Suz(t,i);
            end
            Suz(t,i)=Suz(t,i)-vertical_flow(t,i); 
            if Suz(t,i)<0 
                Suz(t,i)=0;
            end  
        end
        
        % calculate local actual evapotranspiration 
        actual_et(t,i) = 0;
        if pet(t)>0
            actual_et(t,i)=pet(t)*(1-Srz(t,i)/Srz_max); 
            if actual_et(t,i)>Srz_max-Srz(t,i)
                actual_et(t,i)=Srz_max-Srz(t,i);
            end
        end
        Srz(t,i)=Srz(t,i)+actual_et(t,i); 
          
    end
    
    % catchment average vertical flow
    q_vertical(t) = sum((freq').*vertical_flow(t,:))/sum(freq);
    
    % update catchment average sat def through catchment mass balance
    Saverage(t+1)=Saverage(t)+q_subsurface(t)-q_vertical(t);
    
    Srz(t+1,:) = Srz(t,:);
    Suz(t+1,:) = Suz(t,:);
    
    % catchment average overland flow
    q_over(t) = sum((freq').*overland_flow(t,:))/sum(freq);
    
    % catchment average total flow
    q_tot(t) = q_sub(t)+q_over(t);
    
end

% derive time delay histogram based on gamma distribution
size_frac_future=20; 
x = 1:size_frac_future;
frac_future=gampdf(x,a,b); 

% distribute each daily streamflow over predefined number of future time
% steps, according to the predefined fractions 
future = zeros(size_frac_future,1); 
for t=1:T
    
    for j=1:size_frac_future
        future(j)=future(j)+q_tot(t)*frac_future(j); 
    end

    q_rout(t)=future(1);
    future(1:size_frac_future-1)=future(2:size_frac_future);
    future(size_frac_future)=0;

end


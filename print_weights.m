function print_weights(ACFT)
% Function:
%   print_weights
%
% Description: 
%   Prints the weights of the retrofitted additions of the 
%
% Input:
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:

% Unwrap structure
Weight      = ACFT.Weight;
Baseline    = Weight.Baseline;
IntParallel = Weight.IntParallel;
Parallel    = Weight.Parallel;
Series      = Weight.Series;
SeriesT     = Weight.SeriesT;
Hydrogen    = Weight.Hydrogen;

% Start printing
fprintf('%-15s\n','Weight Summary:')
% Header
fprintf('%15s %10s %10s %10s %10s %10s %10s \n','Mass (kg)','Baseline','Int. Para.','Parallel','Series','SeriesT','Hydrogen');
% fuel engine
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','fuel Engine',Baseline.fuelEngine_kg,IntParallel.fuelEngine_kg,Parallel.fuelEngine_kg,...
    Series.fuelEngine_kg,SeriesT.fuelEngine_kg,Hydrogen.fuelEngine_kg);
% electric Motor
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Elec Motor',Baseline.electricMotor_kg,IntParallel.electricMotor_kg,Parallel.electricMotor_kg,...
    Series.electricMotor_kg,SeriesT.electricMotor_kg,Hydrogen.electricMotor_kg);
% battery
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Battery',Baseline.battery_kg,IntParallel.battery_kg,Parallel.battery_kg,...
    Series.battery_kg,SeriesT.battery_kg,Hydrogen.battery_kg);
% Electric systems
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Elec Sys.',Baseline.electricSystems_kg,IntParallel.electricSystems_kg,Parallel.electricSystems_kg,...
    Series.electricSystems_kg,SeriesT.electricSystems_kg,Hydrogen.electricSystems_kg);
% Fuel Tank
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Fuel Tank',Baseline.fuelTank_kg,IntParallel.fuelTank_kg,Parallel.fuelTank_kg,...
    Series.fuelTank_kg,SeriesT.fuelTank_kg,Hydrogen.fuelTank_kg);
% Total weight of propulsion system
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Total',Baseline.totalPropulsion_kg,IntParallel.totalPropulsion_kg,Parallel.totalPropulsion_kg,...
    Series.totalPropulsion_kg,SeriesT.totalPropulsion_kg,Hydrogen.totalPropulsion_kg);
fprintf('\n');
end
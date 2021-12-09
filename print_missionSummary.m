function print_missionSummary(MissionPerformance)
% Function:
%   print_missionSummary
%
% Description: 
%   Prints the total mission fuel consumption and battery energy
%   consumption for the baseline, series and parallel architecture as well
%   as the % reduction in fuel consumption.
%
% Input:
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:
%
% Unwrap structure
baselineFuelConsMission_USgal  = MissionPerformance.Baseline.fuelConsMission_USgal;
parallelFuelConsMission_USgal  = MissionPerformance.Parallel.fuelConsMission_USgal;
seriesFuelConsMission_USgal    = MissionPerformance.Series.fuelConsMission_USgal;
baselineBattEnergyMission_kWh  = 0;
parallelBattEnergyMission_kWh  = MissionPerformance.Parallel.electricEnergyMission_kWh;
seriesBattEnergyMission_kWh    = MissionPerformance.Series.battEnergyMission_kWh;

% Start printing
fprintf('%-15s\n','Mission Summary:')
% Header
fprintf('%15s %10s %10s %10s \n','Retrofit','Baseline','Series','Parallel');
% Mission fuel consumption
fprintf('%15s %10.1f %10.1f %10.1f\n','fuel (US Gal)',baselineFuelConsMission_USgal,seriesFuelConsMission_USgal,...
    parallelFuelConsMission_USgal);
% Percentage change in fuel consumption relative to baseline
fuelConsMission_USgal = [baselineFuelConsMission_USgal,seriesFuelConsMission_USgal,...
    parallelFuelConsMission_USgal];
fuelPercentageDiff    = (fuelConsMission_USgal-baselineFuelConsMission_USgal)./baselineFuelConsMission_USgal*100;
fprintf('%15s %10.0f %10.0f %10.0f \n','Fuel % Diff.',fuelPercentageDiff(1),fuelPercentageDiff(2),...
    fuelPercentageDiff(3));
% Battery Energy consumption
fprintf('%15s %10.1f %10.1f %10.1f\n','BattEnergy(kWh)',baselineBattEnergyMission_kWh,...
   seriesBattEnergyMission_kWh,parallelBattEnergyMission_kWh);
fprintf('\n');
end
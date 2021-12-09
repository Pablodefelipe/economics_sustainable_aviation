function print_missionPerformanceParallel(Mission,MissionPerformance)
% Function:
%   print_missionPerformaceParallel
%
% Description: 
%   Prints the average ICE and electric power, fuel consumption and energy consumption for a particular mission
%   profile 
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:

fprintf('%-15s\n','Parallel:')
% Header
fprintf('%15s %10s %10s %10s %10s %10s %10s\n','Segment','Taxi','TO','Climb','Cruise','Descent','Reserve');
% time and distance
% Pass all times to minutes
timeTaxi_min    = Mission.time.taxi_s/60;
timeTakeOff_min = Mission.time.takeOff_s/60;
timeClimb_min   = Mission.time.climb_s/60;
timeCruise_min  = Mission.time.cruise_s/60;
timeDescent_min = Mission.time.descent_s/60;
timeReserve_min = (Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s)/60;

% Unwrap structure
Power          = MissionPerformance.Parallel.ICE_power_kW;
Fuel           = MissionPerformance.Parallel.fuelCons_USgal;
ElecPower      = MissionPerformance.Parallel.electricPower_kW;
ElecEnergy     = MissionPerformance.Parallel.electricEnergy_kWh;
BattEnergy     = MissionPerformance.Parallel.batteryEnergy_kWh;
Parallel       = MissionPerformance.Parallel;
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','time (min)',timeTaxi_min,timeTakeOff_min,...
    timeClimb_min,timeCruise_min,timeDescent_min,timeReserve_min);

distanceTaxi_nm    = 0;
distanceTakeOff_nm = 0;
distanceClimb_nm   = Mission.distance.climb_m/1852;
distanceCruise_nm  = Mission.distance.cruise_m/1852;
distanceDescent_nm = Mission.distance.descent_m/1852;
distanceReserve_nm = Mission.distance.reserve_m/1852;

fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','distance (nm)',distanceTaxi_nm,distanceTakeOff_nm,...
    distanceClimb_nm,distanceCruise_nm,distanceDescent_nm,distanceReserve_nm);
% Average Power
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Avg.Power(kW)',Power.taxi,Power.takeOff,...
    Power.climb,Power.cruise,Power.descent,Power.reserve);
% Fuel consumption
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Fuel (USgal)',Fuel.taxi,Fuel.takeOff,...
    Fuel.climb,Fuel.cruise,Fuel.descent,MissionPerformance.Parallel.fuelConsReserve_USgal);
% Electric Average Power
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Elec.Power(kW)',ElecPower.taxi,ElecPower.takeOff,...
    ElecPower.climb,ElecPower.cruise,ElecPower.descent,ElecPower.reserve);
% Electric Energy consumption
fprintf('%15s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n','ElecEnergy(kWh)',ElecEnergy.taxi,ElecEnergy.takeOff,...
    ElecEnergy.climb,ElecEnergy.cruise,ElecEnergy.descent,Parallel.electricEnergyReserve_kWh);
% Battery initial energy
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Batt.State(kWh)',BattEnergy.taxi,BattEnergy.takeOff,...
    BattEnergy.climb,BattEnergy.cruise,BattEnergy.descent,BattEnergy.reserveClimb);
fprintf('\n');
end
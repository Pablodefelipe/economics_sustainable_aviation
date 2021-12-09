function print_missionPerformanceSeries(Mission,MissionPerformance)
% Function:
%   print_missionPerformaceSeries
%
% Description: 
%   Prints the average ICE and electric power, fuel consumption and energy consumption for a particular mission
%   profile for the series architecture 
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:

fprintf('%-15s\n','Series:')
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
ICEPower          = MissionPerformance.Series.ICE_power_kW;
Fuel           = MissionPerformance.Series.fuelCons_USgal;
ElecPower      = MissionPerformance.Series.ElectricPower_kW;
BattEnergy     = MissionPerformance.Series.BattEnergy_kWh;
GenEnergy      = MissionPerformance.Series.GenEnergy_kWh;
BatteryEnergy  = MissionPerformance.Series.batteryEnergy_kWh;
Series         = MissionPerformance.Series;
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','time (min)',timeTaxi_min,timeTakeOff_min,...
    timeClimb_min,timeCruise_min,timeDescent_min,timeReserve_min);

distanceTaxi_nm    = 0;
distanceTakeOff_nm = 0;
distanceClimb_nm   = Mission.distance.climb_m/1852;
distanceCruise_nm  = Mission.distance.cruise_m/1852;
distanceDescent_nm = Mission.distance.descent_m/1852;
distanceReserve_nm = '~';

fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10s\n','distance (nm)',distanceTaxi_nm,distanceTakeOff_nm,...
    distanceClimb_nm,distanceCruise_nm,distanceDescent_nm,distanceReserve_nm);
% Average Electric Power
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Elec.Power(kW)',ElecPower.taxi,ElecPower.takeOff,...
    ElecPower.climb,ElecPower.cruise,ElecPower.descent,ElecPower.reserve);
% Average ICE Power
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','ICE. Power(kW)',ICEPower.taxi,ICEPower.takeOff,...
    ICEPower.climb,ICEPower.cruise,ICEPower.descent,ICEPower.reserve);
% Fuel consumption
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Fuel (USgal)',Fuel.taxi,Fuel.takeOff,...
    Fuel.climb,Fuel.cruise,Fuel.descent,MissionPerformance.Series.fuelConsReserve_USgal);
% Generator Energy
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Gen.Energy(kWh)',GenEnergy.taxi,GenEnergy.takeOff,...
    GenEnergy.climb,GenEnergy.cruise,GenEnergy.descent,Series.genEnergyReserve_kWh);
% Battery Energy consumption
fprintf('%15s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n','BattEnergy(kWh)',BattEnergy.taxi,BattEnergy.takeOff,...
    BattEnergy.climb,BattEnergy.cruise,BattEnergy.descent,Series.battEnergyReserve_kWh);
% Battery initial energy
fprintf('%15s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n','Batt.State(kWh)',BatteryEnergy.taxi,BatteryEnergy.takeOff,...
    BatteryEnergy.climb,BatteryEnergy.cruise,BatteryEnergy.descent,BatteryEnergy.reserveClimb);
fprintf('\n');
end
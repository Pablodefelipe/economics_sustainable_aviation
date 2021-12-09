function print_missionPerformanceBaseline(Mission,MissionPerformance)
% Function:
%   print_missionPerformaceBaseline
%
% Description: 
%   Prints the average power and fuel consumption for a particular mission
%   profile 
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:

fprintf('%-15s\n','Baseline:')
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

% Unwrap power
Power = MissionPerformance.Baseline.ICE_power_kW;
Fuel  = MissionPerformance.Baseline.fuelCons_USgal;
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
    Fuel.climb,Fuel.cruise,Fuel.descent,MissionPerformance.Baseline.fuelConsReserve_USgal);
fprintf('\n');
end

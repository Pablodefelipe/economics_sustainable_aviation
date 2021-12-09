function writeExcel_mission(Mission,MissionPerformance)
% Function:
%   writeExcel_mission
%
% Description: 
%   Exports the mission parameters
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:

% time and distance
% Pass all times to minutes
timeTaxi_min    = Mission.time.taxi_s/60;
timeTakeOff_min = Mission.time.takeOff_s/60;
timeClimb_min   = Mission.time.climb_s/60;
timeCruise_min  = Mission.time.cruise_s/60;
timeDescent_min = Mission.time.descent_s/60;
timeReserve_min = (Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s)/60;
distanceTaxi_nm    = 0;
distanceTakeOff_nm = 0;
distanceClimb_nm   = Mission.distance.climb_m/1852;
distanceCruise_nm  = Mission.distance.cruise_m/1852;
distanceDescent_nm = Mission.distance.descent_m/1852;
distanceReserve_nm = Mission.distance.reserve_m/1852;

distanceTaxi_km    = 0;
distanceTakeOff_km = 0;
distanceClimb_km   = Mission.distance.climb_m/1000;
distanceCruise_km  = Mission.distance.cruise_m/1000;
distanceDescent_km = Mission.distance.descent_m/1000;
distanceReserve_km = Mission.distance.reserve_m/1000;

altitudeCruise_ft  = Mission.altitude.cruise_m*3.28;
altitudeReserve_ft = Mission.altitude.reserve_m*3.28;

climbSpeed_kts = MissionPerformance.climbSpeed_ms*3600/1852;
cruiseSpeed_kts= MissionPerformance.cruiseSpeed_ms*3600/1852;
reserveSpeed_kts=MissionPerformance.reserveSpeed_ms*3600/1852;
descentSpeed_kts=MissionPerformance.descentSpeed_ms*3600/1852;

climbRate_fpm = MissionPerformance.climbRate_ms*60*3.28;
descentRate_fpm=MissionPerformance.descentRate_ms*60*3.28;

filename = 'exampleData.xlsx';
Segment  = ["Taxi";"TakeOff";"Climb";"Cruise ";"Descent";"Reserve"];
Time_min    = [timeTaxi_min;timeTakeOff_min;timeClimb_min;timeCruise_min;...
    timeDescent_min;timeReserve_min];
Distance_km  = [distanceTaxi_km;distanceTakeOff_km;distanceClimb_km;...
    distanceCruise_km ;distanceDescent_km;distanceReserve_km];
Altitude_ft    = ["0";"0";"Climbing" ;altitudeCruise_ft;"Descending";...
    altitudeReserve_ft];
Velocity_kts   = ["0";"Accelerating";climbSpeed_kts;cruiseSpeed_kts;...
    descentSpeed_kts;reserveSpeed_kts];
RoC_fpm  = ["0";"0";climbRate_fpm;"0";-descentRate_fpm;"0"];


mission = table(Segment,Time_min,Distance_km,Altitude_ft,Velocity_kts,RoC_fpm);

writetable(mission,filename,'Sheet',1,'Range','A14')
end

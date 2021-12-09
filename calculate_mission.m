function Mission = calculate_mission(Mission,MissionPerformance)
% Function:
%   calculate_missionPerformace
%
% Description: 
%   Calculate distance and time spent in each segment of the mission
%
% Input:
%   Mission             - Deffinition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:
%   Mission            - Mission with all parameters completed

% Taxi and TO already calculated
% Climb
if isfield(Mission,'timeClimb_s') == 1
    % If time prescribed calculate distance
    Mission.distance.climb_m = MissionPerformance.climbSpeed_ms*cos(MissionPerformance.climbAngle_rad)*Mission.time.climb_s;
    Mission.altitude.climb_m = MissionPerformance.climbRate_ms*Mission.time.climb_s;
else
    % The altitude is prescribed
    Mission.time.climb_s     = Mission.altitude.cruise_m/MissionPerformance.climbRate_ms;
    Mission.distance.climb_m = MissionPerformance.climbSpeed_ms*cos(MissionPerformance.climbAngle_rad)*Mission.time.climb_s;
    
end

% Cruise
if isfield(Mission.time,'cruise_s') == 1
    % Calculate the distance
    Mission.distance.cruise_m = MissionPerformance.cruiseSpeed_ms*Mission.time.cruise_s;
else
    Mission.time.cruise_s     = Mission.distance.cruise_m/MissionPerformance.cruiseSpeed_ms;
end

% Descent
Mission.time.descent_s = Mission.altitude.cruise_m/MissionPerformance.descentRate_ms;
Mission.distance.descent_m = MissionPerformance.descentSpeed_ms*cos(asin(MissionPerformance.descentRate_ms/MissionPerformance.descentSpeed_ms))...
    *Mission.time.descent_s;
    
% Reserve climb loiter and descent time and distance
Mission.time.reserveClimb_s   =  Mission.altitude.reserve_m/MissionPerformance.climbRate_ms;
Mission.distance.reserveClimb_m = MissionPerformance.climbSpeed_ms*cos(MissionPerformance.climbAngle_rad)*Mission.time.reserveClimb_s;
Mission.distance.reserveLoiter_m = MissionPerformance.reserveSpeed_ms*Mission.time.reserveLoiter_s;
Mission.time.reserveDescent_s =  Mission.altitude.reserve_m/MissionPerformance.descentRate_ms;
Mission.distance.reserveDescent_m = MissionPerformance.descentSpeed_ms*cos(asin(MissionPerformance.descentRate_ms/MissionPerformance.descentSpeed_ms))...
    *Mission.time.descent_s;
Mission.distance.mission_m = Mission.distance.climb_m+...
    Mission.distance.cruise_m+Mission.distance.descent_m;
Mission.distance.reserve_m = Mission.distance.reserveClimb_m + Mission.distance.reserveLoiter_m + Mission.distance.reserveDescent_m;
% Calculate the density for cruise and reserve loiter segments
[~,~,~,Mission.Rho_kgm3.cruise] = atmosisa(Mission.altitude.cruise_m); 
[~,~,~,Mission.Rho_kgm3.reserveLoiter] = atmosisa(Mission.altitude.reserve_m);
end


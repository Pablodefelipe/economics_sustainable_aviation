function [timeTO_s] = calculate_takeOff(MissionPerformance)
% Function:
%   calculate_takeOff(MissionPerformance)
%
% Description: 
%       Function to calculate the takeoff time assuming constant
%       acceleration.
% Input:
%  MissionPerformance.liftOffSpeed_ms   - Lift Off Speed in m/s
%  MissionPerformance.takeOffDistance_m - Take Off distance to clear 50ft
%
% Outputs:
%  timeTO_m     - Total time in s required to reach 50 ft 

liftOffSpeed_ms = MissionPerformance.liftOffSpeed_ms;
takeOffDistance_m = MissionPerformance.takeOffDistance_m;

timeTO_s = ceil(2*takeOffDistance_m/liftOffSpeed_ms); 

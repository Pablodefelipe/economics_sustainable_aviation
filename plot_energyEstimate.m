function [] = plot_energyEstimate(MissionPerformance)
% Function:
%   calculate_energyEstimate
%
% Description: 
%   Calculate energy estimate for the mission
%
% Input:
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:
%   Energy_kWh              - Mission with all parameters completed

% Unwrap the structure
Energy_kWh = MissionPerformance.EnergyEstimate_kWh;

taxi = Energy_kWh.taxi;       
takeOff = Energy_kWh.takeOff;    
climb  = Energy_kWh.climb;    
cruise  = Energy_kWh.cruise;    
descent = Energy_kWh.descent; 
totalMission = Energy_kWh.totalMission;
reserveClimb = Energy_kWh.reserveClimb; 
reserveLoiter = Energy_kWh.reserveLoiter; 
reserveDescent = Energy_kWh.reserveDescent;


% Plot a pie chart
X = [taxi,takeOff,climb,cruise,descent];
percentage = X/totalMission*100
labels = {'taxi (2%)','TakeOff (4%)','climb (27%)','cruise (59%)','descent (8%)'};
figure()
pie(X,labels) 
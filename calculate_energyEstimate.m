function Energy_kWh = calculate_energyEstimate(Mission,MissionPerformance)
% Function:
%   calculate_energyEstimate
%
% Description: 
%   Calculate energy estimate for the mission
%
% Input:
%   Mission             - Deffinition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:
%   Energy_kWh              - Mission with all parameters completed

% Unwrap power and time
pTaxi       =  MissionPerformance.ShaftPower_kW.taxi;    
pTakeOff    = MissionPerformance.ShaftPower_kW.takeOff;
pClimb      =  MissionPerformance.ShaftPower_kW.climb;
pCruise     = MissionPerformance.ShaftPower_kW.cruise;  
pDescent    = MissionPerformance.ShaftPower_kW.descent;
pReserveClimb   = MissionPerformance.ShaftPower_kW.reserveClimb; 
pReserveLoiter  = MissionPerformance.ShaftPower_kW.reserveLoiter;
pReserveDescent =  MissionPerformance.ShaftPower_kW.reserveDescent; 

tTaxi       =  Mission.time.taxi_s;    
tTakeOff    = Mission.time.takeOff_s;
tClimb      =  Mission.time.climb_s;
tCruise     = Mission.time.cruise_s;  
tDescent    = Mission.time.descent_s;
tReserveClimb   = Mission.time.reserveClimb_s; 
tReserveLoiter  = Mission.time.reserveLoiter_s;
tReserveDescent =  Mission.time.reserveDescent_s;

Energy_kWh.taxi       =  pTaxi*tTaxi/3600;    
Energy_kWh.takeOff    =  pTakeOff*tTakeOff/3600;    
Energy_kWh.climb      =  pClimb*tClimb/3600;    
Energy_kWh.cruise     =  pCruise*tCruise/3600;    
Energy_kWh.descent    =  pDescent*tDescent/3600; 
Energy_kWh.reserveClimb   =  pReserveClimb*tReserveClimb/3600; 
Energy_kWh.reserveLoiter  =  pReserveLoiter*tReserveLoiter/3600; 
Energy_kWh.reserveDescent =  pReserveDescent*tReserveDescent/3600;

Energy_kWh.totalMission = Energy_kWh.taxi+Energy_kWh.takeOff+Energy_kWh.climb+Energy_kWh.cruise+Energy_kWh.descent;
Energy_kWh.totalReserve = Energy_kWh.reserveClimb+Energy_kWh.reserveLoiter+Energy_kWh.reserveDescent;



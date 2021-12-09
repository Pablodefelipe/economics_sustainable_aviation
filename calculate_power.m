function [power_kW] = calculate_power(ACFT,weight_N,airspeed_ms,rho_kgm3)
% Function:
%   calculate_power(ACFT,weight_N,airspeed_ms,rho_kgm3)
%
% Description: 
%       Function to calculate the power required at the shaft for steady levelcruise.
% Input:
%  ACFT            - Aircraft parameters
%  weight_N        - Weight of aircraft in N
%  airspeed_ms     - Airspeed in m
%  rho_kgm3        - Density in kg/m^3
% Outputs:
%  power_kW        - Power required at the shaft in kW

propEff = ACFT.Propulsion.propEff;
wingArea_m2 = ACFT.wingArea_m2;
cd0 = ACFT.Aero.cd0;
k   = ACFT.Aero.k;

lift_N = weight_N;
cl     = lift_N/(0.5*rho_kgm3*airspeed_ms^2*wingArea_m2);
cd     = cd0 + k*cl^2;
LoD    = cl/cd;
drag_N = 0.5*rho_kgm3*airspeed_ms^2*wingArea_m2*cd;
power_kW   = drag_N*airspeed_ms/propEff/1000 ;         
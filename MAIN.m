% Main 
% Script to perform cost estimate analysis for a generic aircraft 
% Date: 10/11/2021

clear
clc
close all

% Fuel densities
rhoJETA_kgm3 = 840;
rhoAVGAS_kgm3= 720;
rhoLH2_kgm3  = 70.85;
%% Define Aircraft and its retrofits
% Could we use a define function? We could
% Inputs from POH
ACFT.name = 'PA-31';
% Define general parameters
ACFT.wingArea_m2 = 229/(3.28^2);      % Reference wing area in m^2
ACFT.aspectRatio = 7.22; 
ACFT.engineNum   = 2;
ACFT.maxFuelVolume_m3  = 190*0.003785;% Maximum Volume inside the original tanks
ACFT.Weight.MTOW_kg    = 6500/2.205;  % MTOW of aircraft
ACFT.Weight.EW_kg      = 4003/2.205;  % Operational empty weight of aircraft, inluding unusable fuel, oil and full operating fluids

% Define performance parameters

% Aero parameters: cd0 and k
ACFT.Aero.cd0 = 0.027;
ACFT.Aero.k   = 1/(pi*ACFT.aspectRatio*0.8);

% Own inputs
% Define weights that are constant for each aircraft and weight parameters
% to calculate weights of the retrofits
ACFT.Weight.electricMotorSpecificPower_kWkg  = 3.1;
ACFT.Weight.batterySpecificEnergy_kWhkg      = 0.25;
ACFT.Weight.icePistonSpecificPower_kWkg      = 0.9;
ACFT.Weight.iceDieselSpecificPower_kWkg      = 0.8;
ACFT.Weight.iceTurbineSpecificPower_kWkg     = 3.0;
ACFT.Weight.gearboxSpecificPower_kWkg        = 3.0;
ACFT.Weight.fuelCellSpecificPower_kWkg       = 1.6;
ACFT.Weight.electricSysSpecificPower_kWkg    = 5;       % Obtained from the MAGNIX Drive
% ACFT.Weight.hydrogenFuelTankSpecificEnergy_kWhkg = 2.38; % This is only the hydrogen fuel tank (Should include the LH2 fuel)
% ACFT.Weight.hydrogenPlusFuelSpecificEnergy_kWhkg = 2.22; % Hydrogen plus fuel
ACFT.Weight.hydrogenFuelTankAreaDensity_kgm2 = 75;

ACFT.Weight.crew_kg                         = 100;
ACFT.Weight.passenger_kg                    = 100;  % Includes luggage
% Prescribe number of passengers for each retrofit
ACFT.Weight.passenger_n                     = 7;    % Same for all retrofits
ACFT.Weight.payload_kg = ACFT.Weight.passenger_kg*ACFT.Weight.passenger_n;

ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg = 12.1;
ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg = 1/0.03;
ACFT.EnergyStorage.Electric.batterySpecificEnergy_kWhkg = [0.4, 0.6];

ACFT.Cost.avgas_usCentskWh = 31.3;
ACFT.Cost.hydrogen_usCentskWh = 50;
ACFT.Cost.electricity_usCentskWh = 10.4;

%% Define the powertrain parameters
ACFT.Propulsion.maxPower_kW = 2*310*0.7457;         % Maximum installed power

ACFT.Propulsion.propEff     = 0.8;
%% Mission Definition
% Define Mission performance parameters based on POH
MissionPerformance.liftOffSpeed_ms   = 83*1852/3600;
MissionPerformance.takeOffDistance_m = 2300/3.28;
MissionPerformance.climbRate_ms      = 1300/(3.28*60);
MissionPerformance.climbSpeed_ms     = 84*1852/3600;
MissionPerformance.climbAngle_rad    = asin(MissionPerformance.climbRate_ms/MissionPerformance.climbSpeed_ms);
MissionPerformance.cruiseSpeed_ms    = 154*1852/3600;
MissionPerformance.reserveSpeed_ms   = 130*1852/3600;
MissionPerformance.descentRate_ms    = 700/(3.28*60);
MissionPerformance.descentSpeed_ms   = 140*1852/3600; % Descent Airspeed
% Define mission 
% Taxi
Mission.time.taxi_s = 5*60;
% TO 
Mission.time.takeOff_s   = calculate_takeOff(MissionPerformance);
% Climb
Mission.altitude.cruise_m = 5000/3.28; 
[~, ~, ~, rhoCruise_kgm3] = atmosisa(Mission.altitude.cruise_m);
% Cruise
Mission.time.cruise_s = 15*60;
% Descent - Calculated with calculate_mission function
% Reserves 
Mission.altitude.reserve_m = 2000/3.28;
[~, ~, ~, rhoReserve_kgm3] = atmosisa(Mission.altitude.reserve_m);
Mission.time.reserveLoiter_s = 45*60;
% Calculate mission parameters required
Mission = calculate_mission(Mission,MissionPerformance); %calculate mission parameters
% Define Mission Performance shaft powers (Only calculate cruise shaft power)
MissionPerformance.ShaftPower_kW.taxi    = 0.1*ACFT.Propulsion.maxPower_kW;
MissionPerformance.ShaftPower_kW.takeOff = ACFT.Propulsion.maxPower_kW;
MissionPerformance.ShaftPower_kW.climb   = ACFT.Propulsion.maxPower_kW;
MissionPerformance.ShaftPower_kW.cruise  = calculate_power(ACFT,ACFT.Weight.MTOW_kg*9.81,MissionPerformance.cruiseSpeed_ms,rhoCruise_kgm3);
MissionPerformance.ShaftPower_kW.descent = 0.15*ACFT.Propulsion.maxPower_kW;
MissionPerformance.ShaftPower_kW.reserveClimb = ACFT.Propulsion.maxPower_kW;
MissionPerformance.ShaftPower_kW.reserveLoiter = calculate_power(ACFT,ACFT.Weight.MTOW_kg*9.81,MissionPerformance.reserveSpeed_ms,rhoReserve_kgm3);
MissionPerformance.ShaftPower_kW.reserveDescent = 0.15*ACFT.Propulsion.maxPower_kW;

MissionPerformance.EnergyEstimate_kWh = calculate_energyEstimate(Mission,MissionPerformance);
%% Size specific powertrain powers
generatorEff       = 0.9;
electricMotorEff   = 0.95;
batteryEff         = 0.9;
ACFT.Propulsion.hybridisationRatioPower = 1-MissionPerformance.ShaftPower_kW.cruise/ACFT.Propulsion.maxPower_kW;
ACFT.Propulsion.hybridisationRatioEnergy = 0.05;     % Use this as a starting point                 
ACFT.Propulsion.Series.ICEsurplusPower_kW = 0;     % Extra Power to size the ICE of the series
ACFT.Propulsion.Hydrogen.fuelCellOversizeFactor = 1/0.75;   % Value to oversize the fuel cell

% Define all the retrofits maximum installed power
ACFT.Propulsion.Baseline.maxICE_Power_kW     = ACFT.Propulsion.maxPower_kW;
% ACFT.Propulsion.Baseline.ICE_n               = 2;

% Parallel 
ACFT.Propulsion.Parallel.maxICE_Power_kW= (1-...
    ACFT.Propulsion.hybridisationRatioPower)*ACFT.Propulsion.maxPower_kW;
    
% ACFT.Propulsion.Parallel.ICE_n               = 2;

ACFT.Propulsion.Parallel.maxElectricPower_kW = ACFT.Propulsion.maxPower_kW-ACFT.Propulsion.Parallel.maxICE_Power_kW;
% ACFT.Propulsion.IntParallel.electric_n          = 2;

% Series
ACFT.Propulsion.Series.maxICE_Power_kW=(1-...
    ACFT.Propulsion.hybridisationRatioPower)*ACFT.Propulsion.maxPower_kW/...
    (electricMotorEff*generatorEff);% Maximum power of the ICE component of the generator
ACFT.Propulsion.Series.maxElectricPower_kW   = ACFT.Propulsion.maxPower_kW;
% ACFT.Propulsion.Series.electric_n              = 2;
% ACFT.Propulsion.Series.electricPropDiameter_m = 2.0;
% ACFT.Propulsion.Series.maxElectricPowerOEI_kW = 0;     % DEPENDS ON BATTERY Assuming 6C discharge rate
% ACFT.Propulsion.Series.electricBladeN   = 2;

% Hydrogen Fuel Cell
% In order to design depends on the mission
ACFT.Propulsion.Hydrogen.maxFC_Power_kW       = (1-...
    ACFT.Propulsion.hybridisationRatioPower)*ACFT.Propulsion.maxPower_kW*...
    ACFT.Propulsion.Hydrogen.fuelCellOversizeFactor; % Fuel cell will only operate at a maximum of 75% its max power
ACFT.Propulsion.Hydrogen.maxElectricPower_kW  = ACFT.Propulsion.maxPower_kW;
% ACFT.Propulsion.Hydrogen.electric_n           = 2;
% ACFT.Propulsion.Hydrogen.electricPropDiameter_m = 2.0;
% ACFT.Propulsion.Hydrogen.maxElectricPowerOEI_kW = 0;     % DEPENDS ON BATTERY Assuming 6C discharge rate
% ACFT.Propulsion.Hydrogen.electricBladeN   = 2;

% Electric retrofit
ACFT.Propulsion.Electric.maxElectricPower_kW   = ACFT.Propulsion.maxPower_kW;
% ACFT.Propulsion.Electric.electric_n              = 2;

% Compute the thrust coefficients for each of the retrofits
% ACFT = calculate_thrustCoefficients(ACFT);


% Define efficiencies and fuel consumptions
% General values
pistonSpecificFuelCons_kgkWh = 0.27;
pistonMaxSpecificFuelCons_kgkWh =0.337;
dieselSpecificFuelCons_kgkWh = 0.23;
dieselMaxSpecificFuelCons_kgkWh = 0.25;
% Baseline
ACFT.Propulsion.Baseline.rhoFuel_kgm3            = rhoAVGAS_kgm3;
ACFT.Propulsion.Baseline.ICE_specificFuelCons_kgkWh = pistonSpecificFuelCons_kgkWh;   % Cruise Specific fuel consumption
ACFT.Propulsion.Baseline.ICE_maxSpecificFuelCons_kgkWh = pistonMaxSpecificFuelCons_kgkWh;% Taxi, take-off and climb specific fuel consumption
ACFT.Propulsion.Baseline.propEff                 = 0.8;
% Parallel
ACFT.Propulsion.Parallel.rhoFuel_kgm3            = rhoAVGAS_kgm3;
ACFT.Propulsion.Parallel.ICE_specificFuelCons_kgkWh = pistonSpecificFuelCons_kgkWh;  % Cruise Specific fuel consumption
ACFT.Propulsion.Parallel.ICE_maxSpecificFuelCons_kgkWh = pistonMaxSpecificFuelCons_kgkWh;
ACFT.Propulsion.Parallel.propEff           = 0.8;

ACFT.Propulsion.Parallel.electricMotorEff = 0.95;
ACFT.Propulsion.Parallel.batteryEff       = 0.9;

% Series 
ACFT.Propulsion.Series.rhoFuel_kgm3            = rhoAVGAS_kgm3;
ACFT.Propulsion.Series.ICE_specificFuelCons_kgkWh =pistonSpecificFuelCons_kgkWh;  % Cruise Specific fuel consumption
ACFT.Propulsion.Series.ICE_maxSpecificFuelCons_kgkWh= pistonMaxSpecificFuelCons_kgkWh;
ACFT.Propulsion.Series.propEff           = 0.8;

ACFT.Propulsion.Series.generatorEff       = 0.9;
ACFT.Propulsion.Series.electricMotorEff   = 0.95;
ACFT.Propulsion.Series.batteryEff         = 0.9;
ACFT.Propulsion.Series.chargeEff          = 0.9;

% Hydrogen (Fuel Cell)
ACFT.Propulsion.Hydrogen.rhoFuel_kgm3      = rhoLH2_kgm3;
ACFT.Propulsion.Hydrogen.fuelCellEff       = 0.55;
ACFT.Propulsion.Hydrogen.FC_specificFuelCons_kgkWh= ...
    1/ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg/...
    ACFT.Propulsion.Hydrogen.fuelCellEff;

ACFT.Propulsion.Hydrogen.electricMotorEff   = 0.95;
ACFT.Propulsion.Hydrogen.batteryEff         = 0.9;
ACFT.Propulsion.Hydrogen.propEff = 0.8;
% 
% % Electric
ACFT.Propulsion.Electric.propEff            = 0.8;
ACFT.Propulsion.Electric.generatorEff       = 0.9;
ACFT.Propulsion.Electric.electricMotorEff   = 0.95;
ACFT.Propulsion.Electric.batteryEff         = 0.9;
ACFT.Propulsion.Electric.chargeEff          = 0.9;

% Define Energy storage for all the retrofits (Important for mission analysis and aircraft design)
% Baseline
ACFT.EnergyStorage.Baseline.maxFuelVolume_m3 = ACFT.maxFuelVolume_m3;
ACFT.EnergyStorage.Baseline.maxFuelWeight_kg = ACFT.EnergyStorage.Baseline.maxFuelVolume_m3...
    *ACFT.Propulsion.Baseline.rhoFuel_kgm3;
% Parallel
ACFT.EnergyStorage.Parallel.maxFuelVolume_m3 = ACFT.maxFuelVolume_m3;
ACFT.EnergyStorage.Parallel.maxFuelWeight_kg = ACFT.EnergyStorage.Parallel.maxFuelVolume_m3...
    *ACFT.Propulsion.Parallel.rhoFuel_kgm3;
% Series
ACFT.EnergyStorage.Series.maxFuelVolume_m3 = ACFT.maxFuelVolume_m3;
ACFT.EnergyStorage.Series.maxFuelWeight_kg = ACFT.EnergyStorage.Series.maxFuelVolume_m3...
    *ACFT.Propulsion.Series.rhoFuel_kgm3;
% Hydrogen
ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg = 1/0.03;
ACFT.EnergyStorage.Hydrogen.tankNumber      = 2;
ACFT.EnergyStorage.Hydrogen.slendernessRatio= 4;

%% Weight Analysis and Sizing of energy storage system
ACFT.Weight = calculate_baselineEmptyWeight(ACFT);
[ACFT] = calculate_weightsSizing(ACFT);
% After weight analysis compute the maximum amount of hydrogen fuel
weight_N = ACFT.Weight.MTOW_kg*9.81;
% Could calculate the hydrogen tank dimensions and compute the drag

%% Mission Performance
% Need to calculate paylaod-range
PayloadRange = calculate_payloadRange(ACFT,Mission,MissionPerformance);

MissionPerformance.Baseline = calculate_missionPerformanceBaseline(Mission,MissionPerformance,ACFT,ACFT.Weight.MTOW_kg*9.81);
hybridCruise = 1;           % Select 1 or 0 for hybrid cruise or not
MissionPerformance.Parallel = calculate_missionPerformanceParallel(Mission,MissionPerformance,ACFT,ACFT.Weight.MTOW_kg*9.81,hybridCruise);
rechargeSeries = 0;
MissionPerformance.Series   = calculate_missionPerformanceSeries(Mission,MissionPerformance,ACFT,ACFT.Weight.MTOW_kg*9.81,hybridCruise,rechargeSeries);
MissionPerformance.Hydrogen = calculate_missionPerformanceHydrogen(Mission,MissionPerformance,ACFT,ACFT.Weight.MTOW_kg*9.81,hybridCruise);
MissionPerformance.Electric = calculate_missionPerformanceElectric(Mission,MissionPerformance,ACFT,ACFT.Weight.MTOW_kg*9.81);

%% Cost Calculation
ACFT.Cost = calculate_cost(ACFT,Mission,MissionPerformance);

%% Plot and Tabulate data
print_missionPerformanceBaseline(Mission,MissionPerformance);
print_missionPerformanceParallel(Mission,MissionPerformance);
print_missionPerformanceSeries(Mission,MissionPerformance);
print_missionPerformanceHydrogen(Mission,MissionPerformance);
print_missionPerformanceElectric(Mission,MissionPerformance);
print_missionSummary(MissionPerformance);

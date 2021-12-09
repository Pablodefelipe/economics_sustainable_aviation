function [ACFT] = calculate_weightsSizing(ACFT)
% Function:
%   calculate_weights
%  
% Description: 
%   Calculates the weights for the baseline, parallel, series and hydrogen
%   retrofits. Computes the weight of the energy storage system available
%
% Input:
%   ACFT        - Structure with aircarft data
% Output:
%   ACFT.Weight - Structure with weights of all aircarft
% Copy ACFT

% Unwrap ACFT structures
maxPower_kW  = ACFT.Propulsion.maxPower_kW;
hRatioEnergy = ACFT.Propulsion.hybridisationRatioEnergy;   
Weight       = ACFT.Weight;
EnergyStorage= ACFT.EnergyStorage;
avgasSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
% Unwrap Structure
MTOW_kg                         = Weight.MTOW_kg;
baselineEmpty_kg                = Weight.baselineEmpty_kg;
crew_kg                         = Weight.crew_kg;
passenger_kg                    = Weight.passenger_kg;
passenger_n                     = Weight.passenger_n;
fuelTankSystem_kg               = Weight.fuelTankSystem_kg;

maxFuelVolume_m3 = ACFT.EnergyStorage.Baseline.maxFuelVolume_m3;
engineNum        = ACFT.engineNum;

electricMotorSpecificPower_kWkg = Weight.electricMotorSpecificPower_kWkg;
batterySpecificEnergy_kWhkg     = Weight.batterySpecificEnergy_kWhkg;
icePistonSpecificPower_kWkg     = Weight.icePistonSpecificPower_kWkg;
iceDieselSpecificPower_kWkg     = Weight.iceDieselSpecificPower_kWkg;
iceTurbineSpecificPower_kWkg    = Weight.iceTurbineSpecificPower_kWkg;
gearboxSpecificPower_kWkg       = Weight.gearboxSpecificPower_kWkg;
fuelCellSpecificPower_kWkg      = Weight.fuelCellSpecificPower_kWkg;
electricSysSpecificPower_kWkg   = Weight.electricSysSpecificPower_kWkg; 
hydrogenFuelTankAreaDensity_kgm2= Weight.hydrogenFuelTankAreaDensity_kgm2;


% Calculate payload weigth due to passengers
payload_kg = passenger_kg*passenger_n;
%------------------------------Baseline------------------------------------
% Calculate the weight chapters of the propulsion system
Baseline.fuelEngine_kg      = maxPower_kW/icePistonSpecificPower_kWkg;
Baseline.electricMotor_kg   = 0; 
Baseline.battery_kg         = 0;
Baseline.electricSystems_kg = 0;
Baseline.fuelTank_kg        = calculate_fuelTankWeight(maxFuelVolume_m3,engineNum);

Baseline.emptyWeight_kg     = baselineEmpty_kg+Baseline.fuelEngine_kg+Baseline.electricMotor_kg+...
    Baseline.battery_kg  + Baseline.electricSystems_kg + Baseline.fuelTank_kg;
Baseline.totalPropulsion_kg = Baseline.fuelEngine_kg+Baseline.electricMotor_kg+...
    Baseline.battery_kg  + Baseline.electricSystems_kg + Baseline.fuelTank_kg;
Baseline.operationalEmptyWeight_kg = Baseline.emptyWeight_kg+crew_kg;
Baseline.zeroFuelWeight_kg  = Baseline.operationalEmptyWeight_kg + passenger_n*passenger_kg; 
Baseline.fuelWeight_kg      = MTOW_kg-Baseline.zeroFuelWeight_kg;
EnergyStorage.Baseline.fuelWeight_kg = Baseline.fuelWeight_kg;
EnergyStorage.Baseline.fuelEnergy_kWh= Baseline.fuelWeight_kg*...
    avgasSpecificEnergy_kWhkg;
EnergyStorage.Baseline.batteryEnergy_kWh = Baseline.battery_kg*...
    batterySpecificEnergy_kWhkg;
%------------------------------Parallel------------------------------------
maxICE_Power_kW             = ACFT.Propulsion.Parallel.maxICE_Power_kW;
maxElectricPower_kW         = ACFT.Propulsion.Parallel.maxElectricPower_kW;

gearboxWeight_kg            = maxPower_kW/gearboxSpecificPower_kWkg;
Parallel.fuelEngine_kg      = maxICE_Power_kW/icePistonSpecificPower_kWkg+gearboxWeight_kg;
Parallel.electricMotor_kg   = maxElectricPower_kW/electricMotorSpecificPower_kWkg; 

Parallel.electricSystems_kg = maxElectricPower_kW/...
    electricSysSpecificPower_kWkg;
Parallel.fuelTank_kg        = fuelTankSystem_kg;

Parallel.totalPropulsion_kg = Parallel.fuelEngine_kg+Parallel.electricMotor_kg+...
    + Parallel.electricSystems_kg;
Parallel.energyStorage_kg   = MTOW_kg - baselineEmpty_kg - payload_kg - ...
    crew_kg - Parallel.totalPropulsion_kg;
batteryEnergy_kWh = (Parallel.energyStorage_kg-fuelTankSystem_kg)/...    
    (1/batterySpecificEnergy_kWhkg+(1/hRatioEnergy-1)/avgasSpecificEnergy_kWhkg);
Parallel.battery_kg         = batteryEnergy_kWh/batterySpecificEnergy_kWhkg;
Parallel.emptyWeight_kg     = baselineEmpty_kg+Parallel.fuelEngine_kg+Parallel.electricMotor_kg+...
    Parallel.battery_kg  + Parallel.electricSystems_kg + Parallel.fuelTank_kg;
Parallel.operationalEmptyWeight_kg = Parallel.emptyWeight_kg+crew_kg;
Parallel.zeroFuelWeight_kg  = Parallel.operationalEmptyWeight_kg+payload_kg; 
Parallel.fuelWeight_kg      = MTOW_kg-Parallel.zeroFuelWeight_kg;
EnergyStorage.Parallel.fuelWeight_kg = Parallel.fuelWeight_kg;
EnergyStorage.Parallel.fuelEnergy_kWh= Parallel.fuelWeight_kg*...
    avgasSpecificEnergy_kWhkg;
EnergyStorage.Parallel.batteryEnergy_kWh = batteryEnergy_kWh;

%------------------------------Series--------------------------------------
maxICE_Power_kW             = ACFT.Propulsion.Series.maxICE_Power_kW;
maxElectricPower_kW         = ACFT.Propulsion.Series.maxElectricPower_kW;

Series.fuelEngine_kg      = maxICE_Power_kW*(1/icePistonSpecificPower_kWkg+1/electricMotorSpecificPower_kWkg); % Including electric generator
Series.electricMotor_kg   = maxElectricPower_kW/electricMotorSpecificPower_kWkg; 

Series.electricSystems_kg = maxElectricPower_kW/...
    electricSysSpecificPower_kWkg;
Series.totalPropulsion_kg = Series.fuelEngine_kg+Series.electricMotor_kg+...
    + Series.electricSystems_kg;
Series.fuelTank_kg        = fuelTankSystem_kg;
Series.energyStorage_kg   = MTOW_kg - baselineEmpty_kg - payload_kg - ...
    crew_kg - Series.totalPropulsion_kg;
batteryEnergy_kWh = (Series.energyStorage_kg-fuelTankSystem_kg)/...    
    (1/batterySpecificEnergy_kWhkg+(1/hRatioEnergy-1)/avgasSpecificEnergy_kWhkg);
Series.battery_kg         = batteryEnergy_kWh/batterySpecificEnergy_kWhkg;
Series.emptyWeight_kg     = baselineEmpty_kg+Series.fuelEngine_kg+Series.electricMotor_kg+...
    Series.battery_kg  + Series.electricSystems_kg + Series.fuelTank_kg;
Series.operationalEmptyWeight_kg = Series.emptyWeight_kg+crew_kg;
Series.zeroFuelWeight_kg  = Series.operationalEmptyWeight_kg + payload_kg; 
Series.fuelWeight_kg      = MTOW_kg-Series.zeroFuelWeight_kg;
EnergyStorage.Series.fuelWeight_kg = Series.fuelWeight_kg;
EnergyStorage.Series.fuelEnergy_kWh= Series.fuelWeight_kg*...
    avgasSpecificEnergy_kWhkg;
EnergyStorage.Series.batteryEnergy_kWh = batteryEnergy_kWh;


%------------------------------Hydrogen------------------------------------
rhoFuel_kgm3                = ACFT.Propulsion.Hydrogen.rhoFuel_kgm3;
maxFC_Power_kW              = ACFT.Propulsion.Hydrogen.maxFC_Power_kW;
maxElectricPower_kW         = ACFT.Propulsion.Hydrogen.maxElectricPower_kW;
LH2SpecificEnergy_kWhkg     = ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg;
tankNumber                  = ACFT.EnergyStorage.Hydrogen.tankNumber; 


Hydrogen.fuelEngine_kg      = maxICE_Power_kW/fuelCellSpecificPower_kWkg;
Hydrogen.electricMotor_kg   = maxElectricPower_kW/electricMotorSpecificPower_kWkg; 
Hydrogen.electricSystems_kg = maxElectricPower_kW/...
    electricSysSpecificPower_kWkg;
Hydrogen.totalPropulsion_kg = Hydrogen.fuelEngine_kg+Hydrogen.electricMotor_kg+...
    Hydrogen.electricSystems_kg;
Hydrogen.energyStorage_kg   = MTOW_kg - baselineEmpty_kg - payload_kg - ...
    crew_kg - Hydrogen.totalPropulsion_kg;
% Before calculating the weight need to find maximum possible
[batteryEnergy_kWh,fuelTank_kg]  = calculate_hydrogenEnergyStorage(ACFT,Hydrogen);
Hydrogen.battery_kg         = batteryEnergy_kWh/batterySpecificEnergy_kWhkg;
Hydrogen.fuelTank_kg        = fuelTank_kg;
Hydrogen.emptyWeight_kg     = baselineEmpty_kg+Hydrogen.fuelEngine_kg+...
    Hydrogen.electricMotor_kg+Hydrogen.battery_kg+...
    Hydrogen.electricSystems_kg + Hydrogen.fuelTank_kg;
Hydrogen.operationalEmptyWeight_kg = Hydrogen.emptyWeight_kg+crew_kg;
Hydrogen.zeroFuelWeight_kg  = Hydrogen.operationalEmptyWeight_kg+payload_kg; 
Hydrogen.fuelWeight_kg      = MTOW_kg-Hydrogen.zeroFuelWeight_kg;
EnergyStorage.Hydrogen.fuelEnergy_kWh= Hydrogen.fuelWeight_kg*...
    LH2SpecificEnergy_kWhkg;
EnergyStorage.Hydrogen.fuelWeight_kg = Hydrogen.fuelWeight_kg;;
EnergyStorage.Hydrogen.maxFuelWeight_kg = Hydrogen.fuelWeight_kg;
EnergyStorage.Hydrogen.batteryEnergy_kWh = batteryEnergy_kWh;
%------------------------------Electric------------------------------------
maxElectricPower_kW         = ACFT.Propulsion.Electric.maxElectricPower_kW;
batterySpecificEnergy_kWhkg = ACFT.EnergyStorage.Electric.batterySpecificEnergy_kWhkg;
Electric.electricMotor_kg   = maxElectricPower_kW/electricMotorSpecificPower_kWkg; 
Electric.electricSystems_kg = maxElectricPower_kW/...
    electricSysSpecificPower_kWkg;
% Calculate the available weight for the battery
Electric.battery_kg         = MTOW_kg-baselineEmpty_kg-crew_kg-...
    Electric.electricMotor_kg-Electric.electricSystems_kg-payload_kg;
Electric.emptyWeight_kg = baselineEmpty_kg+Electric.electricMotor_kg+...
    Electric.electricSystems_kg+Electric.battery_kg;
Electric.operationalEmptyWeight_kg = Electric.emptyWeight_kg+crew_kg;
EnergyStorage.Electric.batteryEnergy_kWh= Electric.battery_kg*...
    batterySpecificEnergy_kWhkg;

% Wrap up all the structures
Weight.Baseline     = Baseline;
Weight.Parallel     = Parallel;
Weight.Series       = Series;
Weight.Hydrogen     = Hydrogen;
Weight.Electric     = Electric;

ACFT.Weight         = Weight;
ACFT.EnergyStorage  = EnergyStorage;



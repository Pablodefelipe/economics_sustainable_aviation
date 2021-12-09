function writeExcel_weight(ACFT)
% Function:
%   writeExcel_weight(ACFT)
%
% Description: 
%   Exports the weight summary to excel
%
% Input:
%   ACFT                - Aircraft data

% Output:

% Unwrap structure
Weight      = ACFT.Weight;
electricMotorSpecificPower_kWkg = Weight.electricMotorSpecificPower_kWkg;
batterySpecificEnergy_kWhkg     = Weight.batterySpecificEnergy_kWhkg;
icePistonSpecificPower_kWkg     = Weight.icePistonSpecificPower_kWkg;
gearboxSpecificPower_kWkg       = Weight.gearboxSpecificPower_kWkg;
fuelCellSpecificPower_kWkg      = Weight.fuelCellSpecificPower_kWkg;
electricSysSpecificPower_kWkg   = Weight.electricSysSpecificPower_kWkg;
hydrogenFuelTankAreaDensity_kgm2= Weight.hydrogenFuelTankAreaDensity_kgm2;

Baseline    = Weight.Baseline;
Parallel    = Weight.Parallel;
Series      = Weight.Series;
Hydrogen    = Weight.Hydrogen;
Electric    = Weight.Electric;

% first table of power densities
filename    = 'exampleData.xlsx';
Component   = ["Electric Motor";"Battery";"ICE";"Gearbox";"Fuel Cell";...
    "Electric Systems";"LH2 Fuel Tank"];
Value   = [electricMotorSpecificPower_kWkg;batterySpecificEnergy_kWhkg;...
    icePistonSpecificPower_kWkg;gearboxSpecificPower_kWkg;...
    fuelCellSpecificPower_kWkg;electricSysSpecificPower_kWkg;...
    hydrogenFuelTankAreaDensity_kgm2];
Units =  ["kW/kg";"kWh/kg";"kW/kg";"kW/kg";"kW/kg";"kW/kg";"kg/m^2"];

densities = table(Component,Value,Units);

writetable(densities,filename,'Sheet',1,'Range','A30')

Mass        = ["fuel engine";"electric Motor";"Battery";"Elec. systems";"Fuel Tank"];
Baseline    = [Baseline.fuelEngine_kg;Baseline.electricMotor_kg;Baseline.battery_kg;Baseline.electricSystems_kg;Baseline.fuelTank_kg];
Parallel    = [Parallel.fuelEngine_kg;Parallel.electricMotor_kg;Parallel.battery_kg;Parallel.electricSystems_kg;Parallel.fuelTank_kg];
Series      = [Series.fuelEngine_kg;Series.electricMotor_kg;Series.battery_kg;Series.electricSystems_kg;Series.fuelTank_kg];
Hydrogen    = [Hydrogen.fuelEngine_kg;Hydrogen.electricMotor_kg;Hydrogen.battery_kg;Hydrogen.electricSystems_kg;Hydrogen.fuelTank_kg];
Electric    = [0;Electric.electricMotor_kg;Electric.battery_kg;Electric.electricSystems_kg;0];
weight = table(Mass,Baseline,Parallel,Series,Hydrogen,Electric);

writetable(weight,filename,'Sheet',1,'Range','A39')
end

function writeExcel_propulsion(ACFT)
% Function:
%   writeExcel_propulsion(ACFT)
%
% Description: 
%   Exports the propulsion system summary to excel
%
% Input:
%   ACFT                - Aircraft data

% Output:

% Unwrap structure
Propulsion  = ACFT.Propulsion;
Baseline    = Propulsion.Baseline;
Parallel    = Propulsion.Parallel;
Series      = Propulsion.Series;
Hydrogen    = Propulsion.Hydrogen;
Electric    = Propulsion.Electric;

a = 'maxICE_Power_kW';
b = 'maxFC_Power_kW';
c = 'maxElectricPower_kW';

filename    = 'exampleData.xlsx';
Retrofit    = ["Baseline";"Parallel";"Series";"Hydrogen";"Electric"];
FuelPower   = [Baseline.(a);Parallel.(a);Series.(a);Hydrogen.(b);0];
ElecPower   = [0;Parallel.(c);Series.(c);Hydrogen.(c);Electric.(c)];

EnergyStorage  = ACFT.EnergyStorage;
Baseline    = EnergyStorage.Baseline;
Parallel    = EnergyStorage.Parallel;
Series      = EnergyStorage.Series;
Hydrogen    = EnergyStorage.Hydrogen;
Electric    = EnergyStorage.Electric;
d = 'fuelEnergy_kWh';
e = 'batteryEnergy_kWh';
FuelEnergy  = [Baseline.(d);Parallel.(d);Series.(d);Hydrogen.(d);0];
BattEnergy  = [Baseline.(e);Parallel.(e);Series.(e);Hydrogen.(e);Electric.(e)(1)];
propulsion  = table(Retrofit,FuelPower,ElecPower,FuelEnergy,BattEnergy);

writetable(propulsion,filename,'Sheet',1,'Range','A23')
end

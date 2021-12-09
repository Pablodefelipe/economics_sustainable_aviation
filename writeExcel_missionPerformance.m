function writeExcel_missionPerformance(ACFT,Mission,MissionPerformance)
% Function:
%   print_missionPerformaceBaseline
%
% Description: 
%   Exports the average power and fuel consumption for a particular mission
%   profile to excel for all the retrofits
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:
filename = 'exampleData.xlsx';

%% Baseline
% Unwrap power
Power = MissionPerformance.Baseline.ICE_power_kW;
Fuel  = MissionPerformance.Baseline.fuelCons_USgal;
Fuel.reserve = MissionPerformance.Baseline.fuelConsReserve_USgal;


Segment  = ["Avg.Power(kW)";"Fuel (USgal)"];
Taxi     = [Power.taxi;Fuel.taxi];
TakeOff  = [Power.takeOff;Fuel.takeOff];
Climb    = [Power.climb;Fuel.climb];
Cruise   = [Power.cruise;Fuel.cruise];
Descent  = [Power.descent;Fuel.descent];
Reserve  = [Power.reserve;Fuel.reserve];

missionPerformance = table(Segment,Taxi,TakeOff,Climb,Cruise,Descent,Reserve);

writetable(missionPerformance,filename,'Sheet',1,'Range','A47')
%% Parallel
% Unwrap power
Power = MissionPerformance.Parallel.ICE_power_kW;
Fuel  = MissionPerformance.Parallel.fuelCons_USgal;
Fuel.reserve = MissionPerformance.Parallel.fuelConsReserve_USgal;
ElecPower    = MissionPerformance.Parallel.electricPower_kW;
ElecEnergy   = MissionPerformance.Parallel.electricEnergy_kWh;
ElecEnergy.reserve = MissionPerformance.Parallel.electricEnergyReserve_kWh;
BattEnergy   = MissionPerformance.Parallel.batteryEnergy_kWh;

Segment  = ["Avg.Power(kW)";"Fuel (USgal)";"Elec.Power(kW)";"ElecEnergy(kWh)";"Batt.State(kWh)"];
Taxi     = [Power.taxi;Fuel.taxi;ElecPower.taxi;ElecEnergy.taxi;BattEnergy.taxi];
TakeOff  = [Power.takeOff;Fuel.takeOff;ElecPower.takeOff;ElecEnergy.takeOff;BattEnergy.takeOff];
Climb    = [Power.climb;Fuel.climb;ElecPower.climb;ElecEnergy.climb;BattEnergy.climb];
Cruise   = [Power.cruise;Fuel.cruise;ElecPower.cruise;ElecEnergy.cruise;BattEnergy.cruise];
Descent  = [Power.descent;Fuel.descent;ElecPower.descent;ElecEnergy.descent;BattEnergy.descent];
Reserve  = [Power.reserve;Fuel.reserve;ElecPower.reserve;ElecEnergy.reserve;BattEnergy.reserveClimb];

missionPerformance = table(Segment,Taxi,TakeOff,Climb,Cruise,Descent,Reserve);

writetable(missionPerformance,filename,'Sheet',1,'Range','A52')

%% Series
% Unwrap power
Power = MissionPerformance.Series.ICE_power_kW;
Fuel  = MissionPerformance.Series.fuelCons_USgal;
Fuel.reserve = MissionPerformance.Series.fuelConsReserve_USgal;
ElecPower    = MissionPerformance.Series.ElectricPower_kW;
GenEnergy   = MissionPerformance.Series.GenEnergy_kWh;
GenEnergy.reserve = MissionPerformance.Series.genEnergyReserve_kWh;
BattEnergy   = MissionPerformance.Series.BattEnergy_kWh;
BattState    = MissionPerformance.Series.batteryEnergy_kWh;

Segment  = ["Avg.E.Power(kW)";"ICE. Power(kWh)";"Fuel (USgal)";"Gen.Energy(kWh)";"BattEnergy(kWh)";"Batt.State(kWh)"];
Taxi     = [ElecPower.taxi;Power.taxi;Fuel.taxi;GenEnergy.taxi;BattEnergy.taxi;BattState.taxi];
TakeOff  = [ElecPower.takeOff;Power.takeOff;Fuel.takeOff;GenEnergy.takeOff;BattEnergy.takeOff;BattState.takeOff];
Climb    = [ElecPower.climb;Power.climb;Fuel.climb;GenEnergy.climb;BattEnergy.climb;BattState.climb];
Cruise   = [ElecPower.cruise;Power.cruise;Fuel.cruise;GenEnergy.cruise;BattEnergy.cruise;BattState.cruise];
Descent  = [ElecPower.descent;Power.descent;Fuel.descent;GenEnergy.descent;BattEnergy.descent;BattState.descent];
Reserve  = [ElecPower.reserve;Power.reserve;Fuel.reserve;GenEnergy.reserve;BattEnergy.reserveClimb;BattState.reserveClimb];

missionPerformance = table(Segment,Taxi,TakeOff,Climb,Cruise,Descent,Reserve);

writetable(missionPerformance,filename,'Sheet',1,'Range','A60')

%% Hydrogen
% Unwrap power
ElecPower    = MissionPerformance.Hydrogen.ElectricPower_kW;
Fuel         = MissionPerformance.Hydrogen.FuelCons_USgal;
Fuel.reserve = MissionPerformance.Hydrogen.fuelConsReserve_USgal;
FCPower      = MissionPerformance.Hydrogen.FC_power_kW; 
BattPower    = MissionPerformance.Hydrogen.BattPower_kW;
BattEnergy   = MissionPerformance.Hydrogen.BattEnergy_kWh;
BattState    = MissionPerformance.Hydrogen.batteryEnergy_kWh;

Segment  = ["Avg.E.Power(kW)";"FC. Power(kW)";"Batt.Power(kW)";"Fuel (USgal)";"BattEnergy(kWh)";"Batt.State(kWh)"];
Taxi     = [ElecPower.taxi;FCPower.taxi;BattPower.taxi;Fuel.taxi;BattEnergy.taxi;BattState.taxi];
TakeOff  = [ElecPower.takeOff;FCPower.takeOff;BattPower.takeOff;Fuel.takeOff;BattEnergy.takeOff;BattState.takeOff];
Climb    = [ElecPower.climb;FCPower.climb;BattPower.climb;Fuel.climb;BattEnergy.climb;BattState.climb];
Cruise   = [ElecPower.cruise;FCPower.cruise;BattPower.cruise;Fuel.cruise;BattEnergy.cruise;BattState.cruise];
Descent  = [ElecPower.descent;FCPower.descent;BattPower.descent;Fuel.descent;BattEnergy.descent;BattState.descent];
Reserve  = [ElecPower.reserve;FCPower.reserve;BattPower.reserve;Fuel.reserve;BattEnergy.reserveClimb;BattState.reserveClimb];

missionPerformance = table(Segment,Taxi,TakeOff,Climb,Cruise,Descent,Reserve);

writetable(missionPerformance,filename,'Sheet',1,'Range','A69')

%% Electric
% Unwrap power
ElecPower    = MissionPerformance.Electric.ElectricPower_kW;
BattEnergy   = MissionPerformance.Electric.BattEnergy_kWh;
BattState    = MissionPerformance.Electric.batteryEnergy_kWh;

Segment  = ["Avg.E.Power(kW)";"BattEnergy(kWh)";"Batt.State(kWh)"];
Taxi     = [ElecPower.taxi;BattEnergy.taxi;BattState.taxi(1)];
TakeOff  = [ElecPower.takeOff;BattEnergy.takeOff;BattState.takeOff(1)];
Climb    = [ElecPower.climb;BattEnergy.climb;BattState.climb(1)];
Cruise   = [ElecPower.cruise;BattEnergy.cruise;BattState.cruise(1)];
Descent  = [ElecPower.descent;BattEnergy.descent;BattState.descent(1)];
Reserve  = [ElecPower.reserve;BattEnergy.reserveClimb;BattState.reserveClimb(1)];

missionPerformance = table(Segment,Taxi,TakeOff,Climb,Cruise,Descent,Reserve);

writetable(missionPerformance,filename,'Sheet',1,'Range','A78')

%% Compare the fuel consumption of the main retrofits against the baseline
% Use litres

fuelBaseline_l = MissionPerformance.Baseline.fuelConsMission_USgal*3.785;
fuelParallel_l = MissionPerformance.Parallel.fuelConsMission_USgal*3.785;
fuelSeries_l   = MissionPerformance.Series.fuelConsMission_USgal*3.785;
percentage = (fuelBaseline_l-[fuelBaseline_l,fuelParallel_l,fuelSeries_l])...
    /fuelBaseline_l*100;
Retrofit = ["Fuel (l)";"Decrease (%)"];
Baseline = [fuelBaseline_l;percentage(1)];
Parallel = [fuelParallel_l;percentage(2)];
Series   = [fuelSeries_l;percentage(3)];

fuelConsumption = table(Retrofit,Baseline,Parallel,Series);

writetable(fuelConsumption,filename,'Sheet',1,'Range','A84')

%% Compare the efficiency of the 
% Get specific energy of fuels
avgasSpecificEnergy_kWhkg  = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
lh2SpecificEnergy_kWhkg    = ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg;
% Get fuel and battery energy for all the retrofits
fuel_kg          = MissionPerformance.Baseline.fuelConsMission_kg;
fuelBaseline_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
baseline         = [fuelBaseline_kWh,0,0,0];

% Parallel
fuel_kg  = MissionPerformance.Parallel.fuelConsMission_kg;
fuelParallel_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
battParallel_kWh = MissionPerformance.Parallel.electricEnergyMission_kWh;
parallel = [fuelParallel_kWh,battParallel_kWh,0,0];
% Series
fuel_kg  = MissionPerformance.Series.fuelConsMission_kg;
fuelSeries_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
battSeries_kWh = MissionPerformance.Series.battEnergyMission_kWh;
series = [fuelSeries_kWh,battSeries_kWh,0,0];
% Hydorgen
fuel_kg  = MissionPerformance.Hydrogen.fuelConsMission_kg;
fuelHydrogen_kWh = fuel_kg*lh2SpecificEnergy_kWhkg;
battHydrogen_kWh = MissionPerformance.Hydrogen.battEnergyMission_kWh;
hydrogen = [0,battHydrogen_kWh,fuelHydrogen_kWh,0];
% Electric
battElectric_kWh = MissionPerformance.Electric.battEnergyMission_kWh;
electric = [0,battElectric_kWh,0,0];
% Theoretical
theoreticalEnergy = MissionPerformance.EnergyEstimate_kWh.totalMission;
theoretical = [0,0,0,theoreticalEnergy];

effBaseline = theoreticalEnergy/fuelBaseline_kWh*100;
effParallel = theoreticalEnergy/(fuelParallel_kWh+battParallel_kWh)*100;
effSeries = theoreticalEnergy/(fuelSeries_kWh+battSeries_kWh)*100;
effHydrogen = theoreticalEnergy/(fuelHydrogen_kWh+battHydrogen_kWh)*100;
effElectric = theoreticalEnergy/(battElectric_kWh)*100;

Retrofit = ["Baseline";"Parallel";"Series";"Hydrogen";"Electric"];
Efficiency = [effBaseline;effParallel;effSeries;effHydrogen;effElectric];

MissionEfficiency = table(Retrofit,Efficiency);

writetable(MissionEfficiency,filename,'Sheet',1,'Range','A89')

end


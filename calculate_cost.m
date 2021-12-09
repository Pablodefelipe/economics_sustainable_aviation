function Cost = calculate_cost(ACFT,Mission,MissionPerformance)
% Function:
%   calculate_missionPerformace
%
% Description: 
%   Calculate distance and time spent in each segment of the mission
%
% Input:
%   ACFT                - Aircraft parameters
%   Mission             - Mission parameters
%   MissionPerformance  - Mission Perfromance data
% Output:
%   Cost                - Structure with costs of fuel and battery for each of
%                           the retrofits

% Unwrap structures
avgasSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
lh2SpecificEnergy_kWhkg    = ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg;

avgas_usCentskWh        = ACFT.Cost.avgas_usCentskWh;
hydrogen_usCentskWh     = ACFT.Cost.hydrogen_usCentskWh;
electricity_usCentskWh  = ACFT.Cost.electricity_usCentskWh;
% Get fuel and battery energy for all the retrofits
fuel_kg  = MissionPerformance.Baseline.fuelConsMission_kg;
fuelBaseline_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
baseline = [fuelBaseline_kWh,0,0,0];
Baseline.fuel_cents = fuelBaseline_kWh*avgas_usCentskWh;
% Parallel
fuel_kg  = MissionPerformance.Parallel.fuelConsMission_kg;
fuelParallel_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
battParallel_kWh = MissionPerformance.Parallel.electricEnergyMission_kWh;
Parallel.fuel_cents = fuelParallel_kWh*avgas_usCentskWh;
Parallel.elec_cents = battParallel_kWh*electricity_usCentskWh;
Parallel.total_cents= Parallel.fuel_cents + Parallel.elec_cents;
parallel = [fuelParallel_kWh,battParallel_kWh,0,0];
% Series
fuel_kg  = MissionPerformance.Series.fuelConsMission_kg;
fuelSeries_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
battSeries_kWh = MissionPerformance.Series.battEnergyMission_kWh;
Series.fuel_cents = fuelSeries_kWh*avgas_usCentskWh;
Series.elec_cents = battSeries_kWh*electricity_usCentskWh;
Series.total_cents= Series.fuel_cents + Series.elec_cents;
series = [fuelSeries_kWh,battSeries_kWh,0,0];
% Hydorgen
fuel_kg  = MissionPerformance.Hydrogen.fuelConsMission_kg;
fuelHydrogen_kWh = fuel_kg*lh2SpecificEnergy_kWhkg;
battHydrogen_kWh = MissionPerformance.Hydrogen.battEnergyMission_kWh;
Hydrogen.fuel_cents = fuelHydrogen_kWh*hydrogen_usCentskWh;
Hydrogen.elec_cents = battHydrogen_kWh*electricity_usCentskWh;
Hydrogen.total_cents= Hydrogen.fuel_cents + Hydrogen.elec_cents;
hydrogen = [0,battHydrogen_kWh,fuelHydrogen_kWh,0];
% Electric
battElectric_kWh = MissionPerformance.Electric.battEnergyMission_kWh;
Electric.elec_cents = battElectric_kWh*electricity_usCentskWh;
electric = [0,battElectric_kWh,0,0];
%% Calulate available seat mile values
asm = Mission.distance.mission_m/1609*ACFT.Weight.passenger_n;
Baseline.total_casm = Baseline.fuel_cents/asm;
Parallel.total_casm = Parallel.total_cents/asm;
Series.total_casm = Series.total_cents/asm;
Hydrogen.total_casm = Hydrogen.total_cents/asm;
Electric.total_casm = Electric.elec_cents/asm;

% Wrap the structure
Cost.Baseline = Baseline;
Cost.Parallel = Parallel;
Cost.Series = Series;
Cost.Hydrogen = Hydrogen;
Cost.Electric = Electric;
end

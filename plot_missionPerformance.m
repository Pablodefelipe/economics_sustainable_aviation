function [] = plot_missionPerformance(ACFT,MissionPerformance)
% Function:
%   plot_missionPerformance
%
% Description: 
%   Plot the mission performance in a bar chart format
% Input:
%   ACFT                    - Aircraft parameters
%   MissionPerformance      - Performance of aircraft defined for all the
%                           missions       
% Output:
%   plot

% Get specific energy of fuels
avgasSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
lh2SpecificEnergy_kWhkg    = ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg;
% Get fuel and battery energy for all the retrofits
fuel_kg  = MissionPerformance.Baseline.fuelConsMission_kg;
fuelBaseline_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
baseline = [fuelBaseline_kWh,0,0,0];
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
theoretical = MissionPerformance.EnergyEstimate_kWh.totalMission;
theoretical = [0,0,0,theoretical];
figure()
X = categorical({'theoretical','baseline','parallel','series','hydrogen','electric'});
X = reordercats(X,{'theoretical','baseline','parallel','series','hydrogen','electric'});
Y = [theoretical;baseline;parallel;series;hydrogen;electric]; 
h = bar(X,Y,0.6,'stacked');
ylim([0 500])
ylabel('Energy (kWh)')
% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'AVGAS Energy','Battery Energy','LH2 Energy',...
    'Theoretical'}')
set(gca,'fontSize',16)
legend('Orientation','horizontal','NumColumns',2,'fontSize',14)
% Legend will show names for each color
h(1).FaceColor = [0.8500, 0.3250, 0.0980];
h(2).FaceColor = [0.4660, 0.6740, 0.1880];
h(3).FaceColor = [0.3010, 0.7450, 0.9330];
h(4).FaceColor = [128,128,128]/255;
end

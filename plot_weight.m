function [] = plot_weight(ACFT)
% Function:
%   plot_weight(ACFT)
%
% Description: 
%   Plot the weights of the retrofits in a bar chart format
% Input:
%   Weight              - Structure with baseline, parallel, series,
%                         weight data       
% Output:
%   plot

% Unwrap structures
Weight      = ACFT.Weight;
Baseline    = Weight.Baseline;
Parallel    = Weight.Parallel;
Series      = Weight.Series;
Hydrogen    = Weight.Hydrogen;
Electric    = Weight.Electric;

Baseline    = [Baseline.fuelEngine_kg;Baseline.electricMotor_kg;Baseline.battery_kg;Baseline.electricSystems_kg;Baseline.fuelTank_kg];
Parallel    = [Parallel.fuelEngine_kg;Parallel.electricMotor_kg;Parallel.battery_kg;Parallel.electricSystems_kg;Parallel.fuelTank_kg];
Series      = [Series.fuelEngine_kg;Series.electricMotor_kg;Series.battery_kg;Series.electricSystems_kg;Series.fuelTank_kg];
Hydrogen    = [Hydrogen.fuelEngine_kg;Hydrogen.electricMotor_kg;Hydrogen.battery_kg;Hydrogen.electricSystems_kg;Hydrogen.fuelTank_kg];
Electric    = [0;Electric.electricMotor_kg;Electric.battery_kg;Electric.electricSystems_kg;0];

figure()
X = categorical({'Baseline','Parallel','Series','Hydrogen','Electric'});
X = reordercats(X,{'Baseline','Parallel','Series','Hydrogen','Electric'});
Y = [Baseline';Parallel';Series';Hydrogen';Electric']; 
h = bar(X,Y,'stacked');
ylim([0 1200])
ylabel('Weight (kg)')
% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'Fuel Engine','Electric Motor','Battery',...
    'Electric Systems','Fuel Tanks'}')
set(gca,'fontSize',18)
% Legend will show names for each color
legend('Orientation','vertical','NumColumns',3,'fontSize',10);
h(1).FaceColor = [0.8500 0.3250 0.0980];
h(2).FaceColor = [0 0.4470 0.7410];
h(3).FaceColor = [0.4660 0.6740 0.1880];
h(4).FaceColor = [0.9290 0.6940 0.1250];
h(5).FaceColor = [0.3010 0.7450 0.9330];


end

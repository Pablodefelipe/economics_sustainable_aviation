function [] = plot_cost(ACFT)
% Function:
%   plot_cost(ACFT)
%
% Description: 
%   Plot the mission performance in a bar chart format
% Input:
%   ACFT                    - Aircraft parameters
%   MissionPerformance      - Performance of aircraft defined for all the
%                           missions       
% Output:
%   plot

% Unwrap the structure
% Wrap the structure
Cost = ACFT.Cost;
Baseline    = Cost.Baseline;
Parallel    = Cost.Parallel;
Series      = Cost.Series;
Hydrogen    = Cost.Hydrogen;
Electric    = Cost.Electric;

% Baseline
baseline = [Baseline.fuel_cents,0,0]/100;
% Parallel
parallel = [Parallel.fuel_cents,Parallel.elec_cents,0]/100;
% Series
series = [Series.fuel_cents,Series.elec_cents,0]/100;
% Hydorgen
hydrogen = [0,Hydrogen.elec_cents,Hydrogen.fuel_cents]/100;
% Electric
electric = [0,Electric.elec_cents,0]/100;

figure()
X = categorical({'baseline','parallel','series','hydrogen','electric'});
X = reordercats(X,{'baseline','parallel','series','hydrogen','electric'});
Y = [baseline;parallel;series;hydrogen;electric]; 
h = bar(X,Y,0.6,'stacked');
%ylim([0 500])
ylabel('Cost ($)')
% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'AVGAS','Electricity','Hydrogen'}')
set(gca,'fontSize',16)
legend('Orientation','horizontal','fontSize',14)
% Legend will show names for each color
h(1).FaceColor = [0.8500, 0.3250, 0.0980];
h(2).FaceColor = [0.4660, 0.6740, 0.1880];
h(3).FaceColor = [0.3010, 0.7450, 0.9330];
end

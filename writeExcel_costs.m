function writeExcel_costs(ACFT)
% Function:
%   writeExcel_costs
%
% Description: 
%   Exports the costs for a particular mission
%   profile to excel for all the retrofits
%
% Input:
%   ACFT        - Aircraft Parameters
% Output:
filename = 'exampleData.xlsx';

% Unwrap the structure
% Wrap the structure
Cost = ACFT.Cost;
baseline    = Cost.Baseline;
parallel    = Cost.Parallel;
series      = Cost.Series;
hydrogen    = Cost.Hydrogen;
electric    = Cost.Electric;

% Calculate % difference
diff = baseline.fuel_cents - [baseline.fuel_cents,parallel.total_cents,...
    series.total_cents,hydrogen.total_cents];
diff = diff/baseline.fuel_cents*100;

Retrofit       = ["CASM";"Diff. (%)"];
Baseline    = [baseline.total_casm;diff(1)];
Parallel    = [parallel.total_casm;diff(2)];
Series      = [series.total_casm;diff(3)];
Hydrogen    = [hydrogen.total_casm;diff(4)];
cost = table(Retrofit,Baseline,Parallel,Series,Hydrogen);

writetable(cost,filename,'Sheet',1,'Range','A100')

Retrofit       = ["Fuel ($)";"Elec. ($)";"Total ($)"];
Baseline    = [baseline.fuel_cents;0;baseline.fuel_cents]/100;
Parallel    = [parallel.fuel_cents;parallel.elec_cents;parallel.total_cents]/100;
Series      = [series.fuel_cents;series.elec_cents;series.total_cents]/100;
Hydrogen    = [hydrogen.fuel_cents;hydrogen.elec_cents;hydrogen.total_cents]/100;
Electric    = [0;electric.elec_cents;electric.elec_cents]/100;
costsDetailed = table(Retrofit,Baseline,Parallel,Series,Hydrogen,Electric);

writetable(costsDetailed,filename,'Sheet',1,'Range','A104')
end

function [] = plot_payloadRange(PayloadRange)
% Function:
%   plot_payloadRange
%
% Description: 
%   Plot the maximum range and payload for each of the retrofits
%
% Input:
%   PayloadRange        - Structure with baseline, parallel, series,
%                         hydrogen and electric payload range data       
% Output:
%   plot

% Unwrap structure
Baseline.payloadWeight_kg   = PayloadRange.Baseline.payloadWeight_kg;
Parallel.payloadWeight_kg   = PayloadRange.Parallel.payloadWeight_kg;
Series.payloadWeight_kg     = PayloadRange.Series.payloadWeight_kg;
Hydrogen.payloadWeight_kg   = PayloadRange.Hydrogen.payloadWeight_kg;
Electric.payloadWeight_kg   = PayloadRange.Electric.payloadWeight_kg;

Baseline.range_km   = PayloadRange.Baseline.range_m/1000;
Parallel.range_km   = PayloadRange.Parallel.range_m/1000;
Series.range_km     = PayloadRange.Series.range_m/1000;
Hydrogen.range_km   = PayloadRange.Hydrogen.range_m/1000;
Electric.range_km   = PayloadRange.Electric.range_m/1000;

figure();
hold on
pBaseline = plot(Baseline.range_km,Baseline.payloadWeight_kg,...
    'LineWidth',2,'Color',	[0.8500, 0.3250, 0.0980]);
pParallel= plot(Parallel.range_km,Parallel.payloadWeight_kg,'-',...
    'LineWidth',2,'Color', 	[0, 0.4470, 0.7410]);
pSeries = plot(Series.range_km,Series.payloadWeight_kg,...
    'LineWidth',2,'Color',[0.4940, 0.1840, 0.5560]);
pHydrogen =  plot(Hydrogen.range_km,Hydrogen.payloadWeight_kg,...
    'LineWidth',2,'Color',[0.3010, 0.7450, 0.9330]);
% Can be improved
pElectric1 = plot(Electric.range_km(1,:),Electric.payloadWeight_kg,...
    'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]);
pElectric2 = plot(Electric.range_km(2,:),Electric.payloadWeight_kg,...
    '--','LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]);
%legend([pBaseline,pParallel,pSeries,pHydrogen,pElectric1,pElectric2],...
    names = {'Baseline','Parallel','Series','Hydrogen','Electric (400Wh/kg)',...
    'Electric(600Wh/kg)'};
columnlegend(2, names, 'location','northEast');
box on
grid on
set(gca,'FontSize',18)
xlabel('Range (km)')
ylabel('Payload (kg)')
axis([0,2500,0,1000])
end
function PayloadRange = calculate_payloadRange(ACFT,Mission,MissionPerformance);
% Function:
%   calculate_payloadRange
%
% Description: 
%   Calculate the maximum range and payload for each of the retrofits
%
% Input:
%   ACFT                - Contains the data on the aircraft
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
% Output:
%   PayloadRange        - Structure with baseline, parallel, series,
%                         hydrogen and electric payload range data

% Unwrap structures
Weight = ACFT.Weight;
EnergyStorage = ACFT.EnergyStorage;
payload_kg = Weight.payload_kg;
MTOW_kg    = Weight.MTOW_kg;
avgasSpecificEnergy_kWhkg = EnergyStorage.avgasSpecificEnergy_kWhkg;
rhoCruise_kgm3 = Mission.Rho_kgm3.cruise;
cruiseSpeed_ms = MissionPerformance.cruiseSpeed_ms;

%% Baseline
maxFuelWeight_kg = EnergyStorage.Baseline.maxFuelWeight_kg;
fuelWeight_kg    = Weight.Baseline.fuelWeight_kg;
fuelAdditionWeight_kg     = maxFuelWeight_kg-fuelWeight_kg;
operationalEmptyWeight_kg = Weight.Baseline.operationalEmptyWeight_kg;
Baseline = calculate_payloadRangeParams(operationalEmptyWeight_kg,...
    maxFuelWeight_kg,fuelWeight_kg,fuelAdditionWeight_kg,payload_kg,MTOW_kg);
Baseline = calculate_rangeBaseline(Baseline,Mission,...
    MissionPerformance,ACFT);

%% Parallel
maxFuelWeight_kg = EnergyStorage.Parallel.maxFuelWeight_kg;
fuelWeight_kg    = Weight.Parallel.fuelWeight_kg;
fuelAdditionWeight_kg     = maxFuelWeight_kg-fuelWeight_kg;
operationalEmptyWeight_kg = Weight.Parallel.operationalEmptyWeight_kg;
Parallel = calculate_payloadRangeParams(operationalEmptyWeight_kg,...
    maxFuelWeight_kg,fuelWeight_kg,fuelAdditionWeight_kg,payload_kg,MTOW_kg);
 Parallel = calculate_rangeParallel(Parallel,Mission,...
    MissionPerformance,ACFT);

%% Series
maxFuelWeight_kg = EnergyStorage.Series.maxFuelWeight_kg;
fuelWeight_kg    = Weight.Series.fuelWeight_kg;
fuelAdditionWeight_kg     = maxFuelWeight_kg-fuelWeight_kg;
operationalEmptyWeight_kg = Weight.Series.operationalEmptyWeight_kg;
Series = calculate_payloadRangeParams(operationalEmptyWeight_kg,...
    maxFuelWeight_kg,fuelWeight_kg,fuelAdditionWeight_kg,payload_kg,MTOW_kg);
Series = calculate_rangeSeries(Series,Mission,MissionPerformance,ACFT);
%% Hydrogen
maxFuelWeight_kg = EnergyStorage.Hydrogen.maxFuelWeight_kg;
fuelWeight_kg    = Weight.Hydrogen.fuelWeight_kg;
fuelAdditionWeight_kg     = maxFuelWeight_kg-fuelWeight_kg;
operationalEmptyWeight_kg = Weight.Hydrogen.operationalEmptyWeight_kg;
Hydrogen = calculate_payloadRangeParams(operationalEmptyWeight_kg,...
    maxFuelWeight_kg,fuelWeight_kg,fuelAdditionWeight_kg,payload_kg,MTOW_kg);
Hydrogen = calculate_rangeHydrogen(Hydrogen,Mission,MissionPerformance,ACFT);
%% Electric
Electric.payloadWeight_kg = [payload_kg,payload_kg,0];
Electric.weight_kg = [MTOW_kg,MTOW_kg,MTOW_kg-payload_kg];
Electric = calculate_rangeElectric(Electric,Mission,MissionPerformance,ACFT);
%% Wrap the structure
PayloadRange.Baseline = Baseline;
PayloadRange.Parallel = Parallel;
PayloadRange.Series   = Series;
PayloadRange.Hydrogen = Hydrogen;
PayloadRange.Electric = Electric;
end

function Retrofit = calculate_payloadRangeParams(operationalEmptyWeight_kg,...
    maxFuelWeight_kg,fuelWeight_kg,fuelAdditionWeight_kg,payload_kg,MTOW_kg)
% Function:
%   calculate_payloadRangeParams
%
% Description: 
%   Calculate the payload, fuel and take-off weight for any general
%   retrofit to be used to calculate the payload range.
%
% Input:
%   operationalEmptyWeight_kg
%   maxFuelWeight_kg           - Maximum fuel weight inside tanks
%   fuelWeight_kg              - Maximum fuel when the aircraft is fully
%                                loaded
%   fuelAdditionWeight_kg      - Fuel that can be added to the aircarft to
%                                reach maximum capacity
%   payload_kg                 - Max. payload 
%   MTOW_kg                    - Maximum take-off weight
% Output:
%   Retrofit                   - Structure wth payloadWeight, fuel weight
%   and total take-off weight

Retrofit.payloadWeight_kg = [payload_kg,payload_kg,...
    payload_kg-fuelAdditionWeight_kg,0];
Retrofit.fuelWeight_kg    = [0,fuelWeight_kg,maxFuelWeight_kg,...
    maxFuelWeight_kg];
Retrofit.weight_kg = [MTOW_kg,operationalEmptyWeight_kg+fuelWeight_kg+...
    payload_kg,...
    operationalEmptyWeight_kg+fuelWeight_kg+Retrofit.payloadWeight_kg(3)+...
    fuelAdditionWeight_kg,...
    operationalEmptyWeight_kg+maxFuelWeight_kg];         % Flying Weight
end

function Baseline = calculate_rangeBaseline(Baseline,Mission,...
    MissionPerformance,ACFT)
% Function:
%   calculate_rangeBaseline(Baseline,Mission,...
%   MissionPerformance,ACFT)
%
% Description: 
%   Calculate the range for all four points on the payload range diagram
%   for the baseline
%
% Input:
%   Baseline            - Strucutre containing TO weight, fuel and payload
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft% Output:
% Ouptut:
%   Baseline            - Added range to the structure

avgasSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
Baseline.range_m = zeros(1,4);
fuelWeight_kg = Baseline.fuelWeight_kg;
weight_kg     = Baseline.weight_kg;
rhoCruise_kgm3 = Mission.Rho_kgm3.cruise;
cruiseSpeed_ms = MissionPerformance.cruiseSpeed_ms;

for i = 2:4
    fuel_kg = fuelWeight_kg(i);
    fuelEnergy_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
    weight_N= weight_kg(i)*9.81;
    cruisePower_kW = calculate_power(ACFT,weight_N,cruiseSpeed_ms,...
        rhoCruise_kgm3);
    % Estimate range
    tEstimate_s = fuelEnergy_kWh*0.3/cruisePower_kW*3600;
    % Define two crusie times, one above and another below
    t1_s = 0;           % Time below the actual max.
    t2_s = 0;           % Time above the actual max.
    % Calculate mission fuel consumption for the estimate
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    GuessBaseline = calculate_missionPerformanceBaseline(Mission,...
        MissionPerformance,ACFT,weight_N);
    fuelCons_kg = GuessBaseline.fuelConsTotal_kg;
    % Iterate until one value is above and another below
    while (t1_s == 0 || t2_s == 0)
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
            % Decrease guess by 1 hour
            tEstimate_s = tEstimate_s-3600;
        else
            t1_s = tEstimate_s;
            % Increase guess by 10 minutes
            tEstimate_s = tEstimate_s+3600;
        end
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessBaseline= calculate_missionPerformanceBaseline(Mission,...
            MissionPerformance,ACFT,weight_N);
        fuelCons_kg = GuessBaseline.fuelConsTotal_kg;       
    end
    % Use bisection method to find the range
    while abs(fuelCons_kg-fuel_kg)>0.01
        tEstimate_s = 0.5*(t1_s+t2_s);
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessBaseline = calculate_missionPerformanceBaseline(Mission,...
            MissionPerformance,ACFT,weight_N);
        fuelCons_kg = GuessBaseline.fuelConsTotal_kg;
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
        else
            t1_s = tEstimate_s;
        end
    end
    
    % Output range and parameters 
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    Mission = calculate_mission(Mission,MissionPerformance);
    Baseline.range_m(i) = Mission.distance.mission_m;

end
end
function Parallel = calculate_rangeParallel(Parallel,Mission,...
    MissionPerformance,ACFT)
% Function:
%   calculate_rangeParallel(Parallel,Mission,...
%   MissionPerformance,ACFT)
%
% Description: 
%   Calculate the range for all four points on the payload range diagram
%   for the parallel
%
% Input:
%   Parallel            - Strucutre containing TO weight, fuel and payload
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft% Output:
% Ouptut:
%   Parallel            - Added range to the structure

avgasSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
Parallel.range_m = zeros(1,4);
fuelWeight_kg = Parallel.fuelWeight_kg;
weight_kg     = Parallel.weight_kg;
rhoCruise_kgm3 = Mission.Rho_kgm3.cruise;
cruiseSpeed_ms = MissionPerformance.cruiseSpeed_ms;
for i = 2:4
    fuel_kg = fuelWeight_kg(i);
    fuelEnergy_kWh = fuel_kg*avgasSpecificEnergy_kWhkg;
    weight_N= weight_kg(i)*9.81;
    cruisePower_kW = calculate_power(ACFT,weight_N,cruiseSpeed_ms,...
        rhoCruise_kgm3);
    % Estimate range
    tEstimate_s = fuelEnergy_kWh*0.3/cruisePower_kW*3600;
    % Define two crusie times, one above and another below
    t1_s = 0;           % Time below the actual max.
    t2_s = 0;           % Time above the actual max.
    % Calculate mission fuel consumption for the estimate
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    GuessParallel = calculate_missionPerformanceParallel(Mission,...
        MissionPerformance,ACFT,weight_N,1);
    fuelCons_kg = GuessParallel.fuelConsTotal_kg;
    % Iterate until one value is above and another below
    while (t1_s == 0 || t2_s == 0)
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
            % Decrease guess by 1 hour
            tEstimate_s = tEstimate_s-3600;
        else
            t1_s = tEstimate_s;
            % Increase guess by 10 minutes
            tEstimate_s = tEstimate_s+3600;
        end
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessParallel = calculate_missionPerformanceParallel(Mission,...
            MissionPerformance,ACFT,weight_N,1);
        fuelCons_kg = GuessParallel.fuelConsTotal_kg;       
    end
    % Use bisection method to find the range
    while abs(fuelCons_kg-fuel_kg)>0.01
        tEstimate_s = 0.5*(t1_s+t2_s);
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessParallel = calculate_missionPerformanceParallel(Mission,...
            MissionPerformance,ACFT,weight_N,1);
        fuelCons_kg = GuessParallel.fuelConsTotal_kg;
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
        else
            t1_s = tEstimate_s;
        end
    end
    
    % Output range and parameters 
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    Mission = calculate_mission(Mission,MissionPerformance);
    Parallel.range_m(i) = Mission.distance.mission_m;

end
end

function Series = calculate_rangeSeries(Series,Mission,...
    MissionPerformance,ACFT)
% Function:
%   calculate_rangeSeries(Series,Mission,...
%   MissionPerformance,ACFT)
%
% Description: 
%   Calculate the range for all four points on the payload range diagram
%   for the Series
%
% Input:
%   Series            - Strucutre containing TO weight, fuel and payload
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft% Output:
% Ouptut:
%   Series            - Added range to the structure

hybridCruise = 1;       % Set hybrid cruise to ON as default
rechargeSeries = 0;     % Set recharge series to OFF as default
fuelSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
Series.range_m = zeros(1,4);
fuelWeight_kg = Series.fuelWeight_kg;
weight_kg     = Series.weight_kg;
rhoCruise_kgm3 = Mission.Rho_kgm3.cruise;
cruiseSpeed_ms = MissionPerformance.cruiseSpeed_ms;
for i = 2:4
    fuel_kg = fuelWeight_kg(i);
    fuelEnergy_kWh = fuel_kg*fuelSpecificEnergy_kWhkg;
    weight_N= weight_kg(i)*9.81;
    cruisePower_kW = calculate_power(ACFT,weight_N,cruiseSpeed_ms,...
        rhoCruise_kgm3);
    % Estimate range
    tEstimate_s = fuelEnergy_kWh*0.3/cruisePower_kW*3600;
    % Define two crusie times, one above and another below
    t1_s = 0;           % Time below the actual max.
    t2_s = 0;           % Time above the actual max.
    % Calculate mission fuel consumption for the estimate
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    GuessSeries = calculate_missionPerformanceSeries(Mission,...
        MissionPerformance,ACFT,weight_N,hybridCruise,rechargeSeries);
    fuelCons_kg = GuessSeries.fuelConsTotal_kg;
    % Iterate until one value is above and another below
    while (t1_s == 0 | t2_s == 0)
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
            % Decrease guess by 1 hour
            tEstimate_s = tEstimate_s-3600;
        else
            t1_s = tEstimate_s;
            % Increase guess by 10 minutes
            tEstimate_s = tEstimate_s+3600;
        end
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessSeries = calculate_missionPerformanceSeries(Mission,...
            MissionPerformance,ACFT,weight_N,hybridCruise,rechargeSeries);
        fuelCons_kg = GuessSeries.fuelConsTotal_kg;       
    end
    % Use bisection method to find the range
    while abs(fuelCons_kg-fuel_kg)>0.01;
        tEstimate_s = 0.5*(t1_s+t2_s);
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessSeries = calculate_missionPerformanceSeries(Mission,...
            MissionPerformance,ACFT,weight_N,hybridCruise,rechargeSeries);
        fuelCons_kg = GuessSeries.fuelConsTotal_kg;
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
        else
            t1_s = tEstimate_s;
        end
    end
    
    % Output range and parameters 
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    Mission = calculate_mission(Mission,MissionPerformance);
    Series.range_m(i) = Mission.distance.mission_m;

end
end

function Hydrogen = calculate_rangeHydrogen(Hydrogen,Mission,...
    MissionPerformance,ACFT)
% Function:
%   calculate_rangeHydrogen(Hydrogen,Mission,...
%   MissionPerformance,ACFT)
%
% Description: 
%   Calculate the range for all four points on the payload range diagram
%   for the Hydrogen retrofit
%
% Input:
%   Hydrogen            - Structure containing TO weight, fuel and payload
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft% Output:
% Ouptut:
%   Hydrogen            - Added range to the structure

hybridCruise = 1;       % Set hybrid cruise to ON as default
fuelSpecificEnergy_kWhkg = ACFT.EnergyStorage.avgasSpecificEnergy_kWhkg;
Hydrogen.range_m = zeros(1,4);
fuelWeight_kg = Hydrogen.fuelWeight_kg;
weight_kg     = Hydrogen.weight_kg;
rhoCruise_kgm3 = Mission.Rho_kgm3.cruise;
cruiseSpeed_ms = MissionPerformance.cruiseSpeed_ms;
for i = 2:4
    fuel_kg = fuelWeight_kg(i);
    fuelEnergy_kWh = fuel_kg*fuelSpecificEnergy_kWhkg;
    weight_N= weight_kg(i)*9.81;
    cruisePower_kW = calculate_power(ACFT,weight_N,cruiseSpeed_ms,...
        rhoCruise_kgm3);
    % Estimate range
    tEstimate_s = fuelEnergy_kWh*0.3/cruisePower_kW*3600;
    % Define two crusie times, one above and another below
    t1_s = 0;           % Time below the actual max.
    t2_s = 0;           % Time above the actual max.
    % Calculate mission fuel consumption for the estimate
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    GuessHydrogen = calculate_missionPerformanceHydrogen(Mission,...
        MissionPerformance,ACFT,weight_N,hybridCruise);
    fuelCons_kg = GuessHydrogen.fuelConsTotal_kg;
    % Iterate until one value is above and another below
    while (t1_s == 0 | t2_s == 0)
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
            % Decrease guess by 1 hour
            tEstimate_s = tEstimate_s-3600;
        else
            t1_s = tEstimate_s;
            % Increase guess by 10 minutes
            tEstimate_s = tEstimate_s+3600;
        end
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessHydrogen = calculate_missionPerformanceHydrogen(Mission,...
            MissionPerformance,ACFT,weight_N,hybridCruise);
        fuelCons_kg = GuessHydrogen.fuelConsTotal_kg;       
    end
    % Use bisection method to find the range
    while abs(fuelCons_kg-fuel_kg)>0.01;
        tEstimate_s = 0.5*(t1_s+t2_s);
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessHydrogen = calculate_missionPerformanceHydrogen(Mission,...
            MissionPerformance,ACFT,weight_N,hybridCruise);
        fuelCons_kg = GuessHydrogen.fuelConsTotal_kg;
        if fuelCons_kg > fuel_kg
            t2_s = tEstimate_s;
        else
            t1_s = tEstimate_s;
        end
    end
    
    % Output range and parameters 
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    Mission = calculate_mission(Mission,MissionPerformance);
    Hydrogen.range_m(i) = Mission.distance.mission_m;

end
end

function Electric = calculate_rangeElectric(Electric,Mission,...
    MissionPerformance,ACFT)
% Function:
%   calculate_rangeElectric(Electric,Mission,...
%   MissionPerformance,ACFT)
%
% Description: 
%   Calculate the range for all four points on the payload range diagram
%   for the Electric retrofit
%
% Input:
%   Electric            - Structure containing TO weight, fuel and payload
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraf
% Ouptut:
%   Electric            - Added range to the structure

batteryEnergy_kWh = ACFT.EnergyStorage.Electric.batteryEnergy_kWh; 
number            = length(batteryEnergy_kWh);
Electric.range_m = zeros(number,3);
weight_kg     = Electric.weight_kg;
rhoCruise_kgm3 = Mission.Rho_kgm3.cruise;
cruiseSpeed_ms = MissionPerformance.cruiseSpeed_ms;
for n = 1:number
    for i = 2:3
        maxBatteryEnergy_kWh = batteryEnergy_kWh(n); 
        weight_N= weight_kg(i)*9.81;
        cruisePower_kW = calculate_power(ACFT,weight_N,cruiseSpeed_ms,...
            rhoCruise_kgm3);
        % Estimate range
        tEstimate_s = maxBatteryEnergy_kWh*0.9*0.95/cruisePower_kW*3600;
        % Define two crusie times, one above and another below
        t1_s = 0;           % Time below the actual max.
        t2_s = 0;           % Time above the actual max.
        % Calculate mission fuel consumption for the estimate
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessElectric = calculate_missionPerformanceElectric(Mission,...
            MissionPerformance,ACFT,weight_N);
        energyCons_kWh = GuessElectric.battEnergyTotal_kWh;
        % Iterate until one value is above and another below
    while (t1_s == 0 | t2_s == 0)
        if  energyCons_kWh > maxBatteryEnergy_kWh
            t2_s = tEstimate_s;
            % Decrease guess by 1 hour
            tEstimate_s = tEstimate_s-3600;
        else
            t1_s = tEstimate_s;
            % Increase guess by 10 minutes
            tEstimate_s = tEstimate_s+3600;
        end
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessElectric = calculate_missionPerformanceElectric(Mission,...
            MissionPerformance,ACFT,weight_N);
        energyCons_kWh = GuessElectric.battEnergyTotal_kWh;       
    end
    % Use bisection method to find the range
    while abs(energyCons_kWh-maxBatteryEnergy_kWh)>0.01;
        tEstimate_s = 0.5*(t1_s+t2_s);
        Mission.time.cruise_s     = tEstimate_s;
        Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
        GuessElectric = calculate_missionPerformanceElectric(Mission,...
            MissionPerformance,ACFT,weight_N);
        energyCons_kWh = GuessElectric.battEnergyTotal_kWh;
        if energyCons_kWh >maxBatteryEnergy_kWh
            t2_s = tEstimate_s;
        else
            t1_s = tEstimate_s;
        end
    end
    
    % Output range and parameters 
    Mission.time.cruise_s     = tEstimate_s;
    Mission.distance.cruise_m = tEstimate_s*cruiseSpeed_ms;
    Mission = calculate_mission(Mission,MissionPerformance);
    Electric.range_m(n,i) = Mission.distance.mission_m;

end
end
end




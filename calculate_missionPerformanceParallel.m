function Parallel = calculate_missionPerformanceParallel(Mission,MissionPerformance,ACFT,weight_N,hybridCruise)
% Function:
%   calculate_missionPerformaceParallel
%
% Description: 
%   Calculate the fuel and electric energy consumption for a particular mission
%   profile for a parallel retrofit with and without hybrid cruise
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
% Output:
%   Parallel            - Structure with fuel, power and energy consumption
%

% Unwrap from Structures required variables
electricMotorEff                = ACFT.Propulsion.Parallel.electricMotorEff;
batteryEff                      = ACFT.Propulsion.Parallel.batteryEff;
maxICE_Power_kW                 = ACFT.Propulsion.Parallel.maxICE_Power_kW;
maxElectricPower_kW             = ACFT.Propulsion.Parallel.maxElectricPower_kW;
ICE_maxSpecificFuelCons_kgkWh   = ACFT.Propulsion.Parallel.ICE_maxSpecificFuelCons_kgkWh;
ICE_specificFuelCons_kgkWh      = ACFT.Propulsion.Parallel.ICE_specificFuelCons_kgkWh;
rhoFuel_kgm3                    = ACFT.Propulsion.Parallel.rhoFuel_kgm3;
ShaftPower_kW                   = MissionPerformance.ShaftPower_kW;
%% Parallel
% Taxi
BatteryEnergy_kWh.taxi          = ACFT.EnergyStorage.Parallel.batteryEnergy_kWh;
ICE_power_kW.taxi= 0;
Weight_N.taxi    = weight_N; 
FuelCons_USgal.taxi = 0;
ElectricPower_kW.taxi   = ShaftPower_kW.taxi;           % Note different nomenclature in code vs in report ElectricPower is that delivered by Electric Motors!
ElectricEnergy_kWh.taxi = ElectricPower_kW.taxi*Mission.time.taxi_s/3600/batteryEff/electricMotorEff; 
% TO
BatteryEnergy_kWh.takeOff  = BatteryEnergy_kWh.taxi-ElectricEnergy_kWh.taxi;
ICE_power_kW.takeOff= maxICE_Power_kW;
Weight_N.takeOff    = Weight_N.taxi-FuelCons_USgal.taxi*3.785/1000*rhoFuel_kgm3*9.81; 
FuelCons_USgal.takeOff = Mission.time.takeOff_s*ICE_maxSpecificFuelCons_kgkWh*ICE_power_kW.takeOff/rhoFuel_kgm3*1000/3.785/3600;

ElectricPower_kW.takeOff   = maxElectricPower_kW;
ElectricEnergy_kWh.takeOff = ElectricPower_kW.takeOff*Mission.time.takeOff_s/3600/batteryEff/electricMotorEff;

% Climb
BatteryEnergy_kWh.climb  = BatteryEnergy_kWh.takeOff-ElectricEnergy_kWh.takeOff;
ICE_power_kW.climb  = maxICE_Power_kW;
Weight_N.climb    = Weight_N.takeOff-FuelCons_USgal.takeOff*3.785/1000*rhoFuel_kgm3*9.81;
FuelCons_USgal.climb = Mission.time.climb_s*ICE_maxSpecificFuelCons_kgkWh*ICE_power_kW.climb/rhoFuel_kgm3*1000/3.785/3600;

ElectricPower_kW.climb     = maxElectricPower_kW;
ElectricEnergy_kWh.climb   = ElectricPower_kW.climb*Mission.time.climb_s/3600/batteryEff/electricMotorEff;
% Cruise
BatteryEnergy_kWh.cruise = BatteryEnergy_kWh.climb-ElectricEnergy_kWh.climb;
Weight_N.cruise   = Weight_N.climb-FuelCons_USgal.climb*3.785/1000*rhoFuel_kgm3*9.81;

[FuelCons_USgal.cruise,ICE_power_kW.cruise,Weight_N.descent,...
    ElectricPower_kW.cruise,ElectricEnergy_kWh.cruise] = ...
    calculate_cruisePerformanceParallel(Mission,MissionPerformance,ACFT,...
    Weight_N.cruise,BatteryEnergy_kWh.cruise,hybridCruise);

% Descent
BatteryEnergy_kWh.descent = BatteryEnergy_kWh.cruise-ElectricEnergy_kWh.cruise;
ICE_power_kW.descent = ShaftPower_kW.descent;
ElectricPower_kW.descent = 0;
FuelCons_USgal.descent = ICE_power_kW.descent*ICE_specificFuelCons_kgkWh*Mission.time.descent_s/rhoFuel_kgm3*1000/3.785/3600;
ElectricEnergy_kWh.descent = ElectricPower_kW.descent*Mission.time.descent_s/3600;

% Reserve Climb
BatteryEnergy_kWh.reserveClimb = BatteryEnergy_kWh.descent-ElectricEnergy_kWh.descent;
ICE_power_kW.reserveClimb  = maxICE_Power_kW;
Weight_N.reserveClimb      = Weight_N.descent-FuelCons_USgal.descent*3.785/1000*rhoFuel_kgm3*9.81;  
FuelCons_USgal.reserveClimb = Mission.time.reserveClimb_s*ICE_maxSpecificFuelCons_kgkWh*ICE_power_kW.climb/rhoFuel_kgm3*1000/3.785/3600;

ElectricPower_kW.reserveClimb     = maxElectricPower_kW;
ElectricEnergy_kWh.reserveClimb   = ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s/3600/batteryEff/electricMotorEff;
% Reserve Loiter
BatteryEnergy_kWh.reserveLoiter   = BatteryEnergy_kWh.reserveClimb-ElectricEnergy_kWh.reserveClimb;
Weight_N.reserveLoiter   = Weight_N.reserveClimb-FuelCons_USgal.reserveClimb*3.785/1000*rhoFuel_kgm3*9.81;
% Use the baseline loiter function for the parallel, this assumes that the
% parallel will not use electric energy during this segment
[FuelCons_USgal.reserveLoiter,ICE_power_kW.reserveLoiter,Weight_N.reserveDescent] = calculate_reserveLoiterPerformanceParallel(Mission,MissionPerformance,ACFT,Weight_N.reserveLoiter);

ElectricPower_kW.reserveLoiter = 0;
ElectricEnergy_kWh.reserveLoiter = 0;

% Reserve Descent
BatteryEnergy_kWh.reserveDescent = BatteryEnergy_kWh.reserveLoiter-ElectricEnergy_kWh.reserveLoiter;
ICE_power_kW.reserveDescent = ShaftPower_kW.descent;
FuelCons_USgal.reserveDescent = ICE_power_kW.reserveDescent*ICE_specificFuelCons_kgkWh*Mission.time.reserveDescent_s/rhoFuel_kgm3*1000/3.785/3600;

ElectricPower_kW.reserveDescent = 0;
ElectricEnergy_kWh.reserveDescent = 0;
% Average reserve power
ICE_power_kW.reserve  = (ICE_power_kW.reserveClimb*Mission.time.reserveClimb_s+ICE_power_kW.reserveLoiter*Mission.time.reserveLoiter_s)...
    /(Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);
ElectricPower_kW.reserve = (ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s+ElectricPower_kW.reserveLoiter*Mission.time.reserveLoiter_s)...
    /(Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);

fuelConsMission_USgal = FuelCons_USgal.taxi+FuelCons_USgal.takeOff+FuelCons_USgal.climb...
    +FuelCons_USgal.cruise+FuelCons_USgal.descent;
fuelConsReserve_USgal = FuelCons_USgal.reserveClimb+FuelCons_USgal.reserveLoiter+FuelCons_USgal.reserveDescent;

elecEnergyMission_kWh = ElectricEnergy_kWh.taxi + ElectricEnergy_kWh.takeOff + ElectricEnergy_kWh.climb + ElectricEnergy_kWh.cruise + ElectricEnergy_kWh.descent;
elecEnergyReserve_kWh = ElectricEnergy_kWh.reserveClimb + ElectricEnergy_kWh.reserveLoiter + ElectricEnergy_kWh.reserveDescent;

% Pack into a structure
Parallel.ICE_power_kW   = ICE_power_kW;
Parallel.weight_N       = Weight_N;
Parallel.fuelCons_USgal = FuelCons_USgal; 
Parallel.fuelConsMission_USgal = fuelConsMission_USgal;
Parallel.fuelConsReserve_USgal = fuelConsReserve_USgal;

Parallel.electricPower_kW     = ElectricPower_kW;
Parallel.electricEnergy_kWh    = ElectricEnergy_kWh;
Parallel.batteryEnergy_kWh     = BatteryEnergy_kWh;
Parallel.electricEnergyMission_kWh = elecEnergyMission_kWh;
Parallel.electricEnergyReserve_kWh = elecEnergyReserve_kWh;

Parallel.fuelConsMission_kg    = Parallel.fuelConsMission_USgal*0.003785*...
    rhoFuel_kgm3;
Parallel.fuelConsReserve_kg    = Parallel.fuelConsReserve_USgal*0.003785*...
    rhoFuel_kgm3;
Parallel.fuelConsTotal_kg = Parallel.fuelConsMission_kg+...
    Parallel.fuelConsReserve_kg;
end
function [fuelCruise_USgal,avgICEPower_kW,weightFinal_N,avgElecPower_kW,elecEnergy_kWh] =  calculate_cruisePerformanceParallel(Mission,MissionPerformance,ACFT,weight_N,initialBattEnergy_kWh,hybridCruise);
% Function:
%   calculate_missionPerformaceParallel
%
% Description: 
%   Calculate the fuel and electric energy consumption for a particular mission
%   profile for a parallel retrofit with and without hybrid cruise
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
%   hybridCruise        - Option to choose whether electric energy is used
%                         during cruise or not
% Output:
%   Parallel            - Structure with fuel, power and energy consumption
%


% Unwrap the mission and mission performance structures
timeCruise_s     = Mission.time.cruise_s;
altitudeCruise_m = Mission.altitude.cruise_m;
cruiseSpeed_ms   = MissionPerformance.cruiseSpeed_ms;
distanceCruise_m = Mission.distance.cruise_m;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Baseline.propEff;
ICE_specificFuelCons_kgkWh = ACFT.Propulsion.Parallel.ICE_specificFuelCons_kgkWh;
maxElectricPower_kW        = ACFT.Propulsion.Parallel.maxElectricPower_kW;
electricMotorEff           = ACFT.Propulsion.Parallel.electricMotorEff;
batteryEff                 = ACFT.Propulsion.Parallel.batteryEff;
rhoFuel_kgm3               = ACFT.Propulsion.Parallel.rhoFuel_kgm3;
batteryEnergy_kWh          = ACFT.EnergyStorage.Parallel.batteryEnergy_kWh;


% Define number of descritizations
n = 2;
% Start the matrices
timeSegment_s = ones(1,n)*timeCruise_s/n;       % Allows for variable time segments
[~, ~, ~, rho_kgm3] = atmosisa(altitudeCruise_m);
weightSegment_N = ones(1,n)*weight_N;
lift_N        = ones(1,n);
cl            = ones(1,n);
cd            = ones(1,n);
drag_N        = ones(1,n);
LoD           = ones(1,n);
fuel_USgal    = ones(1,n); 
fuel_N        = ones(1,n);
power_kW      = ones(1,n);      % Shaft power output (Needs multiplying by Prop eff)
powerICE_kW   = ones(1,n);      % Power delivered at the shaft by ICE engine (needs multiplying by prop eff)
powerElec_kW  = ones(1,n);      % Power delivered to electric motor (needs multiplying by motor eff and prop eff to get electric motor shaft power)
electricEnergy_kWh = ones(1,n);
if hybridCruise == 1
        reserveClimbEnergy_kWh = maxElectricPower_kW*...
        Mission.time.reserveClimb_s/3600/batteryEff/electricMotorEff;
    energyLeftOver_kWh=initialBattEnergy_kWh-0.05*batteryEnergy_kWh-...
        reserveClimbEnergy_kWh;  
    powerElec_kW = (energyLeftOver_kWh)*batteryEff*electricMotorEff*3600/timeCruise_s*ones(1,n);
    if powerElec_kW > maxElectricPower_kW
        powerElec_kW = maxElectricPower_kW*ones(1,n);
    end
    for i = 1:n
        if i ~= 1
            weightSegment_N(i) = weightSegment_N(i-1)-fuel_N(i-1);
        end
        lift_N(i) = weightSegment_N(i);
        cl(i)     = lift_N(i)/(0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2);
        cd(i)     = cd0 + k*cl(i)^2;
        LoD(i)    = cl(i)/cd(i);
        drag_N(i) = 0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2*cd(i);
        power_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;
        powerICE_kW(i)= power_kW(i)-powerElec_kW(i);
        fuel_N(i) = power_kW(i)*ICE_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
        fuel_USgal(i) = fuel_N(i)/9.81/rhoFuel_kgm3*1000/3.785;
        electricEnergy_kWh(i) = powerElec_kW(i)*timeSegment_s(i)/3600/batteryEff/electricMotorEff;
        
    end
    
else
    powerElec_kW = 0*powerElec_kW;
    for i = 1:n
        if i ~= 1
            weightSegment_N(i) = weightSegment_N(i-1)-fuel_N(i-1);
        end
        lift_N(i) = weightSegment_N(i);
        cl(i)     = lift_N(i)/(0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2);
        cd(i)     = cd0 + k*cl(i)^2;
        LoD(i)    = cl(i)/cd(i);
        drag_N(i) = 0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2*cd(i);
        power_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;
        powerICE_kW(i)= power_kW(i)-powerElec_kW(i);
        fuel_N(i) = powerICE_kW(i)*ICE_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
        fuel_USgal(i) = fuel_N(i)/9.81/rhoFuel_kgm3*1000/3.785;
        electricEnergy_kWh(i) = powerElec_kW(i)*timeSegment_s(i)/3600/batteryEff/electricMotorEff;
        
    end
end
fuelCruise_N      = sum(fuel_N);
fuelCruise_USgal  = sum(fuel_USgal);
weightFinal_N     = weight_N-fuelCruise_N;
avgICEPower_kW    = sum(powerICE_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
elecEnergy_kWh    = sum(electricEnergy_kWh);
end
function [fuelCruise_USgal,avgPower_kW,weightFinal_N] = calculate_reserveLoiterPerformanceParallel(Mission,MissionPerformance,ACFT,weight_N)
% Function:
%   calculate_cruisePerformanceBaseline
%
% Description: 
%   Calculate the fuel consumption for a cruise segment of specified lenght
%   and airspeed for a baseline conventional aircraft.
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight of cruise in N
% Output:
%   fuelCruise_USGal    - fuel consumption in US Gal
%   avgPower_kW         - average power output from engine in kW
%   weightFinal_N       - weight after end of cruise

checkBreguet = 0;       % 1 for Breguet check; 0 for no Breguet
% Unwrap the mission and mission performance structures
timeCruise_s     = Mission.time.reserveLoiter_s;
altitudeCruise_m = Mission.altitude.reserve_m;
cruiseSpeed_ms   = MissionPerformance.reserveSpeed_ms;
distanceCruise_m = cruiseSpeed_ms*timeCruise_s;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Parallel.propEff;
rhoFuel_kgm3     = ACFT.Propulsion.Parallel.rhoFuel_kgm3;
ICE_specificFuelCons_kgkWh = ACFT.Propulsion.Parallel.ICE_specificFuelCons_kgkWh;
% Define number of descritizations
n = 1;
% Start the matrices
timeSegment_s = ones(1,n)*timeCruise_s/n;       % Allows for variable time segments
[~, ~, ~, rho_kgm3] = atmosisa(altitudeCruise_m);
weightSegment_N = ones(1,n)*weight_N;
lift_N        = ones(1,n);
cl            = ones(1,n);
cd            = ones(1,n);
drag_N        = ones(1,n);
LoD           = ones(1,n);
fuel_USgal    = ones(1,n); 
fuel_N        = ones(1,n);
power_kW      = ones(1,n);
for i = 1:n
    if i ~= 1
        weightSegment_N(i) = weightSegment_N(i-1)-fuel_N(i-1);
    end
    lift_N(i) = weightSegment_N(i);
    cl(i)     = lift_N(i)/(0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2);
    cd(i)     = cd0 + k*cl(i)^2;
    LoD(i)    = cl(i)/cd(i); 
    drag_N(i) = 0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2*cd(i);
    power_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;
    fuel_N(i) = power_kW(i)*ICE_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
    fuel_USgal(i) = fuel_N(i)/9.81/rhoFuel_kgm3*1000/3.785;
end
fuelCruise_N  = sum(fuel_N);
fuelCruise_USgal  = sum(fuel_USgal);
weightFinal_N = weight_N-fuelCruise_N;
avgPower_kW   = sum(power_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
% Validate with Breguet Range Equation
if checkBreguet == 1
    weightRatio = exp(-distanceCruise_m*9.81*ICE_specificFuelCons_kgkWh/1000/3600/(propEff*LoD(1)));
    breguetFuel_N = weight_N*(1-weightRatio);
    breguetFuel_USgal = breguetFuel_N/9.81/rhoFuel_kgm3*1000/3.785;

    percentageDiff = (breguetFuel_N-fuelCruise_N)/fuelCruise_N*100;
    fprintf('%18s %18s %18s\n','Discrete (Gal)','Breguet (Gal)','% Diff');
    fprintf('%18.1f %18.1f %18.1f\n',fuelCruise_USgal,breguetFuel_USgal,percentageDiff);
end
end


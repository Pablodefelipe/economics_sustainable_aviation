function Baseline = calculate_missionPerformanceBaseline(Mission,MissionPerformance,ACFT,weight_N)
% Function:
%   calculate_missionPerformaceBaseline
%
% Description: 
%   Calculate the fuel and energy consumption for a particular mission
%   profile for a baseline retrofit (ICE Engines only)
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
% Output:
%   MissionPerformance  - fuel and energy consumption

% Unwrap structures;
rhoFuel_kgm3                = ACFT.Propulsion.Baseline.rhoFuel_kgm3;
ShaftPower_kW               = MissionPerformance.ShaftPower_kW;
ICE_specificFuelCons_kgkWh  = ACFT.Propulsion.Baseline.ICE_specificFuelCons_kgkWh;
ICE_maxSpecificFuelCons_kgkWh=ACFT.Propulsion.Baseline.ICE_maxSpecificFuelCons_kgkWh;

%% Baseline
% Taxi
Baseline.ICE_power_kW.taxi= 0;
Baseline.weight_N.taxi    = weight_N; 
%Baseline.fuelCons_USgal.taxi = Mission.time.taxi_s*ACFT.Propulsion.Baseline.ICE_idleFuelCons_USgalh/3600;
Baseline.fuelCons_USgal.taxi = Mission.time.taxi_s*ShaftPower_kW.taxi*....
    ICE_maxSpecificFuelCons_kgkWh/rhoFuel_kgm3*1000/3.785/3600;

% TO
Baseline.ICE_power_kW.takeOff= ACFT.Propulsion.Baseline.maxICE_Power_kW;
Baseline.weight_N.takeOff    = Baseline.weight_N.taxi-Baseline.fuelCons_USgal.taxi*3.785/1000*rhoFuel_kgm3 *9.81; 
Baseline.fuelCons_USgal.takeOff = Mission.time.takeOff_s*ACFT.Propulsion.Baseline.ICE_maxSpecificFuelCons_kgkWh*Baseline.ICE_power_kW.takeOff/rhoFuel_kgm3 *1000/3.785/3600;

% Climb
Baseline.ICE_power_kW.climb  = ShaftPower_kW.climb;
Baseline.weight_N.climb    = Baseline.weight_N.takeOff-Baseline.fuelCons_USgal.takeOff*3.785/1000*rhoFuel_kgm3 *9.81;
Baseline.fuelCons_USgal.climb = Mission.time.climb_s*ACFT.Propulsion.Baseline.ICE_maxSpecificFuelCons_kgkWh*Baseline.ICE_power_kW.climb/rhoFuel_kgm3 *1000/3.785/3600;

% Cruise
Baseline.weight_N.cruise   = Baseline.weight_N.climb-Baseline.fuelCons_USgal.climb*3.785/1000*rhoFuel_kgm3 *9.81;
[Baseline.fuelCons_USgal.cruise,Baseline.ICE_power_kW.cruise,Baseline.weight_N.descent] = calculate_cruisePerformanceBaseline(Mission,MissionPerformance,ACFT,Baseline.weight_N.cruise);

% Descent
Baseline.ICE_power_kW.descent = ShaftPower_kW.descent;
Baseline.fuelCons_USgal.descent = Baseline.ICE_power_kW.descent*ICE_specificFuelCons_kgkWh*Mission.time.descent_s/rhoFuel_kgm3*1000/3.785/3600;

% Reserve Climb
Baseline.ICE_power_kW.reserveClimb  = ACFT.Propulsion.Baseline.maxICE_Power_kW;
Baseline.weight_N.reserveClimb    = Baseline.weight_N.descent-Baseline.fuelCons_USgal.descent*3.785/1000*rhoFuel_kgm3*9.81;  
Baseline.fuelCons_USgal.reserveClimb = Mission.time.reserveClimb_s*ACFT.Propulsion.Baseline.ICE_maxSpecificFuelCons_kgkWh*Baseline.ICE_power_kW.climb/rhoFuel_kgm3*1000/3.785/3600;

% Reserve Loiter
Baseline.weight_N.reserveLoiter   = Baseline.weight_N.reserveClimb-Baseline.fuelCons_USgal.reserveClimb*3.785/1000*rhoFuel_kgm3*9.81;
[Baseline.fuelCons_USgal.reserveLoiter,Baseline.ICE_power_kW.reserveLoiter,Baseline.weight_N.reserveDescent] = calculate_reserveLoiterPerformanceBaseline(Mission,MissionPerformance,ACFT,Baseline.weight_N.reserveLoiter);

% Reserve Descent
Baseline.ICE_power_kW.reserveDescent = ShaftPower_kW.descent;
Baseline.fuelCons_USgal.reserveDescent = Baseline.ICE_power_kW.reserveDescent*ICE_specificFuelCons_kgkWh*Mission.time.reserveDescent_s/rhoFuel_kgm3*1000/3.785/3600;
% Average reserve power
Baseline.ICE_power_kW.reserve  = (Baseline.ICE_power_kW.reserveClimb*Mission.time.reserveClimb_s+Baseline.ICE_power_kW.reserveLoiter*Mission.time.reserveLoiter_s)...
    /(Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);

Baseline.fuelConsMission_USgal = Baseline.fuelCons_USgal.taxi+Baseline.fuelCons_USgal.takeOff+Baseline.fuelCons_USgal.climb...
    +Baseline.fuelCons_USgal.cruise+Baseline.fuelCons_USgal.descent;
Baseline.fuelConsReserve_USgal = Baseline.fuelCons_USgal.reserveClimb+Baseline.fuelCons_USgal.reserveLoiter+Baseline.fuelCons_USgal.reserveDescent;
Baseline.fuelConsMission_kg    = Baseline.fuelConsMission_USgal*0.003785*...
    rhoFuel_kgm3;
Baseline.fuelConsReserve_kg    = Baseline.fuelConsReserve_USgal*0.003785*...
    rhoFuel_kgm3;
Baseline.fuelConsTotal_kg = Baseline.fuelConsMission_kg+...
    Baseline.fuelConsReserve_kg;
end

function [fuelCruise_USgal,avgPower_kW,weightFinal_N] = calculate_reserveLoiterPerformanceBaseline(Mission,MissionPerformance,ACFT,weight_N)
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
rhoFuel_kgm3     = ACFT.Propulsion.Baseline.rhoFuel_kgm3;
propEff          = ACFT.Propulsion.Baseline.propEff;
ICE_specificFuelCons_kgkWh = ACFT.Propulsion.Baseline.ICE_specificFuelCons_kgkWh;
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

function [fuelCruise_USgal,avgPower_kW,weightFinal_N] = calculate_cruisePerformanceBaseline(Mission,MissionPerformance,ACFT,weight_N)
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
timeCruise_s     = Mission.time.cruise_s;
altitudeCruise_m = Mission.altitude.cruise_m;
cruiseSpeed_ms   = MissionPerformance.cruiseSpeed_ms;
distanceCruise_m = Mission.distance.cruise_m;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;
rhoFuel_kgm3     = ACFT.Propulsion.Baseline.rhoFuel_kgm3;

propEff          = ACFT.Propulsion.Baseline.propEff;
ICE_specificFuelCons_kgkWh = ACFT.Propulsion.Baseline.ICE_specificFuelCons_kgkWh;
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

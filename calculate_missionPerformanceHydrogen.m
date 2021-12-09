function Hydrogen = calculate_missionPerformanceHydrogen(Mission,MissionPerformance,ACFT,weight_N,hybridCruise)
% Function:
%   calculate_missionPerformaceHydrogen
%
% Description: 
%   Calculate the fuel and electric energy consumption for a particular mission
%   profile for a hydrogen retrofit with and without hybrid cruise
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
electricMotorEff                = ACFT.Propulsion.Hydrogen.electricMotorEff;
batteryEff                      = ACFT.Propulsion.Hydrogen.batteryEff;
maxFC_Power_kW                  = ACFT.Propulsion.Hydrogen.maxFC_Power_kW;
maxElectricPower_kW             = ACFT.Propulsion.Hydrogen.maxElectricPower_kW;
FC_specificFuelCons_kgkWh       = ACFT.Propulsion.Hydrogen.FC_specificFuelCons_kgkWh;
rhoFuel_kgm3                    = ACFT.Propulsion.Hydrogen.rhoFuel_kgm3;
ShaftPower_kW                   = MissionPerformance.ShaftPower_kW;
time                            = Mission.time;
%% Hydrogen
% Taxi
BatteryEnergy_kWh.taxi          = ACFT.EnergyStorage.Parallel.batteryEnergy_kWh; % Remaining Battery energy at start of segment
FC_power_kW.taxi= 0;
Weight_N.taxi    = weight_N; 
FuelCons_kg.taxi = 0;
FCEnergy_kWh.taxi    = FC_power_kW.taxi*Mission.time.taxi_s;

ElectricPower_kW.taxi   = ShaftPower_kW.taxi;
BattPower_kW.taxi       = ElectricPower_kW.taxi/electricMotorEff-FC_power_kW.taxi;
BattEnergy_kWh.taxi = (ElectricPower_kW.taxi*Mission.time.taxi_s/3600/electricMotorEff-FCEnergy_kWh.taxi)/batteryEff; 
% TO
BatteryEnergy_kWh.takeOff  = BatteryEnergy_kWh.taxi-BattEnergy_kWh.taxi;
FC_power_kW.takeOff= 0.75*maxFC_Power_kW;
Weight_N.takeOff    = Weight_N.taxi-FuelCons_kg.taxi*9.81; 
FuelCons_kg.takeOff = Mission.time.takeOff_s*FC_specificFuelCons_kgkWh*FC_power_kW.takeOff/3600;
FCEnergy_kWh.takeOff= FC_power_kW.takeOff*Mission.time.takeOff_s/3600;

ElectricPower_kW.takeOff   = ShaftPower_kW.takeOff;
BattPower_kW.takeOff   = ElectricPower_kW.takeOff/electricMotorEff-FC_power_kW.takeOff;
BattEnergy_kWh.takeOff = (ElectricPower_kW.takeOff*Mission.time.takeOff_s/3600/electricMotorEff-FCEnergy_kWh.takeOff)/batteryEff;

% Climb
BatteryEnergy_kWh.climb  = BatteryEnergy_kWh.takeOff-BattEnergy_kWh.takeOff;
FC_power_kW.climb  = 0.75*maxFC_Power_kW;
Weight_N.climb    = Weight_N.takeOff-FuelCons_kg.takeOff*9.81;
FuelCons_kg.climb = Mission.time.climb_s*FC_specificFuelCons_kgkWh*FC_power_kW.climb/3600;
FCEnergy_kWh.climb   = FC_power_kW.climb*Mission.time.climb_s/3600;

ElectricPower_kW.climb     = ShaftPower_kW.climb;
BattPower_kW.climb   = ElectricPower_kW.climb/electricMotorEff-FC_power_kW.climb;
BattEnergy_kWh.climb   = (ElectricPower_kW.climb*Mission.time.climb_s/3600/electricMotorEff-FCEnergy_kWh.climb)/batteryEff;
% Cruise
BatteryEnergy_kWh.cruise = BatteryEnergy_kWh.climb-BattEnergy_kWh.climb;
Weight_N.cruise   = Weight_N.climb-FuelCons_kg.climb*9.81;
[FuelCons_kg.cruise,FC_power_kW.cruise,Weight_N.descent,...
    ElectricPower_kW.cruise,BattPower_kW.cruise,FCEnergy_kWh.cruise,...
    BattEnergy_kWh.cruise] = calculate_cruisePerformanceHydrogen(Mission,...
    MissionPerformance,ACFT,Weight_N.cruise,BatteryEnergy_kWh.cruise,hybridCruise);

% Descent
BatteryEnergy_kWh.descent = BatteryEnergy_kWh.cruise-BattEnergy_kWh.cruise;
FC_power_kW.descent = ShaftPower_kW.descent/electricMotorEff;
ElectricPower_kW.descent = 0;
FuelCons_kg.descent = FC_power_kW.descent*FC_specificFuelCons_kgkWh*Mission.time.descent_s/3600;
FCEnergy_kWh.descent  = FC_power_kW.descent*Mission.time.descent_s/3600;

ElectricPower_kW.descent   = ShaftPower_kW.descent;
BattPower_kW.descent       = 0;
BattEnergy_kWh.descent = (BattPower_kW.descent*Mission.time.descent_s/3600)/batteryEff; 

% Reserve Climb
BatteryEnergy_kWh.reserveClimb = BatteryEnergy_kWh.descent-BattEnergy_kWh.descent;
FC_power_kW.reserveClimb  = 0.75*maxFC_Power_kW;
Weight_N.reserveClimb      = Weight_N.descent-FuelCons_kg.descent*9.81;  
FuelCons_kg.reserveClimb = Mission.time.reserveClimb_s*FC_specificFuelCons_kgkWh*FC_power_kW.reserveClimb/3600;
FCEnergy_kWh.reserveClimb   = FC_power_kW.reserveClimb*Mission.time.reserveClimb_s/3600;

ElectricPower_kW.reserveClimb     = ShaftPower_kW.climb;
BattPower_kW.reserveClimb         = ElectricPower_kW.reserveClimb/electricMotorEff-FC_power_kW.reserveClimb;
BattEnergy_kWh.reserveClimb   = (ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s/3600/electricMotorEff-FCEnergy_kWh.reserveClimb)/batteryEff;

% Reserve Loiter
BatteryEnergy_kWh.reserveLoiter   = BatteryEnergy_kWh.reserveClimb-BattEnergy_kWh.reserveClimb;
Weight_N.reserveLoiter   = Weight_N.reserveClimb-FuelCons_kg.reserveClimb*9.81;
% Loiter function assumes that no battery energy will be used during this segment
[FuelCons_kg.reserveLoiter,FC_power_kW.reserveLoiter,...
    Weight_N.reserveDescent,ElectricPower_kW.reserveLoiter,...
    BattPower_kW.reserveLoiter,FCEnergy_kWh.reserveLoiter,...
    BattEnergy_kWh.reserveLoiter] = ...
    calculate_reserveLoiterPerformanceHydrogen(Mission,MissionPerformance,...
    ACFT,Weight_N.reserveLoiter);

% Reserve Descent
BatteryEnergy_kWh.reserveDescent = BatteryEnergy_kWh.reserveLoiter-BattEnergy_kWh.reserveLoiter;
FC_power_kW.reserveDescent = ShaftPower_kW.descent/electricMotorEff;
FuelCons_kg.reserveDescent = FC_power_kW.reserveDescent*FC_specificFuelCons_kgkWh*Mission.time.reserveDescent_s/3600;
FCEnergy_kWh.reserveDescent  = FC_power_kW.reserveDescent*Mission.time.reserveDescent_s/3600;

ElectricPower_kW.reserveDescent = ShaftPower_kW.descent;
BattPower_kW.reserveDescent     = ElectricPower_kW.reserveDescent/electricMotorEff-FC_power_kW.reserveDescent;
BattEnergy_kWh.reserveDescent = (ElectricPower_kW.reserveDescent*Mission.time.reserveDescent_s/3600/electricMotorEff-FCEnergy_kWh.reserveDescent)/batteryEff; 

% Average reserve power
FC_power_kW.reserve  = (FC_power_kW.reserveClimb*Mission.time.reserveClimb_s+FC_power_kW.reserveLoiter*Mission.time.reserveLoiter_s+FC_power_kW.reserveDescent*Mission.time.reserveDescent_s)...
    /(Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);
ElectricPower_kW.reserve = (ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s+ElectricPower_kW.reserveLoiter*Mission.time.reserveLoiter_s+ElectricPower_kW.reserveDescent*Mission.time.reserveDescent_s)...
    /(Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);

fuelConsMission_kg = FuelCons_kg.taxi+FuelCons_kg.takeOff+FuelCons_kg.climb...
    +FuelCons_kg.cruise+FuelCons_kg.descent;
fuelConsReserve_kg = FuelCons_kg.reserveClimb+FuelCons_kg.reserveLoiter+FuelCons_kg.reserveDescent;

battEnergyMission_kWh = BattEnergy_kWh.taxi + BattEnergy_kWh.takeOff + BattEnergy_kWh.climb + BattEnergy_kWh.cruise + BattEnergy_kWh.descent;
battEnergyReserve_kWh = BattEnergy_kWh.reserveClimb + BattEnergy_kWh.reserveLoiter + BattEnergy_kWh.reserveDescent;
BattPower_kW.reserve = battEnergyReserve_kWh/(time.reserveClimb_s+time.reserveLoiter_s+time.reserveDescent_s);
FCEnergyMission_kWh  = FCEnergy_kWh.taxi + FCEnergy_kWh.takeOff  + FCEnergy_kWh.climb + FCEnergy_kWh.cruise;
FCEnergyReserve_kWh  = FCEnergy_kWh.reserveClimb + FCEnergy_kWh.reserveLoiter+FCEnergy_kWh.reserveDescent;

% Change fuel consumption from kg to US gal
FuelCons_USgal.taxi    = FuelCons_kg.taxi/rhoFuel_kgm3*1000/3.785;
FuelCons_USgal.takeOff = FuelCons_kg.takeOff/rhoFuel_kgm3*1000/3.785;
FuelCons_USgal.climb   = FuelCons_kg.climb/rhoFuel_kgm3*1000/3.785;
FuelCons_USgal.cruise  = FuelCons_kg.cruise/rhoFuel_kgm3*1000/3.785;
FuelCons_USgal.descent = FuelCons_kg.descent/rhoFuel_kgm3*1000/3.785;
% Pack into a structure
Hydrogen.FC_power_kW   = FC_power_kW;
Hydrogen.weight_N       = Weight_N;
Hydrogen.FuelCons_USgal = FuelCons_USgal;
Hydrogen.FuelCons_kg    = FuelCons_kg; 
Hydrogen.fuelConsMission_kg = fuelConsMission_kg;
Hydrogen.fuelConsMission_USgal=fuelConsMission_kg/rhoFuel_kgm3*1000/3.785;
Hydrogen.fuelConsReserve_kg = fuelConsReserve_kg;
Hydrogen.fuelConsReserve_USgal = fuelConsReserve_kg/rhoFuel_kgm3*1000/3.785;
Hydrogen.fuelConsTotal_kg   = fuelConsMission_kg+fuelConsReserve_kg;

Hydrogen.ElectricPower_kW     = ElectricPower_kW;
Hydrogen.FCEnergy_kWh         = FCEnergy_kWh;
Hydrogen.BattPower_kW         = BattPower_kW;
Hydrogen.BattEnergy_kWh       = BattEnergy_kWh;
Hydrogen.batteryEnergy_kWh     = BatteryEnergy_kWh;
Hydrogen.battEnergyMission_kWh = battEnergyMission_kWh;
Hydrogen.battEnergyReserve_kWh = battEnergyReserve_kWh;
Hydrogen.FCEnergyMission_kWh  = FCEnergyMission_kWh;
Hydrogen.FCEnergyReserve_kWh  = FCEnergyReserve_kWh;


end

function [fuelCruise_kg,avgFCPower_kW,weightFinal_N,avgElecPower_kW,...
avgBattPower_kW,FCEnergy_kWh,BattEnergy_kWh] = ...
 calculate_cruisePerformanceHydrogen(Mission,MissionPerformance,ACFT,...
weight_N,initialBattEnergy_kWh,hybridCruise)
% Function:
%   calculate_missionPerformaceHydrogen 
%
% Description: 
%   Calculate the fuel and electric energy consumption for a particular mission
%   profile for a hyfrogen fuel cell retrofit
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N

% Output:
%   Hydrogen            - Structure with fuel, power and energy consumption
%


% Unwrap the mission and mission performance structures

timeCruise_s     = Mission.time.cruise_s;
altitudeCruise_m = Mission.altitude.cruise_m;
cruiseSpeed_ms   = MissionPerformance.cruiseSpeed_ms;
distanceCruise_m = Mission.distance.cruise_m;
wingArea_m2      = ACFT.wingArea_m2;
%hydrogenTankDrag = ACFT.Aero.hydrogenTankDrag;
cd0              = ACFT.Aero.cd0;%+hydrogenTankDrag;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Hydrogen.propEff;
FC_specificFuelCons_kgkWh = ACFT.Propulsion.Hydrogen.FC_specificFuelCons_kgkWh;
maxFC_Power_kW            = ACFT.Propulsion.Hydrogen.maxFC_Power_kW;
maxElectricPower_kW        = ACFT.Propulsion.Hydrogen.maxElectricPower_kW;
electricMotorEff           = ACFT.Propulsion.Hydrogen.electricMotorEff;
batteryEff                 = ACFT.Propulsion.Hydrogen.batteryEff;
rhoFuel_kgm3               = ACFT.Propulsion.Hydrogen.rhoFuel_kgm3;
batteryEnergy_kWh          = ACFT.EnergyStorage.Hydrogen.batteryEnergy_kWh;


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
fuel_kg       = ones(1,n); 
fuel_N        = ones(1,n);
powerFC_kW    = ones(1,n);      % Power delivered as electric power to the electric motor (needs multiplying by prop eff and motor eff)
powerElec_kW  = ones(1,n);      % Power delivered at shaft by electric motors (needs multiplying by prop eff to useful thrust power)
FCEnergy_kWh = ones(1,n);       % Energy supplied by fuel cell

if hybridCruise == 1
            reserveClimbEnergy_kWh = (maxElectricPower_kW/electricMotorEff-...
                0.75*maxFC_Power_kW)*Mission.time.reserveClimb_s/...
                3600/batteryEff;
            energyLeftOver_kWh=initialBattEnergy_kWh-0.05*batteryEnergy_kWh-...
                reserveClimbEnergy_kWh;
            if energyLeftOver_kWh > 0
                powerBatt_kW = (energyLeftOver_kWh)*batteryEff*...
                    electricMotorEff*3600/timeCruise_s; % Power delivered as electric power to the electric motor by the battery (needs multiplying by motor Eff)
            else
                disp('CAUTION: Series architecture does not have enough battery')
                powerBatt_kW = 0;
            end
else 
    powerBatt_kW = 0;
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
        powerElec_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;          % Electric Motor shaft power
        powerFC_kW(i)     = (powerElec_kW(i)-powerBatt_kW)/electricMotorEff;
        FCEnergy_kWh(i) = powerFC_kW(i)*timeSegment_s(i)/3600;
        fuel_N(i) = powerFC_kW(i)*FC_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
        fuel_kg(i) = fuel_N(i)/9.81;
    end
    
fuelCruise_N      = sum(fuel_N);
fuelCruise_kg  = sum(fuel_kg);
weightFinal_N     = weight_N-fuelCruise_N;
avgFCPower_kW    = sum(powerFC_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
avgBattPower_kW   = powerBatt_kW;
FCEnergy_kWh     = sum(FCEnergy_kWh);
BattEnergy_kWh    = avgBattPower_kW/batteryEff*timeCruise_s/3600;
end

function [fuelCruise_kg,avgFCPower_kW,weightFinal_N,avgElecPower_kW,...
avgBattPower_kW,FCEnergy_kWh,BattEnergy_kWh] = ...
 calculate_reserveLoiterPerformanceHydrogen(Mission,MissionPerformance,...
ACFT,weight_N,initialBattEnergy_kWh)
% Function:
%   calculate_reserveLoiterHydrogen
%
% Description: 
%   Calculate the fuel and electric energy consumption for the reserve
%   loiter segment for a hydrogen fuel cell aircraft

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
%   Hydrogen            - Structure with fuel, power and energy consumption
%


% Unwrap the mission and mission performance structures
timeCruise_s     = Mission.time.reserveLoiter_s;
altitudeCruise_m = Mission.altitude.reserve_m;
cruiseSpeed_ms   = MissionPerformance.reserveSpeed_ms;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Hydrogen.propEff;
FC_specificFuelCons_kgkWh = ACFT.Propulsion.Hydrogen.FC_specificFuelCons_kgkWh;
maxFC_Power_kW            = ACFT.Propulsion.Hydrogen.maxFC_Power_kW;
maxElectricPower_kW        = ACFT.Propulsion.Hydrogen.maxElectricPower_kW;
electricMotorEff           = ACFT.Propulsion.Hydrogen.electricMotorEff;
batteryEff                 = ACFT.Propulsion.Hydrogen.batteryEff;
rhoFuel_kgm3               = ACFT.Propulsion.Hydrogen.rhoFuel_kgm3;
batteryEnergy_kWh          = ACFT.EnergyStorage.Hydrogen.batteryEnergy_kWh;

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
fuel_kg       = ones(1,n); 
fuel_N        = ones(1,n);
powerFC_kW    = ones(1,n);      % Power delivered at the shaft by ICE engine (needs multiplying by prop eff)
powerElec_kW  = ones(1,n);      % Power delivered to electric motor (needs multiplying by motor eff and prop eff to get electric motor shaft power)
powerBatt_kW  = ones(1,n);
FCEnergy_kWh  = ones(1,n);
    for i = 1:n
        if i ~= 1
            weightSegment_N(i) = weightSegment_N(i-1)-fuel_N(i-1);
        end
        lift_N(i) = weightSegment_N(i);
        cl(i)     = lift_N(i)/(0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2);
        cd(i)     = cd0 + k*cl(i)^2;
        LoD(i)    = cl(i)/cd(i);
        drag_N(i) = 0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2*cd(i);
        powerElec_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;          % Electric Motor shaft power
        powerFC_kW(i)= powerElec_kW(i)/electricMotorEff;
        powerBatt_kW(i)   = powerElec_kW(i)/electricMotorEff-powerFC_kW(i);
        FCEnergy_kWh(i) = powerFC_kW(i)*timeSegment_s(i)/3600;
        fuel_N(i) = powerFC_kW(i)*FC_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
        fuel_kg(i) = fuel_N(i)/9.81;
   end
    
fuelCruise_N      = sum(fuel_N);
fuelCruise_kg     = sum(fuel_kg);
weightFinal_N     = weight_N-fuelCruise_N;
avgFCPower_kW     = sum(powerFC_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
avgBattPower_kW   = sum(powerBatt_kW.*timeSegment_s)/timeCruise_s;
FCEnergy_kWh     = sum(FCEnergy_kWh);
BattEnergy_kWh    = avgBattPower_kW/batteryEff*timeCruise_s;
end
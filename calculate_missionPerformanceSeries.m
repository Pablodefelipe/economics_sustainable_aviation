function Series = calculate_missionPerformanceSeries(Mission,...
    MissionPerformance,ACFT,weight_N,hybridCruise,rechargeSeries)
% Function:
%   calculate_missionPerformaceSeries
%
% Description: 
%   Calculate the fuel and electric energy consumption for a particular mission
%   profile for a series retrofit with and without hybrid cruise
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
%   hybridCruise        - 1 or 0 to use battery power during cruise or not
%   rechargeSeries      - 1 or 0 to charge the battery during cruise or not
% Output:
%   Series              - Structure with fuel, power and energy consumption
%

% Unwrap from Structures required variables
maxPower_kW                     = ACFT.Propulsion.maxPower_kW;
electricMotorEff                = ACFT.Propulsion.Series.electricMotorEff;
batteryEff                      = ACFT.Propulsion.Series.batteryEff;
generatorEff                    = ACFT.Propulsion.Series.generatorEff;
maxICE_Power_kW                 = ACFT.Propulsion.Series.maxICE_Power_kW;
maxElectricPower_kW             = ACFT.Propulsion.Series.maxElectricPower_kW;
ICE_maxSpecificFuelCons_kgkWh   = ACFT.Propulsion.Series.ICE_maxSpecificFuelCons_kgkWh;
ICE_specificFuelCons_kgkWh      = ACFT.Propulsion.Series.ICE_specificFuelCons_kgkWh;
rhoFuel_kgm3                    = ACFT.Propulsion.Series.rhoFuel_kgm3;
ShaftPower_kW                   = MissionPerformance.ShaftPower_kW;
%% Series
% Taxi
BatteryEnergy_kWh.taxi     = ACFT.EnergyStorage.Series.batteryEnergy_kWh; % Remaining Battery energy at start of segment
ICE_power_kW.taxi= 0;
GenPower_kW.taxi= ICE_power_kW.taxi*generatorEff;
Weight_N.taxi    = weight_N; 
FuelCons_USgal.taxi = 0;
GenEnergy_kWh.taxi=ICE_power_kW.taxi*Mission.time.taxi_s/3600*generatorEff;

ElectricPower_kW.taxi   = ShaftPower_kW.taxi;
BattPower_kW.taxi   = ElectricPower_kW.taxi/electricMotorEff-...
    GenPower_kW.taxi;
BattEnergy_kWh.taxi = (ElectricPower_kW.taxi*Mission.time.taxi_s/3600/electricMotorEff-GenEnergy_kWh.taxi)/batteryEff; 
% TO
BatteryEnergy_kWh.takeOff  = BatteryEnergy_kWh.taxi-BattEnergy_kWh.taxi;
ICE_power_kW.takeOff= maxICE_Power_kW;
GenPower_kW.takeOff = ICE_power_kW.takeOff*generatorEff;
Weight_N.takeOff    = Weight_N.taxi-FuelCons_USgal.taxi*3.785/1000*rhoFuel_kgm3*9.81; 
FuelCons_USgal.takeOff = Mission.time.takeOff_s*ICE_maxSpecificFuelCons_kgkWh*ICE_power_kW.takeOff/rhoFuel_kgm3*1000/3.785/3600;
GenEnergy_kWh.takeOff    = ICE_power_kW.takeOff*Mission.time.takeOff_s*generatorEff/3600;

ElectricPower_kW.takeOff   = maxElectricPower_kW;
BattPower_kW.takeOff   = ElectricPower_kW.takeOff/electricMotorEff-...
    GenPower_kW.takeOff;
BattEnergy_kWh.takeOff = (ElectricPower_kW.takeOff*Mission.time.takeOff_s/3600/electricMotorEff-GenEnergy_kWh.takeOff)/batteryEff;

% Climb
BatteryEnergy_kWh.climb  = BatteryEnergy_kWh.takeOff-BattEnergy_kWh.takeOff;
ICE_power_kW.climb  = maxICE_Power_kW;
GenPower_kW.climb = ICE_power_kW.climb*generatorEff;
Weight_N.climb    = Weight_N.takeOff-FuelCons_USgal.takeOff*3.785/1000*rhoFuel_kgm3*9.81;
FuelCons_USgal.climb = Mission.time.climb_s*ICE_maxSpecificFuelCons_kgkWh*ICE_power_kW.climb/rhoFuel_kgm3*1000/3.785/3600;
GenEnergy_kWh.climb   = ICE_power_kW.climb*Mission.time.climb_s*generatorEff/3600;

ElectricPower_kW.climb     = maxElectricPower_kW;
BattPower_kW.climb  = ElectricPower_kW.climb/electricMotorEff-...
    GenPower_kW.climb;
BattEnergy_kWh.climb   = (ElectricPower_kW.climb*Mission.time.climb_s/3600/electricMotorEff-GenEnergy_kWh.climb)/batteryEff;
% Cruise
BatteryEnergy_kWh.cruise = BatteryEnergy_kWh.climb-BattEnergy_kWh.climb;
Weight_N.cruise   = Weight_N.climb-FuelCons_USgal.climb*3.785/1000*rhoFuel_kgm3*9.81;
[FuelCons_USgal.cruise,ICE_power_kW.cruise,GenPower_kW.cruise,...
    Weight_N.descent,ElectricPower_kW.cruise,BattPower_kW.cruise,...
    GenEnergy_kWh.cruise,BattEnergy_kWh.cruise] = ...
    calculate_cruisePerformanceSeries(Mission,MissionPerformance,ACFT,...
    Weight_N.cruise,BatteryEnergy_kWh.cruise,hybridCruise,rechargeSeries);

% Descent
BatteryEnergy_kWh.descent = BatteryEnergy_kWh.cruise-BattEnergy_kWh.cruise;
ICE_power_kW.descent = ShaftPower_kW.descent/generatorEff/electricMotorEff;
GenPower_kW.descent = ICE_power_kW.descent*generatorEff;
FuelCons_USgal.descent = ICE_power_kW.descent*ICE_specificFuelCons_kgkWh*Mission.time.descent_s/rhoFuel_kgm3*1000/3.785/3600;
GenEnergy_kWh.descent  = ICE_power_kW.descent*generatorEff*Mission.time.descent_s/3600;

ElectricPower_kW.descent = ShaftPower_kW.descent;
BattPower_kW.descent  = ElectricPower_kW.descent/electricMotorEff-GenPower_kW.descent;
BattEnergy_kWh.descent = BattPower_kW.descent*Mission.time.descent_s...
    /3600/batteryEff;

% Reserve Climb
BatteryEnergy_kWh.reserveClimb = BatteryEnergy_kWh.descent-BattEnergy_kWh.descent;
ICE_power_kW.reserveClimb  = maxICE_Power_kW;
GenPower_kW.reserveClimb = ICE_power_kW.reserveClimb*generatorEff;
Weight_N.reserveClimb      = Weight_N.descent-FuelCons_USgal.descent*3.785/1000*rhoFuel_kgm3*9.81;  
FuelCons_USgal.reserveClimb = Mission.time.reserveClimb_s*ICE_maxSpecificFuelCons_kgkWh*ICE_power_kW.climb/rhoFuel_kgm3*1000/3.785/3600;
GenEnergy_kWh.reserveClimb   = ICE_power_kW.reserveClimb*Mission.time.reserveClimb_s*generatorEff/3600;

ElectricPower_kW.reserveClimb     = maxElectricPower_kW;
BattPower_kW.reserveClimb  = ElectricPower_kW.reserveClimb/electricMotorEff...
    -GenPower_kW.reserveClimb;
BattEnergy_kWh.reserveClimb   = (ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s/3600/electricMotorEff-GenEnergy_kWh.reserveClimb)/batteryEff;

% Reserve Loiter
BatteryEnergy_kWh.reserveLoiter   = BatteryEnergy_kWh.reserveClimb-BattEnergy_kWh.reserveClimb;
Weight_N.reserveLoiter   = Weight_N.reserveClimb-FuelCons_USgal.reserveClimb*3.785/1000*rhoFuel_kgm3*9.81;
% Loiter function assumes that no battery energy will be used during this segment
[FuelCons_USgal.reserveLoiter,ICE_power_kW.reserveLoiter,...
    Weight_N.reserveDescent,ElectricPower_kW.reserveLoiter,...
    BattEnergy_kWh.reserveLoiter] = ...
    calculate_reserveLoiterPerformanceSeries(Mission,MissionPerformance,...
    ACFT,Weight_N.reserveLoiter);
GenPower_kW.reserveLoiter = ICE_power_kW.reserveLoiter*generatorEff;
GenEnergy_kWh.reserveLoiter = ICE_power_kW.reserveLoiter*Mission.time.reserveLoiter_s/3600;
BattPower_kW.reserveLoiter  = BattEnergy_kWh.reserveLoiter*3600/...
    Mission.time.reserveLoiter_s;

% Reserve Descent
BatteryEnergy_kWh.reserveDescent = BatteryEnergy_kWh.reserveLoiter-BattEnergy_kWh.reserveLoiter;
ICE_power_kW.reserveDescent = ShaftPower_kW.descent/generatorEff/...
    electricMotorEff;
GenPower_kW.reserveDescent = ICE_power_kW.reserveDescent*generatorEff;
FuelCons_USgal.reserveDescent = ICE_power_kW.reserveDescent*ICE_specificFuelCons_kgkWh*Mission.time.reserveDescent_s/rhoFuel_kgm3*1000/3.785/3600;
GenEnergy_kWh.reserveDescent  = ICE_power_kW.reserveDescent*Mission.time.reserveDescent_s*generatorEff/3600;

ElectricPower_kW.reserveDescent =  ShaftPower_kW.descent;
BattPower_kW.reserveDescent  = ElectricPower_kW.reserveDescent/...
    electricMotorEff-GenPower_kW.reserveDescent;
BattEnergy_kWh.reserveDescent = 0;
% Average reserve power
ICE_power_kW.reserve=(ICE_power_kW.reserveClimb*Mission.time.reserveClimb_s...
    +ICE_power_kW.reserveLoiter*Mission.time.reserveLoiter_s)/...
    (Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+...
    Mission.time.reserveDescent_s);
ElectricPower_kW.reserve=(ElectricPower_kW.reserveClimb*...
    Mission.time.reserveClimb_s+ElectricPower_kW.reserveLoiter*...
    Mission.time.reserveLoiter_s)/(Mission.time.reserveClimb_s+...
    Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);

fuelConsMission_USgal = FuelCons_USgal.taxi+FuelCons_USgal.takeOff+FuelCons_USgal.climb...
    +FuelCons_USgal.cruise+FuelCons_USgal.descent;
fuelConsReserve_USgal = FuelCons_USgal.reserveClimb+...
    FuelCons_USgal.reserveLoiter+FuelCons_USgal.reserveDescent;

battEnergyMission_kWh = BattEnergy_kWh.taxi + BattEnergy_kWh.takeOff + ...
    BattEnergy_kWh.climb + BattEnergy_kWh.cruise + BattEnergy_kWh.descent;
battEnergyReserve_kWh = BattEnergy_kWh.reserveClimb +...
    BattEnergy_kWh.reserveLoiter + BattEnergy_kWh.reserveDescent;
genEnergyMission_kWh  = GenEnergy_kWh.taxi + GenEnergy_kWh.takeOff  + ...
    GenEnergy_kWh.climb + GenEnergy_kWh.cruise;
genEnergyReserve_kWh  = GenEnergy_kWh.reserveClimb + GenEnergy_kWh.reserveLoiter;
% Pack into a structure
Series.ICE_power_kW   = ICE_power_kW;
Series.weight_N       = Weight_N;
Series.fuelCons_USgal = FuelCons_USgal; 
Series.fuelConsMission_USgal = fuelConsMission_USgal;
Series.fuelConsReserve_USgal = fuelConsReserve_USgal;
Series.fuelConsMission_kg    = Series.fuelConsMission_USgal*0.003785*...
    rhoFuel_kgm3;
Series.fuelConsReserve_kg    = Series.fuelConsReserve_USgal*0.003785*...
    rhoFuel_kgm3;
Series.fuelConsTotal_kg = Series.fuelConsMission_kg+...
    Series.fuelConsReserve_kg;

Series.ElectricPower_kW     = ElectricPower_kW;
Series.GenPower_kW          = GenPower_kW;
Series.BattPower_kW         = BattPower_kW;
Series.GenEnergy_kWh        = GenEnergy_kWh;
Series.BattEnergy_kWh       = BattEnergy_kWh;
Series.batteryEnergy_kWh     = BatteryEnergy_kWh;
Series.battEnergyMission_kWh = battEnergyMission_kWh;
Series.battEnergyReserve_kWh = battEnergyReserve_kWh;
Series.genEnergyMission_kWh  = genEnergyMission_kWh;
Series.genEnergyReserve_kWh  = genEnergyReserve_kWh;
end

function [fuelCruise_USgal,avgICEPower_kW,avgGenPower_kW,weightFinal_N,...
    avgElecPower_kW,avgBatPower_kW,genEnergy_kWh,battEnergy_kWh] =  ...
    calculate_cruisePerformanceSeries(Mission,MissionPerformance,ACFT,...
    weight_N,initialBattEnergy_kWh,hybridCruise,rechargeSeries)
% Function:
%   calculate_missionPerformaceSeries
%
% Description: 
%   Calculate the fuel and electric energy consumption for a particular mission
%   profile for a series retrofit with and without charging
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
%   hybridCruise        - Option to choose whether electric energy is used
%                         during cruise or not
%   rechargeSeries      - Recharge battery during cruise option
% Output:
%   Series              - Structure with fuel, power and energy consumption
%

% Unwrap the mission and mission performance structures
timeCruise_s     = Mission.time.cruise_s;
altitudeCruise_m = Mission.altitude.cruise_m;
cruiseSpeed_ms   = MissionPerformance.cruiseSpeed_ms;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Series.propEff;
ICE_specificFuelCons_kgkWh = ACFT.Propulsion.Series.ICE_specificFuelCons_kgkWh;
maxICE_Power_kW            = ACFT.Propulsion.Series.maxICE_Power_kW;
maxElectricPower_kW        = ACFT.Propulsion.Series.maxElectricPower_kW;
electricMotorEff           = ACFT.Propulsion.Series.electricMotorEff;
generatorEff               = ACFT.Propulsion.Series.generatorEff;
batteryEff                 = ACFT.Propulsion.Series.batteryEff;
chargeEff                  = ACFT.Propulsion.Series.chargeEff;
rhoFuel_kgm3               = ACFT.Propulsion.Series.rhoFuel_kgm3;
batteryEnergy_kWh          = ACFT.EnergyStorage.Series.batteryEnergy_kWh;

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
powerICE_kW   = ones(1,n);      % Power delivered at the shaft by ICE engine (needs multiplying by prop eff)
powerElec_kW  = ones(1,n);      % Power delivered to electric motor (needs multiplying by motor eff and prop eff to get electric motor shaft power)
genEnergy_kWh = ones(1,n);      % Energy supplied by generator
electricEnergy_kWh = ones(1,n); % Energy supplied by battery
    for i = 1:n
        if i ~= 1
            weightSegment_N(i) = weightSegment_N(i-1)-fuel_N(i-1);
        end
        if hybridCruise == 1
            reserveClimbEnergy_kWh = (maxElectricPower_kW/electricMotorEff-...
                maxICE_Power_kW*generatorEff)*Mission.time.reserveClimb_s/...
                3600/batteryEff;
            energyLeftOver_kWh=initialBattEnergy_kWh-0.05*batteryEnergy_kWh-...
                reserveClimbEnergy_kWh;
            if energyLeftOver_kWh > 0
                battPower_kW = (energyLeftOver_kWh)*batteryEff*...
                    electricMotorEff*3600/timeCruise_s;
            else
                disp('CAUTION: Series architecture does not have enough battery')
                battPower_kW = 0;
            end
        end
        lift_N(i) = weightSegment_N(i);
        cl(i)     = lift_N(i)/(0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2);
        cd(i)     = cd0 + k*cl(i)^2;
        LoD(i)    = cl(i)/cd(i);
        drag_N(i) = 0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2*cd(i);
        powerElec_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;          % Electric Motor shaft power
        powerICE_kW(i)= (powerElec_kW(i)/electricMotorEff-battPower_kW)/generatorEff;
       
        if rechargeSeries == 1
            energyRequired_kWh = (0.9*batteryEnergy_kWh-initialBattEnergy_kWh)/chargeEff;    % Battery energy required to reach 90% capacity  
            surplusPowerICE_kW = energyRequired_kWh/(timeCruise_s/3600)*ones(1,n);
            powerICE_kW(i) = powerICE_kW(i)+surplusPowerICE_kW(i);
            if powerICE_kW(i) > maxICE_Power_kW
                surplusPowerICE_kW(i) = (maxICE_Power_kW-powerElec_kW(i)...
                    /electricMotorEff)/generatorEff;
                powerICE_kW(i)        = maxICE_Power_kW;
            end
            electricEnergy_kWh(i)     = -surplusPowerICE_kW/chargeEff*timeSegment_s(i)/3600;
        else
            electricEnergy_kWh(i)     = 0;
        end
        genEnergy_kWh(i) = powerICE_kW(i)*timeSegment_s(i)/3600*generatorEff;
        fuel_N(i) = powerICE_kW(i)*ICE_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
        fuel_USgal(i) = fuel_N(i)/9.81/rhoFuel_kgm3*1000/3.785;
    end
    
fuelCruise_N      = sum(fuel_N);
fuelCruise_USgal  = sum(fuel_USgal);
weightFinal_N     = weight_N-fuelCruise_N;
avgICEPower_kW    = sum(powerICE_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW
avgGenPower_kW    = avgICEPower_kW*generatorEff;
avgBatPower_kW    = battPower_kW;
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
genEnergy_kWh     = sum(genEnergy_kWh);
elecEnergy_kWh    = sum(electricEnergy_kWh);
battEnergy_kWh    = avgBatPower_kW*timeCruise_s/3600;
end

function [fuelCruise_USgal,avgICEPower_kW,weightFinal_N,avgElecPower_kW,...
    elecEnergy_kWh] =  calculate_reserveLoiterPerformanceSeries(...
    Mission,MissionPerformance,ACFT,weight_N)
% Function:
%   calculate_reserveLoiterSeries
%
% Description: 
%   Calculate the fuel and electric energy consumption for the reserve
%   loiter segment
%   profile for a series 

%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
% Output:
%   Series              - Structure with fuel, power and energy consumption
%

% Unwrap the mission and mission performance structures
timeCruise_s     = Mission.time.reserveLoiter_s;
altitudeCruise_m = Mission.altitude.reserve_m;
cruiseSpeed_ms   = MissionPerformance.reserveSpeed_ms;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Series.propEff;
ICE_specificFuelCons_kgkWh = ACFT.Propulsion.Series.ICE_specificFuelCons_kgkWh;
maxICE_Power_kW            = ACFT.Propulsion.Series.maxICE_Power_kW;
maxElectricPower_kW        = ACFT.Propulsion.Series.maxElectricPower_kW;
electricMotorEff           = ACFT.Propulsion.Series.electricMotorEff;
generatorEff               = ACFT.Propulsion.Series.generatorEff;
batteryEff                 = ACFT.Propulsion.Series.batteryEff;
chargeEff                  = ACFT.Propulsion.Series.chargeEff;
rhoFuel_kgm3               = ACFT.Propulsion.Series.rhoFuel_kgm3;
batteryEnergy_kWh          = ACFT.EnergyStorage.Series.batteryEnergy_kWh;


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
powerICE_kW   = ones(1,n);      % Power delivered at the shaft by ICE engine (needs multiplying by prop eff)
powerElec_kW  = ones(1,n);      % Power delivered to electric motor (needs multiplying by motor eff and prop eff to get electric motor shaft power)
electricEnergy_kWh = ones(1,n);
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
        powerICE_kW(i)= powerElec_kW(i)/electricMotorEff/generatorEff;
        electricEnergy_kWh(i)     = 0;
        
        fuel_N(i) = powerICE_kW(i)*ICE_specificFuelCons_kgkWh*timeSegment_s(i)/3600*9.81;
        fuel_USgal(i) = fuel_N(i)/9.81/rhoFuel_kgm3*1000/3.785;
   end
    
fuelCruise_N      = sum(fuel_N);
fuelCruise_USgal  = sum(fuel_USgal);
weightFinal_N     = weight_N-fuelCruise_N;
avgICEPower_kW    = sum(powerICE_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
elecEnergy_kWh    = sum(electricEnergy_kWh);
end
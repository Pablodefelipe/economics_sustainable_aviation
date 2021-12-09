function Electric = calculate_missionPerformanceElectric(Mission,MissionPerformance,ACFT,weight_N,rechargeSeries)
% Function:
%   calculate_missionPerformaceElectric
%
% Description: 
%   Calculate the electric energy consumption for a particular mission
%   profile for an electric retrofit
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
% Output:
%   Electric            - Structure with power and energy consumption
%

% Unwrap from Structures required variables
electricMotorEff                = ACFT.Propulsion.Electric.electricMotorEff;
batteryEff                      = ACFT.Propulsion.Electric.batteryEff;
maxElectricPower_kW             = ACFT.Propulsion.Electric.maxElectricPower_kW;
ShaftPower_kW                   = MissionPerformance.ShaftPower_kW;
%% Electric
% Taxi
BatteryEnergy_kWh.taxi          = ACFT.EnergyStorage.Electric.batteryEnergy_kWh; % Remaining Battery energy at start of segment

ElectricPower_kW.taxi   = ShaftPower_kW.taxi;
BattEnergy_kWh.taxi = (ElectricPower_kW.taxi*Mission.time.taxi_s/3600)/batteryEff/electricMotorEff; 
% TO
BatteryEnergy_kWh.takeOff  = BatteryEnergy_kWh.taxi-BattEnergy_kWh.taxi;

ElectricPower_kW.takeOff   = maxElectricPower_kW;
BattEnergy_kWh.takeOff = (ElectricPower_kW.takeOff*Mission.time.takeOff_s/3600/electricMotorEff)/batteryEff;

% Climb
BatteryEnergy_kWh.climb  = BatteryEnergy_kWh.takeOff-BattEnergy_kWh.takeOff;
ElectricPower_kW.climb     = maxElectricPower_kW;
BattEnergy_kWh.climb   = (ElectricPower_kW.climb*Mission.time.climb_s/3600/electricMotorEff)/batteryEff;
% Cruise
BatteryEnergy_kWh.cruise = BatteryEnergy_kWh.climb-BattEnergy_kWh.climb;
[ElectricPower_kW.cruise,BattEnergy_kWh.cruise] = calculate_cruisePerformanceElectric(Mission,MissionPerformance,ACFT,weight_N);

% Descent
BatteryEnergy_kWh.descent = BatteryEnergy_kWh.cruise-BattEnergy_kWh.cruise;
ElectricPower_kW.descent = ShaftPower_kW.descent;
BattEnergy_kWh.descent = ElectricPower_kW.descent*Mission.time.descent_s/3600;

% Reserve Climb
BatteryEnergy_kWh.reserveClimb = BatteryEnergy_kWh.descent-BattEnergy_kWh.descent;
ElectricPower_kW.reserveClimb  = maxElectricPower_kW;
BattEnergy_kWh.reserveClimb    = (ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s/3600/electricMotorEff)/batteryEff;

% Reserve Loiter
BatteryEnergy_kWh.reserveLoiter   = BatteryEnergy_kWh.reserveClimb-BattEnergy_kWh.reserveClimb;
[ElectricPower_kW.reserveLoiter,BattEnergy_kWh.reserveLoiter] = calculate_reserveLoiterPerformanceElectric(Mission,MissionPerformance,ACFT,weight_N);
% Reserve Descent
BatteryEnergy_kWh.reserveDescent = BatteryEnergy_kWh.reserveLoiter-BattEnergy_kWh.reserveLoiter;

ElectricPower_kW.reserveDescent = ShaftPower_kW.descent;
BattEnergy_kWh.reserveDescent = ElectricPower_kW.reserveDescent*Mission.time.reserveDescent_s/electricMotorEff/batteryEff/3600;
% Average reserve power
ElectricPower_kW.reserve = (ElectricPower_kW.reserveClimb*Mission.time.reserveClimb_s+ElectricPower_kW.reserveLoiter*Mission.time.reserveLoiter_s)...
    /(Mission.time.reserveClimb_s+Mission.time.reserveLoiter_s+Mission.time.reserveDescent_s);

battEnergyMission_kWh = BattEnergy_kWh.taxi + BattEnergy_kWh.takeOff + BattEnergy_kWh.climb + BattEnergy_kWh.cruise + BattEnergy_kWh.descent;
battEnergyReserve_kWh = BattEnergy_kWh.reserveClimb + BattEnergy_kWh.reserveLoiter + BattEnergy_kWh.reserveDescent;

% Pack into a structure
Electric.ElectricPower_kW     = ElectricPower_kW;
Electric.BattEnergy_kWh       = BattEnergy_kWh;     
Electric.batteryEnergy_kWh     = BatteryEnergy_kWh;  
Electric.battEnergyMission_kWh = battEnergyMission_kWh;
Electric.battEnergyReserve_kWh = battEnergyReserve_kWh;
Electric.battEnergyTotal_kWh   = battEnergyMission_kWh+battEnergyReserve_kWh;
end

function [avgElecPower_kW,BattEnergy_kWh] =  calculate_cruisePerformanceElectric(Mission,MissionPerformance,ACFT,weight_N)
% Function:
%   calculate_cruisePerformaceElectric 
%
% Description: 
%   Calculate the electric energy consumption during for an electric retrofit
%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N

% Output:
%   avgElecPower_kW     - Average power delivered at the shaft of the
%                         electric motor
%   BattEnergy_kWh      - Battery energy consumed during the segment


% Unwrap the mission and mission performance structures

timeCruise_s     = Mission.time.cruise_s;
altitudeCruise_m = Mission.altitude.cruise_m;
cruiseSpeed_ms   = MissionPerformance.cruiseSpeed_ms;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Electric.propEff;
electricMotorEff = ACFT.Propulsion.Electric.electricMotorEff;
batteryEff       = ACFT.Propulsion.Electric.batteryEff;


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
powerElec_kW  = ones(1,n);      % Power delivered at shaft by electric motors (needs multiplying by prop eff to useful thrust power)
powerBatt_kW  = ones(1,n);      % Power delivered as electric power to the electric motor by the battery (needs multiplying by motor Eff)

    for i = 1:n
        if i ~= 1
            weightSegment_N(i) = weightSegment_N(i-1);
        end
        lift_N(i) = weightSegment_N(i);
        cl(i)     = lift_N(i)/(0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2);
        cd(i)     = cd0 + k*cl(i)^2;
        LoD(i)    = cl(i)/cd(i);
        drag_N(i) = 0.5*rho_kgm3*cruiseSpeed_ms^2*wingArea_m2*cd(i);
        powerElec_kW(i)   = drag_N(i)*cruiseSpeed_ms/propEff/1000;          % Electric Motor shaft power
        powerBatt_kW(i)   = powerElec_kW(i)/electricMotorEff;
    end
    
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW
avgBattPower_kW   = sum(powerBatt_kW.*timeSegment_s)/timeCruise_s;
BattEnergy_kWh    = avgBattPower_kW/batteryEff*timeCruise_s/3600;
end

function [avgElecPower_kW,elecEnergy_kWh] =  calculate_reserveLoiterPerformanceElectric(Mission,MissionPerformance,ACFT,weight_N);
% Function:
%   calculate_reserveLoiterElectric
%
% Description: 
%   Calculate the electric energy consumption for the reserve
%   loiter segment for an electric retrofit

%
% Input:
%   Mission             - Definition of mission
%   MissionPerformance  - Performance of aircraft defined for all the
%                           missions
%   ACFT                - Contains the data on the aircraft
%   weight_N            - Starting weight in N
%
% Output:
%   avgElecPower_kW     - shaft power output of electric engines
%   elecEnergy_kWh      - Battery energy consumed


% Unwrap the mission and mission performance structures
timeCruise_s     = Mission.time.reserveLoiter_s;
altitudeCruise_m = Mission.altitude.reserve_m;
cruiseSpeed_ms   = MissionPerformance.reserveSpeed_ms;
wingArea_m2      = ACFT.wingArea_m2;
cd0              = ACFT.Aero.cd0;
k                = ACFT.Aero.k;

propEff          = ACFT.Propulsion.Electric.propEff;
maxElectricPower_kW        = ACFT.Propulsion.Electric.maxElectricPower_kW;
electricMotorEff           = ACFT.Propulsion.Electric.electricMotorEff;
batteryEff                 = ACFT.Propulsion.Electric.batteryEff;


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
        electricEnergy_kWh(i)     = powerElec_kW(i)*timeSegment_s(i)/3600/electricMotorEff/batteryEff;
   end
    
avgElecPower_kW   = sum(powerElec_kW.*timeSegment_s)/timeCruise_s; % Average Power in kW 
elecEnergy_kWh    = sum(electricEnergy_kWh);
end

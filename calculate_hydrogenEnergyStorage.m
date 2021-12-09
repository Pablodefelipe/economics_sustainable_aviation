function [batteryEnergy_kWh,fuelTank_kg] = calculate_hydrogenEnergyStorage(ACFT,Hydrogen)
% Function:
%   calculate_HydrogenBattery
% Description: 
%   Calculates the battery energy which the hydrogen retrofit can carry
%   given the tank weight depends on the hydrogen volume.
% Input:
%   ACFT         - Structure with aircarft data
%   Hydrogen     - Structure with Hydrogen data -> Need available energy
%                  storage mass
%  Parameters unwrapped from structures:
%   rhoFuel_kgm3 - LH2 fuel density
%   LH2SpecificEnergy_kWhkg - Specific energy of hydrogen fuel
%   tankNumber              - Number of cylindrical tanks
%   hRatioEnergy            - Hybridisation ratio of energy
% Output:
%   AvailableWeight_kg - Available weight for fuel tank + hydrogen

hRatioEnergy = ACFT.Propulsion.hybridisationRatioEnergy;   
Weight       = ACFT.Weight;
LH2SpecificEnergy_kWhkg = ACFT.EnergyStorage.Hydrogen.specificEnergy_kWhkg;
tankNumber              = ACFT.EnergyStorage.Hydrogen.tankNumber;
slendernessRatio        = ACFT.EnergyStorage.Hydrogen.slendernessRatio;
rhoFuel_kgm3            = ACFT.Propulsion.Hydrogen.rhoFuel_kgm3;
% Unwrap the structure
energyStorage_kg                  = Hydrogen.energyStorage_kg;
hydrogenFuelTankAreaDensity_kgm2  = Weight.hydrogenFuelTankAreaDensity_kgm2;
batterySpecificEnergy_kWhkg       = Weight.batterySpecificEnergy_kWhkg;

% Solve the system A*batteryEnergy_kWh + B*batteryEnergy_kWh^(2/3) =
% energyStorage_kg;

% Compute A and B values
A = 1/batterySpecificEnergy_kWhkg+(1/hRatioEnergy-1)/...
    LH2SpecificEnergy_kWhkg;
B = hydrogenFuelTankAreaDensity_kgm2*tankNumber^(1/3)*slendernessRatio*...
    pi*(12/(pi*(3*slendernessRatio-1)))^(2/3)*(1/(rhoFuel_kgm3*...
    LH2SpecificEnergy_kWhkg)*(1/hRatioEnergy-1))^(2/3);
options = optimset('Display','off');
fun = @(batteryEnergy_kWh)[A*batteryEnergy_kWh+B*batteryEnergy_kWh^(2/3)...
    - energyStorage_kg];
x0  = 39;
batteryEnergy_kWh = fsolve(fun,x0,options);
hydrogenEnergy_kWh= (1/hRatioEnergy-1)*batteryEnergy_kWh;
hydrogenWeight_kg = hydrogenEnergy_kWh/LH2SpecificEnergy_kWhkg;
hydrogenVolume_m3 = hydrogenWeight_kg/rhoFuel_kgm3;
fuelTank_kg       = slendernessRatio*pi*tankNumber^(1/3)*(12/(pi*...
    (3*slendernessRatio-1))*hydrogenVolume_m3)^(2/3)*...
    hydrogenFuelTankAreaDensity_kgm2;
end






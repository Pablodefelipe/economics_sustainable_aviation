function Weight = calculate_baselineEmptyWeight(ACFT)

% Function:
%   calculate_baselineEmptyWeight(ACFT)
%
% Description: 
%   Calculates the empty weight of the aircraft without propulsion and fuel
%   system
%   Also calculates the fuelTankSystem weight for the baseline, parallel
%   and series hybrids
%
% Input:
%   ACFT                - In particular need EW, engine specific power,
%                         power and fuel quantity
% Output:
%   Weight              - Baseline Empty Weight in kg

% Copy the weight structure
Weight = ACFT.Weight;
% Unwrap the structure ACFT
EW_kg = ACFT.Weight.EW_kg;
icePistonSpecificPower_kWkg = ACFT.Weight.icePistonSpecificPower_kWkg;
maxPower_kW = ACFT.Propulsion.maxPower_kW;
maxFuelVolume_m3 = ACFT.EnergyStorage.Baseline.maxFuelVolume_m3;
engineNum        = ACFT.engineNum;
% Calculate the engine and tank weight
engineWeight_kg = maxPower_kW/icePistonSpecificPower_kWkg;
fuelTankSystem_kg= calculate_fuelTankWeight(maxFuelVolume_m3,engineNum);

% Calculate the BEW
baselineEmpty_kg = EW_kg - engineWeight_kg - fuelTankSystem_kg;

Weight.baselineEmpty_kg = baselineEmpty_kg;
Weight.fuelTankSystem_kg = fuelTankSystem_kg;

end

function weight_kg = calculate_fuelTankWeight(fuelVolume_m3,engineNum)
% Function:
%   calculate_fuelTankWeight
%
% Description: 
%   Calculate fuel tank volume for the retrofits. Note the following
%   assumptions have been used in the equation:
%   1. All fuel is stored in integral tanks
%   2. Two tanks are used
%
% Input:
%   fuelVolume_m3       - Volume of fuel
%                         
% Output:
%   weight_kg           - Weight of fuel tank

% Convert volume into US Gal
fuelVolume_USgal = fuelVolume_m3/0.003785;

weight_lbs = 2.49*(fuelVolume_USgal)^0.726*(0.5)^0.363*2^0.242*engineNum^0.157;
weight_kg  = weight_lbs/2.205;

end
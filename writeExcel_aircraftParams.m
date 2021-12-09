function writeExcel_aircraftParams(ACFT)
% Function:
%   writeExcel_aircraftParams(ACFT)
%
% Description: 
%   Exports the aircraft parameter summary to excel
%
% Input:
%   ACFT                - Aircraft data

% Output:

% Unwrap structure
wingArea_m2  = ACFT.wingArea_m2 ;      % Reference wing area in m^2
aspectRatio  = ACFT.aspectRatio; 
MTOW_kg      = ACFT.Weight.MTOW_kg;  % MTOW of aircraft
EW_kg        = ACFT.Weight.EW_kg;  % Operational empty weight of aircraft, inluding unusable fuel, oil and full operating fluids
BEW_kg       = ACFT.Weight.baselineEmpty_kg;

% Aero parameters: cd0 and k
cd0          = ACFT.Aero.cd0;
k            = ACFT.Aero.k;
pax          = ACFT.Weight.passenger_n; 

filename    = 'exampleData.xlsx';
Parameter   = ["MTOW";"EW";"BEW";"Wing Area";"AR";"pax";"Cd0";"k"];
Value       = [MTOW_kg;EW_kg;BEW_kg;wingArea_m2;aspectRatio;pax;cd0;k];
Units       = ["kg";    "kg";  "kg";      "m^2";        "-";"-";"-";"-"];

params = table(Parameter,Value,Units);

writetable(params,filename,'Sheet',1,'Range','A4')
end

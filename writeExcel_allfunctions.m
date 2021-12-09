function writeExcel_allfunctions(ACFT,Mission,MissionPerformance)
% Function:
%   writeExcel_allfunctions(ACFT)

% Description: 
%   writes an excel with all the data required
%
% Input:
%   ACFT                - Aircraft data
%   Mission             - Mission data
%   MissionPerformnace  - MissionPerformance

% Output:

writeExcel_aircraftParams(ACFT);
writeExcel_mission(Mission,MissionPerformance);
writeExcel_propulsion(ACFT);
writeExcel_weight(ACFT);
end


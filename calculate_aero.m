function [acf] = calculate_aero(ACFT)
% Function to calculate the Lift over drag ratio for cruise or loiter
% Calculate L/D ratio for cruise
% Method 1 - Roskam's Drag polar method to estimate Cd0
% Cd = Cd0 + CL^2/(pi*AR*e)
if acf.type == "Twin"
    % For twin engine constants for empirical equations are:
    a = -2.2218;
    b = 1.0;
    c = 0.8635;
    d = 0.5632;
end
if acf.type == "Single"
    % For single engine constants for empirical equations are:
    a = -2.2218;
    b = 1.0;
    c = 1.0892;
    d = 0.5147;
end
if acf.type == "Turboprop"
    % For regional turboprop engine constants for empirical equations are:
    a = -2.2218;
    b = 1.0;
    c = -0.0866;
    d = 0.8099;
end
acf.Swet_ft2 = 10^(c+d*log10(acf.Wo_lb));
acf.f = 10^(a+b*log10(acf.Swet_ft2));
acf.Cd0 = acf.f/acf.S_ft2;              % Min. drag coefficient
acf.e  = 0.8;                           % Oswald efficiency
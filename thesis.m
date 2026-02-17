%% MD Validation using governing equations
clear; clc; close all;

%% ---------------- Constants ----------------
Hfg    = 2.406e6;       % J/kg (latent heat of vaporization)
Rconst = 461.5;         % J/kg·K (specific gas constant for water vapor)
Thi    = 66;            % °C, feed hot side
Tci    = 44;            % °C, cold side
Tm     = (Thi + Tci)/2; % mean temperature
dT     = Thi - Tci;     % temperature difference
Pm     = 1.013e5;       % Pa (mean pore pressure, ~1 atm)
gamma  = 0.067;         % N/m (surface tension of water)
tau    = 2;             % tortuosity (assumed)

%% ---------------- Parameter sweeps ----------------
th_um  = [70.31 80 90 110 130 160 190.18];    % thickness µm
epsil  = [0.50 0.60 0.70 0.80 0.90];          % porosity (-)
d_um   = [0.22 0.30 0.40 0.50 0.60 0.70 0.80];% pore diameter µm
theta  = [90 100 110 120 130 140 150];        % contact angle deg

% Baseline
eps0 = 0.80; d0 = 0.22; th0 = 70.31; theta0 = 120;

%% ---------------- Mass transfer coefficient function ----------------
% Overall mass transfer coefficient Cw ∝ (ε * d) / (τ * δ)
Cw_fun = @(eps, d_um_loc, th_um_loc) ...
    (eps .* d_um_loc*1e-6) ./ (tau .* th_um_loc*1e-6);

% Baseline scaling to match Jp = 66.79 kg/m^2·h
Cw_base   = Cw_fun(eps0, d0, th0);
scaleFact = 66.79/3600 / ( Cw_base * (Hfg/(Rconst*(Tm+273)*Pm)) * dT );

% Final calibrated Cw function
Cw_fun_cal = @(eps, d_um_loc, th_um_loc) ...
    scaleFact * ( (eps .* d_um_loc*1e-6) ./ (tau .* th_um_loc*1e-6) );

%% ---------------- Governing equations ----------------
% Permeate flux
Jp_fun = @(eps, d_um_loc, th_um_loc) ...
    Cw_fun_cal(eps, d_um_loc, th_um_loc) .* ...
    (Hfg ./ (Rconst*(Tm+273).*Pm)) .* dT;   % kg/m2·s
% Convert to kg/m2·h
Jp_fun_h = @(eps, d_um_loc, th_um_loc) 3600 * Jp_fun(eps, d_um_loc, th_um_loc);

% Mass transfer coefficient
Km_fun = @(Jp_h) ( (Jp_h/3600) * Hfg ) / dT / 1000; % kJ/m2·s·°C

% LEP
LEP_fun = @(d_um_loc, theta_deg) ...
    (2*gamma*abs(cosd(theta_deg))) ./ (d_um_loc*1e-6) / 1000; % kPa

%% ---------------- Model predictions ----------------
% Thickness sweep
Jp_th_mod  = Jp_fun_h(eps0, d0, th_um);
Km_th_mod  = Km_fun(Jp_th_mod);
LEP_th_mod = LEP_fun(d0, theta0) * ones(size(th_um));

% Porosity sweep
Jp_eps_mod  = Jp_fun_h(epsil, d0, th0);
Km_eps_mod  = Km_fun(Jp_eps_mod);
LEP_eps_mod = LEP_fun(d0, theta0) * ones(size(epsil));

% Pore size sweep
Jp_d_mod  = Jp_fun_h(eps0, d_um, th0);
Km_d_mod  = Km_fun(Jp_d_mod);
LEP_d_mod = LEP_fun(d_um, theta0);

% Contact angle sweep
Jp_tht_mod  = Jp_fun_h(eps0, d0, th0) * ones(size(theta));
Km_tht_mod  = Km_fun(Jp_tht_mod);
LEP_tht_mod = LEP_fun(d0, theta);

%% ---------------- Plotting ----------------
plotOv(th_um,  Jp_th_mod,  'Thickness (\mum)',    'Permeate flux (kg m^{-2} h^{-1})', 'Permeate flux vs Thickness (\mum)');
plotOv(epsil,  Jp_eps_mod, 'Porosity (-)',        'Permeate flux (kg m^{-2} h^{-1})', 'Permeate flux vs Porosity (-)');
plotOv(d_um,   Jp_d_mod,   'Pore diameter (\mum)','Permeate flux (kg m^{-2} h^{-1})', 'Permeate flux vs Pore diameter (\mum)');
plotOv(theta,  Jp_tht_mod, 'Contact angle (deg)', 'Permeate flux (kg m^{-2} h^{-1})', 'Permeate flux vs Contact angle (deg)');

plotOv(th_um,  Km_th_mod,  'Thickness (\mum)',    'Mass transfer coefficient (kJ m^{-2} s^{-1} ^\circC)', 'Mass transfer coefficient vs Thickness (\mum)');
plotOv(epsil,  Km_eps_mod, 'Porosity (-)',        'Mass transfer coefficient (kJ m^{-2} s^{-1} ^\circC)', 'Mass transfer coefficient vs Porosity (-)');
plotOv(d_um,   Km_d_mod,   'Pore diameter (\mum)','Mass transfer coefficient (kJ m^{-2} s^{-1} ^\circC)', 'Mass transfer coefficient vs Pore diameter (\mum)');
plotOv(theta,  Km_tht_mod, 'Contact angle (deg)', 'Mass transfer coefficient (kJ m^{-2} s^{-1} ^\circC)', 'Mass transfer coefficient vs Contact angle (deg)');

plotOv(th_um,  LEP_th_mod,  'Thickness (\mum)',    'LEP (kPa)', 'LEP vs Thickness (\mum)');
plotOv(epsil,  LEP_eps_mod, 'Porosity (-)',        'LEP (kPa)', 'LEP vs Porosity (-)');
plotOv(d_um,   LEP_d_mod,   'Pore diameter (\mum)','LEP (kPa)', 'LEP vs Pore diameter (\mum)');
plotOv(theta,  LEP_tht_mod, 'Contact angle (deg)', 'LEP (kPa)', 'LEP vs Contact angle (deg)');

%% ---------------- Helper function ----------------
function plotOv(x, y, xl, yl, ttl)
    figure('Color','w');
    plot(x, y, 'bo-', 'LineWidth',1.6, 'MarkerFaceColor','b');
    grid on;
    xlabel(xl, 'Interpreter','tex');
    ylabel(yl, 'Interpreter','tex');
    title(ttl, 'Interpreter','tex');
end

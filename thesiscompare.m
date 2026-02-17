%% ========= MD: Compare Optimized Design vs PDF Reference (units fixed) =========
clear; clc; close all;
%% ---- Physics constants ----
Hfg    = 2.406e6;          % J/kg (latent heat)
Rconst = 461.5;            % J/(kg·K) (water vapor)
Thi    = 66;               % °C (feed hot)
Tci    = 44;               % °C (permeate cold)
%% ========= MD: Compare Optimized Design vs Journal Article Data (units fixed) =========
clear; clc; close all;
%% ---- Physics constants ----
Hfg    = 2.406e6;          % J/kg (latent heat)
Rconst = 461.5;            % J/(kg·K) (water vapor)
Thi    = 66;               % °C (feed hot)
Tci    = 44;               % °C (permeate cold)
Tm     = (Thi + Tci)/2;    % °C (mean temperature)
dT     = Thi - Tci;        % °C (temperature difference)
Pm     = 1.013e5;          % Pa (mean pore pressure)
gamma  = 0.067;            % N/m (surface tension)
tau    = 2;                % (-) tortuosity
%% ---- Governing relations ----
Cw_raw = @(eps, d_um, th_um) (eps .* (d_um*1e-6)) ./ (tau .* (th_um*1e-6));
Jp_h   = @(Cw) 3600 .* Cw .* (Hfg./(Rconst*(Tm+273).*Pm)) .* dT;
Km_fun = @(Jp_h_val) ((Jp_h_val/3600) .* Hfg ./ dT) ./ 1000;  % kJ·m^-2·s^-1·°C^-1
LEP_fun = @(d_um, theta_deg) (2*gamma*abs(cosd(theta_deg))) ./ (d_um*1e-6) ./ 1000; % kPa
%% ---- Calibration ----
eps0 = 0.80; d0 = 0.22; th0 = 70.31; Jp_baseline_h = 66.79; % kg·m^-2·h^-1
Cw_base   = Cw_raw(eps0, d0, th0);
scaleFact = (Jp_baseline_h/3600) / ( Cw_base * (Hfg/(Rconst*(Tm+273)*Pm)) * dT );
Cw = @(eps, d_um, th_um) scaleFact .* Cw_raw(eps, d_um, th_um);
%% ====== 1) Optimized Design ======
eps_opt   = 0.9021;
d_opt     = 0.2496;
th_opt    = 76.5204;
theta_opt = 146.5;
Jp_opt  = Jp_h(Cw(eps_opt, d_opt, th_opt));
Km_opt  = Km_fun(Jp_opt);
LEP_opt = LEP_fun(d_opt, theta_opt);
%% ====== 2) Journal Article Data ======
eps_journal   = 0.80;
d_journal     = 0.30;
th_journal    = 70.0;
theta_journal = 130;
Jp_journal  = 75.34;                        % updated as requested
Km_journal  = Km_fun(Jp_journal);
LEP_journal = LEP_fun(d_journal, theta_journal);
%% ---- Print results ----
fprintf('\n================ Comparison (Optimized vs Journal Article Data) ================\n');
fprintf('Porosity ε              : %7.4f  vs  %7.4f  (-)\n',  eps_opt, eps_journal);
fprintf('Pore diameter d         : %7.4f  vs  %7.4f  µm\n',  d_opt, d_journal);
fprintf('Thickness δ             : %7.4f  vs  %7.4f  µm\n',  th_opt, th_journal);
fprintf('Contact angle θ         : %7.1f  vs  %7.1f  deg\n',  theta_opt, theta_journal);
fprintf('Permeate flux Jp        : %7.2f  vs  %7.2f  kg m^-2 h^-1\n', Jp_opt, Jp_journal);
fprintf('Mass transfer coeff K_m : %7.2f  vs  %7.2f  kJ m^-2 s^-1 °C^-1\n', Km_opt, Km_journal);
fprintf('LEP                     : %7.2f  vs  %7.2f  kPa\n', LEP_opt, LEP_journal);
%% ---- Plot 1: Design variables ----
figure('Color','w');
cats = categorical({'\epsilon (-)','d (\mum)','\delta (\mum)','\theta (^\circ)'});
bar(cats, [eps_opt d_opt th_opt theta_opt; eps_journal d_journal th_journal theta_journal],'grouped');
ylabel('Value','Interpreter','tex'); grid on;
legend({'Optimized','Journal article data'},'Location','best');
title('Membrane properties: Optimized vs Journal article data','Interpreter','tex');
%% ---- Plot 2: Performance metrics ----
figure('Color','w');
cats2 = categorical({'J_p (kg m^{-2} h^{-1})','K_m (kJ m^{-2} s^{-1} ^\circC)','LEP (kPa)'});
bar(cats2, [Jp_opt Km_opt LEP_opt; Jp_journal Km_journal LEP_journal],'grouped');
ylabel('Value','Interpreter','tex'); grid on;
legend({'Optimized','Journal article data'},'Location','best');
title('Performance: Optimized vs Journal article data (units fixed)','Interpreter','tex');

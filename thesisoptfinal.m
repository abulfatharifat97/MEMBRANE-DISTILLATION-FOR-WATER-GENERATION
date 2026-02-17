%% ===================== MD Optimization via Desirability (Finalized Physics) =====================
clear; clc; close all; rng(1);  % for reproducibility
%% ---------------- Constants (unchanged physics) ----------------
Hfg    = 2.406e6;        % J/kg (latent heat)
Rconst = 461.5;          % J/(kg·K) (water vapor)
Thi    = 66;             % °C (feed hot)
Tci    = 44;             % °C (permeate cold)
Tm     = (Thi + Tci)/2;  % °C (mean)
dT     = Thi - Tci;      % °C (bulk ΔT)
Pm     = 1.013e5;        % Pa (mean pore pressure ~ 1 atm)
gamma  = 0.067;          % N/m (surface tension)
tau    = 2;              % (-) tortuosity (assumed constant)
%% ---------------- Design bounds (from your PDFs) ----------------
eps_lb = 0.48; eps_ub = 0.91;             % porosity (-)
d_lb   = 0.21; d_ub   = 0.98;             % pore DIAMETER (µm)
th_lb  = 70.31; th_ub = 190.18;           % thickness (µm)
th0    = th_lb;                           % baseline thickness used in calibration
th_units = '\mum';
theta_lb = 90; theta_ub = 150;            % contact angle (deg)
%% ---------------- Baseline (for Cw scaling to match your PDF) ----------------
% Baseline from your validated figures:
eps0 = 0.80; d0 = 0.22; theta0 = 120;     % (deg); LEP uses theta only
Jp_baseline_h = 66.79;                    % kg m^-2 h^-1
% Membrane coefficient Cw ∝ (ε * d) / (τ * δ)
Cw_raw  = @(eps, d_um, th_um) (eps .* (d_um*1e-6)) ./ (tau .* (th_um*1e-6));
% Calibrate a scale so baseline reproduces Jp = 66.79 kg m^-2 h^-1
Cw_base    = Cw_raw(eps0, d0, th0);
scaleFact  = (Jp_baseline_h/3600) / ( Cw_base * (Hfg/(Rconst*(Tm+273)*Pm)) * dT );
Cw         = @(eps, d_um, th_um) scaleFact .* Cw_raw(eps, d_um, th_um);
%% ---------------- Finalized performance models (unchanged) ----------------
% Permeate flux Jp (kg m^-2 h^-1)
Jp_h = @(eps, d_um, th_um) 3600 .* Cw(eps, d_um, th_um) .* (Hfg./(Rconst*(Tm+273).*Pm)) .* dT;
% Mass transfer coefficient (kJ m^-2 s^-1 °C^-1)
Km   = @(Jp_h_val) ((Jp_h_val/3600) .* Hfg ./ dT) ./ 1000;
% Liquid Entry Pressure (kPa) using pore DIAMETER (µm)
LEP  = @(d_um, theta_deg) (2*gamma*abs(cosd(theta_deg))) ./ (d_um*1e-6) ./ 1000;
%% ---------------- Desirability setup (Derringer–Suich, larger-is-better) ----------------
% Choose engineering lower bounds (acceptable) and targets (desired)
% If you have specific specs, set them here. Otherwise auto-set from achievable min/max.
% Auto compute achievable extrema at corners (Jp and Km are monotone in eps,d,th; LEP in d,theta)
Jp_min_auto = Jp_h(eps_lb, d_lb, th_ub);
Jp_max_auto = Jp_h(eps_ub, d_ub, th_lb);
Km_min_auto = Km(Jp_min_auto);
Km_max_auto = Km(Jp_max_auto);
LEP_min_auto = LEP(d_ub, theta_lb);       % worst LEP at largest d, smallest |cosθ| (near 90°)
LEP_max_auto = LEP(d_lb, theta_ub);       % best LEP at smallest d, largest |cosθ| (150° in your range)
% ---- User-specifiable desirability thresholds (edit if you have specs) ----
L1 = Jp_min_auto;    T1 = Jp_max_auto;    r1 = 1;   % Jp  (kg m^-2 h^-1)
L2 = Km_min_auto;    T2 = Km_max_auto;    r2 = 1;   % Km  (kJ m^-2 s^-1 °C^-1)
% Safety requirement example: enforce LEP >= 300 kPa by setting L3 = 300
LEP_min_req = 300;                          % kPa  (set [] or NaN to disable)
L3 = max(LEP_min_auto, LEP_min_req);        T3 = LEP_max_auto;  r3 = 1;  % LEP
% Weights (importance) for overall desirability (equal by default)
w1 = 1; w2 = 1; w3 = 1;
% Individual desirability, larger-is-better
d_lob = @(y, L, T, r) (y<=L).*0 + (y>=T).*1 + ...
                      (y>L & y<T).*(((y-L)./(T-L)).^r);
% Overall desirability (weighted geometric mean)
D_overall = @(d1,d2,d3) (d1.^w1 .* d2.^w2 .* d3.^w3).^(1/(w1+w2+w3));
%% ---------------- Optimization (Monte-Carlo search + optional local refine) ----------------
Nsamp = 25000;                                  % number of random samples
eps_s = eps_lb + (eps_ub-eps_lb).*rand(Nsamp,1);
d_s   = d_lb   + (d_ub-d_lb).*rand(Nsamp,1);
th_s  = th_lb  + (th_ub-th_lb).*rand(Nsamp,1);
tht_s = theta_lb + (theta_ub-theta_lb).*rand(Nsamp,1);
% Evaluate performance
Jp_s  = Jp_h(eps_s, d_s, th_s);
Km_s  = Km(Jp_s);
LEP_s = LEP(d_s, tht_s);
% Desirabilities
d1_s = d_lob(Jp_s,  L1, T1, r1);
d2_s = d_lob(Km_s,  L2, T2, r2);
d3_s = d_lob(LEP_s, L3, T3, r3);
D_s   = D_overall(d1_s, d2_s, d3_s);
% Pick best sample
[~, idx_best] = max(D_s);
eps_opt   = eps_s(idx_best);
d_opt     = d_s(idx_best);
th_opt    = th_s(idx_best);
theta_opt = tht_s(idx_best);
% Compute optimal performance
Jp_opt  = Jp_h(eps_opt, d_opt, th_opt);
Km_opt  = Km(Jp_opt);
LEP_opt = LEP(d_opt, theta_opt);
D_opt   = D_s(idx_best);
%% ---------------- Results printout ----------------
fprintf('\n================= Desirability Optimization Result =================\n');
fprintf('Porosity (ε)                 = %.4f (-)\n',   eps_opt);
fprintf('Pore diameter (d)            = %.4f %s\n',     d_opt, th_units);
fprintf('Thickness (δ)                = %.4f %s\n',     th_opt, th_units);
fprintf('Contact angle (θ)            = %.1f °\n',      theta_opt);
fprintf('Permeate flux                = %.2f kg m^-2 h^-1\n', Jp_opt);
fprintf('Mass transfer coefficient    = %.2f kJ m^-2 s^-1 °C^-1\n', Km_opt);
fprintf('LEP                          = %.2f kPa\n', LEP_opt);
fprintf('Overall desirability (D)     = %.4f (0–1)\n', D_opt);
%% ---------------- Visualizations ----------------
% 1) Trade-off across pore diameter at optimal ε, δ, θ
d_grid = linspace(d_lb, d_ub, 400);
Jp_g   = Jp_h(eps_opt, d_grid, th_opt);
Km_g   = Km(Jp_g);
LEP_g  = LEP(d_grid, theta_opt);
d1_g   = d_lob(Jp_g,  L1, T1, r1);
d2_g   = d_lob(Km_g,  L2, T2, r2);
d3_g   = d_lob(LEP_g, L3, T3, r3);
D_g    = D_overall(d1_g, d2_g, d3_g);
figure('Color','w');
yyaxis left;  plot(d_grid, Jp_g, 'b-', 'LineWidth',1.6); ylabel('Permeate flux (kg m^{-2} h^{-1})', 'Interpreter','tex');
yyaxis right; plot(d_grid, LEP_g,'r--','LineWidth',1.6); ylabel('LEP (kPa)', 'Interpreter','tex');
xlabel('Pore diameter (\mum)', 'Interpreter','tex'); grid on; hold on;
yyaxis left;  plot(d_opt, Jp_opt, 'bo', 'MarkerFaceColor','b');
yyaxis right; plot(d_opt, LEP_opt,'ro', 'MarkerFaceColor','r');
title('Trade-off across pore diameter at optimal \epsilon,\delta,\theta', 'Interpreter','tex');
legend({'Permeate flux','LEP','Optimum'},'Location','best');
% 2) Overall desirability vs pore diameter (others fixed at optimum)
figure('Color','w');
plot(d_grid, D_g, 'k-', 'LineWidth',1.6); grid on; hold on;
plot(d_opt, D_overall(d_lob(Jp_opt,L1,T1,r1), d_lob(Km_opt,L2,T2,r2), d_lob(LEP_opt,L3,T3,r3)), ...
     'ko', 'MarkerFaceColor','k');
xlabel('Pore diameter (\mum)', 'Interpreter','tex');
ylabel('Overall desirability, D (0–1)', 'Interpreter','tex');
title('Overall desirability across pore diameter', 'Interpreter','tex');
% 3) Bar chart of optimal design variables
figure('Color','w');
cats = categorical({'\epsilon (-)','d (\mum)','\delta (\mum)','\theta (^\circ)'});
vals = [eps_opt, d_opt, th_opt, theta_opt];
bar(cats, vals); grid on;
ylabel('Value', 'Interpreter','tex');
title('Optimized membrane properties', 'Interpreter','tex');
% 4) Desirability components at optimum
figure('Color','w');
d1_opt = d_lob(Jp_opt,  L1, T1, r1);
d2_opt = d_lob(Km_opt,  L2, T2, r2);
d3_opt = d_lob(LEP_opt, L3, T3, r3);
bar(categorical({'d_1 (Flux)','d_2 (K_m)','d_3 (LEP)','D overall'}), [d1_opt, d2_opt, d3_opt, D_opt]);
ylim([0,1]); grid on;
ylabel('Desirability (0–1)', 'Interpreter','tex');
title('Desirability at the optimum', 'Interpreter','tex');
%% ---------------- (Optional) Local refinement using fmincon ----------------
% Convert to maximization by minimizing -D; keep bounds.
% Requires Optimization Toolbox. Comment out if not available.
try
    x0 = [eps_opt, d_opt, th_opt, theta_opt];
    lb = [eps_lb, d_lb, th_lb, theta_lb];
    ub = [eps_ub, d_ub, th_ub, theta_ub];
    D_neg = @(x) -D_overall( ...
        d_lob(Jp_h(x(1), x(2), x(3)), L1, T1, r1), ...
        d_lob(Km(Jp_h(x(1), x(2), x(3))), L2, T2, r2), ...
        d_lob(LEP(x(2), x(4)), L3, T3, r3) );
    opts = optimoptions('fmincon','Display','none','Algorithm','sqp');
    [x_ref, fval] = fmincon(D_neg, x0, [],[],[],[], lb, ub, [], opts);
    % Update if improved
    if -fval > D_opt
        eps_opt = x_ref(1); d_opt = x_ref(2); th_opt = x_ref(3); theta_opt = x_ref(4);
        Jp_opt  = Jp_h(eps_opt, d_opt, th_opt);
        Km_opt  = Km(Jp_opt);
        LEP_opt = LEP(d_opt, theta_opt);
        D_opt   = -fval;
        fprintf('\n(Local refine) Improved overall desirability D = %.4f\n', D_opt);
    end
catch
    % If no Optimization Toolbox, skip refinement
end
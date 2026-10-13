%% Generate reference data for validating the Julia port against the MATLAB model.
%
%  Run from MATLAB (any working directory):
%     run('D:/bile-acid-julia/test/reference/generate_reference.m')
%
%  Writes CSV files next to this script:
%    rhs_x.csv, rhs_t.csv, rhs_p.csv, rhs_c.csv, rhs_dx.csv   - right-hand-side checks
%    sim_tight_t.csv, sim_tight_x.csv                           - tight-tolerance simulation
%    sim_matlab_t.csv, sim_matlab_x.csv, sim_matlab_time.csv     - simulate_human_day as shipped

here = fileparts(mfilename('fullpath'));
pkg  = fullfile(here, '..', '..', 'BA_model_package');
addpath(genpath(pkg));
load(fullfile(pkg, 'Results', 'refined_result_MST.mat'));
mi    = result.model_info;
p_opt = result.p_opt(:);
c_opt = mi.c_vec(:);
f     = @BA_human_model_v9_to_fit_full_ode;
w     = @(name, M) dlmwrite(fullfile(here, name), M, 'delimiter', ',', 'precision', '%.17g');

%% 1. Right-hand side at random states -------------------------------------
rng(1234);
N  = 40;
X  = zeros(N, 402); T = zeros(N, 2); P = zeros(N, 33); C = zeros(N, 128); DX = zeros(N, 402);
for i = 1:N
    p = p_opt; c = c_opt;
    if i > N/2
        % Perturb parameters and switch on every constant that is zero in the
        % calibrated set, so that each flux term in the model is exercised.
        p      = p_opt .* exp(0.3 * randn(33, 1));
        c(2:7) = 1e-4 * (1 + rand(6, 1));          % passive uptake of conjugated BA
        c(8)   = c(8) * (0.8 + 0.4 * rand);        % V_P
        c(9)   = c(9) * (0.8 + 0.4 * rand);        % Q_L
        c(14)  = 0.05 * rand;                      % k_bact_psi
        c(15)  = 0.5 + rand;                       % k_bact_co
        c([20 24 25 26 30]) = 0.02 * rand(5, 1);   % synthesis fractions
        c(32)  = 0.01 * rand;                      % k_pl
    end
    x        = 10.^(-2 + 5 * rand(402, 1));
    tau      = 1000 * rand;
    t        = tau + (-60 + 700 * rand);           % includes pre-meal times (reflex clamped at 1)
    c(1)     = tau;
    X(i, :)  = x'; T(i, :) = [t tau]; P(i, :) = p'; C(i, :) = c';
    DX(i, :) = f(t, x, p, c)';
end
w('rhs_x.csv', X); w('rhs_t.csv', T); w('rhs_p.csv', P); w('rhs_c.csv', C); w('rhs_dx.csv', DX);

%% 2. Tight-tolerance simulation of the standard protocol ---------------------
%  20 days without meals (tau = 0), then 35 days of meals with gaps [6 6 12] h.
opts    = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
x0      = mi.x0(:);
t_ss    = mi.doss * 24 * 60;
c       = c_opt; c(1) = 0;
[~, xs] = ode15s(f, [0 t_ss/2 t_ss], x0, opts, p_opt, c);
Tout    = [0; t_ss]; Xout = [x0'; xs(end, :)];
gaps    = repmat(mi.gap(:)', 1, mi.dors) * 60;
t0      = t_ss; x = xs(end, :)';
for k = 1:numel(gaps)
    t1 = t0 + gaps(k);
    if k > numel(gaps) - numel(mi.gap)              % last day: output every 30 min
        tspan = t0:30:t1;
    else
        tspan = [t0 (t0 + t1)/2 t1];
    end
    c(1)    = t0;
    [tt, xx] = ode15s(f, tspan, x, opts, p_opt, c);
    if k > numel(gaps) - numel(mi.gap)
        Tout = [Tout; tt(2:end)]; Xout = [Xout; xx(2:end, :)];
    else
        Tout = [Tout; tt(end)];   Xout = [Xout; xx(end, :)];
    end
    t0 = t1; x = xx(end, :)';
end
w('sim_tight_t.csv', Tout); w('sim_tight_x.csv', Xout);

%% 3. simulate_human_day exactly as shipped (ode15s, default tolerances) -----
tic;
[static, days, day] = simulate_human_day(p_opt, mi.doss, mi.dors, mi.x0, mi, 0, 1);
elapsed = toc;
meal_ends = t_ss + cumsum(gaps);
keep      = ismember(days.t, [t_ss meal_ends]);
w('sim_matlab_t.csv', days.t(keep)); w('sim_matlab_x.csv', days.x(keep, :));
w('sim_matlab_time.csv', elapsed);
fprintf('simulate_human_day: %.2f s, %d output points\n', elapsed, numel(days.t));

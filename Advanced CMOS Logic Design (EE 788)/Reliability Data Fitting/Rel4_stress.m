clear all; 
close all;

filename = "data_assignment_Rel-04.xlsx";
sheet = 1;

data = readmatrix(filename,'Sheet',sheet);
% data(any(isnan(data),2),:) = [];

% From data set

T_vals = [25, 65, 100, 125] + 273;
Vg_vals = [1.3, 1.5, 1.7, 1.9];

stress_times = rmmissing(data(:,1));
del_VT = rmmissing(data(:, 2:end)) * 1e-3; % Data is given in mV 
 
%% Parameters

t_il = 0.3e-9; 
t_hk = 2.3e-9;

eps_0 = 8.854e-12;

eps_r_ox = 3.9;
eps_ox = eps_0 * eps_r_ox;

eps_il = 3.9;
eps_hk = 25;

EOT = eps_r_ox*(t_il/eps_il + t_hk/eps_hk);

q = 1.6e-19; 
Kb = 1.380649e-23;
C_ox = eps_ox/EOT;   

%% n value

n_vals = zeros(1,size(del_VT,2));
for i = 1:size(n_vals,2)
    fit1 = polyfit(log(stress_times), log(del_VT(:,i)),  1);
    n_vals(i) = fit1(1);
end
n_avg = mean(n_vals);

%% Ea (at fixed Vg)

Ea_vals = zeros(size(Vg_vals,2),length(stress_times));

for Vg_index = 1:size(Vg_vals,2)
    for j = 1:length(stress_times)

        del_VT_arr = del_VT(j, [Vg_index, Vg_index + size(Vg_vals,2), Vg_index + 2*size(Vg_vals,2), Vg_index + 3*size(Vg_vals,2)]);
        if (~any(isnan(del_VT_arr)))

            fit_Ea = polyfit(-1*q./(Kb*T_vals), log(del_VT_arr),  1);
            Ea_vals(Vg_index, j) = fit_Ea(1);

        end
    end
end

Ea_vals(Ea_vals ==  0) = nan;
Ea_avg = mean(Ea_vals, 'all', 'omitnan');


%% VAF (at fixed T)

gamma_vals = zeros(size(T_vals,2),length(stress_times));

for T_idx = 1:size(T_vals,2)
    for i = 1:length(stress_times)
        del_VT_arr = del_VT(i, (T_idx - 1)*size(Vg_vals,2) + (1:size(Vg_vals,2)));
        if (~any(isnan(del_VT_arr)))
            fit_gamma = polyfit(log(Vg_vals), log(del_VT_arr),  1);
            gamma_vals(T_idx, i) = fit_gamma(1);
        end
    end
end

gamma_vals(gamma_vals ==  0) = nan;
gamma_avg = mean(gamma_vals, 'all', 'omitnan');


%% A 

del_VT_func = @(A, Vg, gamma, Ea, T, n, stress_times) A .* kron(exp(-q * Ea./(Kb*T)) , Vg.^gamma) .*  stress_times.^n;

optim_func = @(params) norm((del_VT_func(params(1), Vg_vals, gamma_avg, Ea_avg, T_vals, n_avg, stress_times) - del_VT));

param_guess = 2e-1;
% init_guess = optim_func(param_guess)

%
% options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);
options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals', 1e30, 'Maxiter', 1e30);
A_avg = fminsearch(optim_func, param_guess, options);


disp("Fit Parameters found: ");
A_avg
gamma_avg 
Ea_avg
n_avg

%% Plotting

data_fitted = del_VT_func(A_avg, Vg_vals, gamma_avg, Ea_avg, T_vals, n_avg, stress_times);

figure;

Vg_vals = [Vg_vals Vg_vals];
for i = 1:16
    Vg = Vg_vals(4 + mod(i, 4));
    if i <= 4
        T  = T_vals(1);
    elseif i > 4 && i <= 8
        T  = T_vals(2);
    elseif i > 8 && i <= 12
        T  = T_vals(3);
    else
        T  = T_vals(4);
    end

    color = [0.1 * randi([0,10]), 0.2 * randi([0,5]), 0.1 * randi([0,10])];
    not_nan_idx = ~isnan(del_VT(:,i));
    scatter(stress_times(not_nan_idx), del_VT(not_nan_idx,i),[], color, 'filled','DisplayName', '');
    hold on;
    grid on;
    set(gca,'xscale','log'); 
    set(gca,'yscale','log');
    plot(stress_times, data_fitted(:,i), '-','Color', color, 'DisplayName',...
        "Fitted Data - V_g = " + Vg + ", Temp (in C): " + (T - 273));
    legend('-DynamicLegend', 'location', 'best','FontSize',7);
end   

xlabel('Stress time(s)', 'interpreter', 'latex');
ylabel("$\Delta V_{IT}(V)$",  'interpreter', 'latex');
title("$\Delta V_{IT}$ vs Stress time (Stress data) : Slope = " + n_avg,  'interpreter', 'latex');
xlim padded;



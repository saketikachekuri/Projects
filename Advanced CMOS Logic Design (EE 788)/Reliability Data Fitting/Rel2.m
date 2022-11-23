clear all; close all;

filename = "data_assignment_Rel-02.xlsx";
sheet = 1;

data = readmatrix(filename,'Sheet',sheet);

%% Split data into arrays

stress_times = data(:,1);
delay1 = 10;
delay2 = 100;

del_NIT1 = data(:,3)*1e4; %Since given in /cm^2
del_NIT2 = data(:,2)*1e4;

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


del_VT_1 = del_NIT1*q/C_ox;
del_VT_2 = del_NIT2*q/C_ox;
%% Fit given data

fit_data1 = polyfit(log10(stress_times), log10(del_VT_1), 1);
fit_data2 = polyfit(log10(stress_times), log10(del_VT_2), 1);

slope_data1 = fit_data1(1);
int_data1 = fit_data1(2);
% Estimates according to polyfit
est_del_VT_1 = 10.^(slope_data1*log10(stress_times) + int_data1);

slope_data2 = fit_data2(1);
int_data2 = fit_data2(2);
% Estimates according to polyfit
est_del_VT_2 = 10.^(slope_data2*log10(stress_times) + int_data2);


% For Rel4 - find sep A, gamma, k, m for 16, and also 'universal' using ML


% Extra:
% Use (VG- Vt)^gamma, with Vt from prev time point (at t=0, Vt = 0), you will be able to model the
% field reduction. When we extrapolate, field reduction also into account
%% Optimization

t_arr1 = delay1./stress_times;
t_arr2 = delay2./stress_times;


VT_EOS_func = @(int) (10.^((1/6)*log10(stress_times) + int));
rec_factor_func = @(k, m, t_arr) (1 + k*(t_arr).^m);
opt_func = @(int, k, m, t_arr1, del_VT_1, t_arr2, del_VT_2)(abs(VT_EOS_func(int)./rec_factor_func(k, m, t_arr1) - del_VT_1) + abs(VT_EOS_func(int)./rec_factor_func(k, m, t_arr2) - del_VT_2));

norm_func = @(var_arr) norm(opt_func(var_arr(1), var_arr(2), var_arr(3), t_arr1, del_VT_1, t_arr2, del_VT_2));

avg_int = mean([int_data1, int_data2]);

%% fminsearch 

% var_arr_init = [avg_int, 0.5, 0.5];
% norm_func(var_arr_init)
% options = optimset('PlotFcns', @optimplotfval, 'TolFun',1e-6,'TolX',1e-6);
% var_arr_opt = fminsearch(norm_func, var_arr_init, options)

%% use fmincon instead of fminsearch, so that you can put bounds on parameters

var_arr_init = [avg_int, 0.5, 0.5]';
options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);
[var_arr_opt, minval] = fmincon(norm_func, var_arr_init, [],[],[],[],[-1e3, 0, 0]',[1e3, 1, 1]',[],options);

%% Plotting
slope = 1/6
est_int = var_arr_opt(1)
est_A = 10^est_int
est_k = var_arr_opt(2)
est_m = var_arr_opt(3)

est_del_VT_EOS = 10.^(slope*log10(stress_times) + est_int)


figure(1);
grid on;
loglog(stress_times, del_VT_1, 'bx');
hold on;
loglog(stress_times, del_VT_2, 'rx');
hold on;
loglog(stress_times, est_del_VT_1, 'b-');
hold on;
loglog(stress_times, est_del_VT_2, 'r-');
hold on;
loglog(stress_times, est_del_VT_EOS, 'k-');
hold on;

xlabel('Stress time(s)', 'interpreter', 'latex');
ylabel('$\Delta V_{T}(V)$', 'interpreter', 'latex');
title("$\Delta V_{T}$ vs Stress time",'interpreter', 'latex');
str_1 = sprintf('Data fit for 10s, Slope = %.4f', slope_data1);
str_2 = sprintf('Data fit for 100s, Slope = %.4f', slope_data2);

str_corr = sprintf('Corrected data fit, Slope = %.4f', slope);

legend("Data for delay = 10s", "Data for delay = 100s", str_1, str_2, str_corr, 'interpreter', 'latex', 'Location', 'southeast');

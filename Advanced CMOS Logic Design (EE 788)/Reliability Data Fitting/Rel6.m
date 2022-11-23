clear all; 
close all;

filename = "data_assignment_Rel-05.xlsx";
sheet = 1;

data = readmatrix(filename,'Sheet',sheet);
% data(any(isnan(data),2),:) = [];

% From data set

T_vals = [100, 130, 155] + 273;
Vg_vals = [1.3, 1.5, 1.7];

stress_times = rmmissing(data(:,1));
del_VT = rmmissing(data(:, 2:end)) * 1e-3; % Data is given in mV 
 
%% Parameters
q = 1.6e-19; 
Kb = 1.380649e-23;


%% A 

del_VT_func_IT = @(A_IT, gamma_IT, Ea_IT)...
    A_IT .* kron(exp(-q * Ea_IT./(Kb*T_vals)) , Vg_vals.^gamma_IT) .*  stress_times.^(1/6);

del_VT_func_HT = @(A_HT, gamma_HT, Ea_HT)...
    A_HT .* kron(exp(-q * Ea_HT./(Kb*T_vals)) , Vg_vals.^gamma_HT) .*  stress_times.^(0);

del_VT_func_OT = @(A_OT, gamma_OT, Ea_OT, n_OT)...
    A_OT .* kron(exp(-q * Ea_OT./(Kb*T_vals)) , Vg_vals.^gamma_OT) .*  stress_times.^(n_OT);


optim_func = @(params) norm(del_VT_func_IT(params(1),params(2), params(3)) + del_VT_func_HT(params(4),params(5), params(6)) + del_VT_func_OT(params(7), params(8), params(9), params(10)) - del_VT);

param_guess = [30 , 4.5, 0.08, 3 , 4.5, 0.05, 20, 9, 0.35, 0.3]';

param_lower = [1e-4, 3.5, 0.06, 1e-4, 3.5, 0.06, 15, 4, 0.2, 1/4]';
param_upper = [1e4, 5.5, 0.14, 1e4, 5.5, 0.14, 25, 13, 0.5, 1/3]';

options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);

[param_vals, minval] = fmincon(optim_func, param_guess, [], [], [], [], param_lower, param_upper, [], options);


fprintf ("Cost function: %f\n\n\nFit Parameters:\n\nA_IT = %f \ngamma_IT = %f \nEa_IT = %f \n\n" + ...
    "A_HT = %f \ngamma_HT = %f \nEa_HT = %f \n \n\nA_OT = %f \ngamma_OT = %f \nEa_OT = %f \n",minval, param_vals(1), param_vals(2), param_vals(3), param_vals(4), param_vals(5), param_vals(6), param_vals(7), param_vals(8), param_vals(9));

fprintf("\n\n\nSlopes\n\nn_IT = %f\nn_HT = %f\nn_OT = %f\n", 1/6, 0, param_vals(10))
%% Plotting

data_fitted_IT = del_VT_func_IT(param_vals(1), param_vals(2), param_vals(3));
data_fitted_HT = del_VT_func_HT(param_vals(4), param_vals(5), param_vals(6));
data_fitted_OT = del_VT_func_OT(param_vals(7), param_vals(8), param_vals(9), param_vals(10));

data_fitted = data_fitted_IT + data_fitted_HT + data_fitted_OT;

contrib_IT = zeros(1,9);
contrib_HT = zeros(1,9);
contrib_OT = zeros(1,9);

Vg_vals = [Vg_vals Vg_vals];

figure;

for i = 1:9

    Vg = Vg_vals(3 + mod(i, 3));
    if i <= 3
        T  = T_vals(1);
    elseif i > 3 && i <= 6
        T  = T_vals(2);
    else
        T  = T_vals(3);
    end

    contrib_IT(i) = mean(abs(data_fitted_IT(:,i)./data_fitted(:,i)));
    contrib_HT(i) = mean(abs(data_fitted_HT(:,i)./data_fitted(:,i)));
    contrib_OT(i) = mean(abs(data_fitted_OT(:,i)./data_fitted(:,i)));

    color = [0.1 * randi([0,10]), 0.2 * randi([0,5]), 0.1 * randi([0,10])];

    scatter(stress_times, del_VT(:,i),[], color, 'filled','DisplayName', '');
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
title("$\Delta V_{IT}$ vs Stress time with IT, HT, OT contributions separate",  'interpreter', 'latex');
xlim padded;


%% Bar graph
figure;

subplot(1,3,1);
y1 = reshape(contrib_IT, 3, 3)';
opt = reordercats(categorical({'T = 100 C','T = 130C','T = 150C'}),{'T = 100 C','T = 130C','T = 150C'});
b1 = bar(opt, y1*100);

xlabel('Temperature (C)', 'interpreter', 'latex')
ylabel('Percentage Contribution of $\Delta V_{IT}$', 'interpreter', 'latex')
title("Contribution of $\Delta V_{IT}$ component", 'interpreter', 'latex');
legend('$V_G = 1.3 V$', '$V_G = 1.5 V$', '$V_G = 1.7 V$', 'interpreter', 'latex', 'Location', 'best');
xlim padded;



subplot(1,3,2);
y2 = reshape(contrib_HT, 3, 3)';
opt = reordercats(categorical({'T = 100 C','T = 130C','T = 150C'}),{'T = 100 C','T = 130C','T = 150C'});
b2 = bar(opt, y2*100);

xlabel('Temperature (C)', 'interpreter', 'latex')
ylabel('Percentage Contribution of $\Delta V_{HT}$', 'interpreter', 'latex')
title("Contribution of $\Delta V_{HT}$ component", 'interpreter', 'latex');
legend('$V_G = 1.3 V$', '$V_G = 1.5 V$', '$V_G = 1.7 V$', 'interpreter', 'latex', 'Location', 'best');
xlim padded;

subplot(1,3,3);
y3 = reshape(contrib_OT, 3, 3)';
opt = reordercats(categorical({'T = 100 C','T = 130C','T = 150C'}),{'T = 100 C','T = 130C','T = 150C'});
b3 = bar(opt, y3*100);

xlabel('Temperature (C)', 'interpreter', 'latex')
ylabel('Percentage Contribution of $\Delta V_{OT}$', 'interpreter', 'latex')
title("Contribution of $\Delta V_{OT}$ component", 'interpreter', 'latex');
legend('$V_G = 1.3 V$', '$V_G = 1.5 V$', '$V_G = 1.7 V$', 'interpreter', 'latex', 'Location', 'best');
xlim padded;

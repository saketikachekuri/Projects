clear all; 
close all;

stress = 1000;
filename = "data_assignment_Rel-04.xlsx";
% From data set

T_vals = [25, 65, 100, 125] + 273;
Vg_vals = [1.3, 1.5, 1.7, 1.9];

sheet = 2;
data = readmatrix(filename,'Sheet',sheet);
delay_times = rmmissing(data(:,1));
del_VT = rmmissing(data(:, 2:end)) * 1e-3; % Data is given in mV 

sheet = 1;
data1 = readmatrix(filename,'Sheet',sheet);
del_VT_EOS = rmmissing(data1(end, 2:end))* 1e-3;    % Data is given in mV 
del_VT_EOS_mat = repelem(del_VT_EOS, size(delay_times,1), 1);

%% Find universal k and m

t_arr = delay_times/stress; %7x1
t_mat = kron(t_arr,ones(1,16));
%%
rec_func = @(k, m)...
    (sum(abs(del_VT_EOS_mat./(1 + k*(t_mat).^m) - del_VT), 'all')); 

norm_func = @(var_arr) (rec_func(var_arr(1), var_arr(2)));

var_arr_init = [0.5, 0.5]';

options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);
[var_arr_opt, minval] = fmincon(norm_func, var_arr_init, [],[],[],[],[0, 0]',[2, 1]',[],options);

var_arr_opt

del_VT_final = del_VT_EOS./(1+ var_arr_opt(1)*(t_mat).^var_arr_opt(2));
%% Plot data
% % plot given recovery data and fitted recovery data (from formula) to
% observe accuracy of fit
figure;
for i = 1:16
    Vg = Vg_vals(1+mod(i-1, 4));
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
    % plot actual recovery data
    scatter(delay_times, del_VT(:,i), '<','MarkerEdgeColor', color, 'DisplayName',...
        "Given Data - V_g(V) = " + Vg + ", Temp(C) = " + (T - 273));
    hold on;
    set(gca,'xscale','log'); 
    set(gca,'yscale','log');
    grid on;
    % plot theoretical recovery estimates
    scatter(delay_times, del_VT_final(:,i), '>','MarkerEdgeColor', color, 'DisplayName',...
        "Fitted Data - V_g(V) = " + Vg + ", Temp(C) = " + (T - 273));
    hold  on;
    legend('-DynamicLegend', 'location', 'best','FontSize',7)
end   
xlabel('Recovery time (s)', 'interpreter', 'latex')
ylabel('$\Delta V_{T}$ $(V)$', 'interpreter', 'latex')
title("$\Delta V_{T}$ $(V)$ vs Recovery time", 'interpreter', 'latex');
xlim padded;

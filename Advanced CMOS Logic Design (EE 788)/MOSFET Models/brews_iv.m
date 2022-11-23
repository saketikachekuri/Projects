% All calculations are assuming type = NMOS
%% Change parameters here
idx = 3;

u_n_fixed = 200e-4; % Ball-park figure for constant mu calculations

% Taking width as 1 um so current is obtained per unit um (other values
% would just scale it by W, if width is W um

W = 1e-6;
kB = 1.380649e-23;
T = 300;

%% Values
L_arr = [1, 0.5, 0.35];
tox_arr = [20, 10, 7];
% From graph
VDD_arr = [5, 3.5, 3];
Vth_arr = [0.8, 0.55, 0.5];

% Converting to given units
L = L_arr(idx) * 10^(-6);
tox = tox_arr(idx) * 10^(-9);
VDD = VDD_arr(idx);
Vth = Vth_arr(idx);

%% Constants for reuse
q = 1.6e-19;
VT = kB*T/q;
eps_0 = 8.854e-12;

eps_r_ox = 3.7;
eps_ox = eps_r_ox * eps_0;

eps_r_si = 11.7;
eps_si = eps_r_si * eps_0;

Cox = eps_ox/tox;

n_i = (1.5e10) * 1e6;  % SI units
Eg = 1.12*q;
Xs = 4.05;
phi_m = Xs; % Conduction band since polysilicon

%% Find NA for given Vt

NA_arr = 10.^(linspace(21, 26, 1000)); %In SI units
phi_b_arr = VT * log(NA_arr/n_i);
phi_ms_arr = phi_m - (Xs + Eg/(2*q) + phi_b_arr);

% Assuming no oxide/interface charges, V_fb is phi_ms 

Vth_vals = phi_ms_arr + sqrt(4*q*eps_si*(NA_arr.*phi_b_arr))/(Cox) + 2*phi_b_arr;

NA = interp1(Vth_vals, NA_arr, Vth);

phi_b = VT * log(NA/n_i);
phi_ms = phi_m - (Xs + Eg/(2*q) + phi_b);

% Sanity check; Vth_calc matched Vth from graph, can use Vth in further
% calcs
Vth_calc = phi_ms + sqrt(4*q*eps_si*(NA*phi_b))/(Cox) + 2*phi_b;

V_FB = phi_ms;



%% Code for Id-Vd plots at different Vg vals

Vgs_arr = [2.5 3.5 4.5];
Vds_arr = 0 : 0.01 : VDD;

figure;

for Vgs = Vgs_arr
    Id_vals = [];
    for Vds = Vds_arr
        psi_s_arr = 0: 0.01 : 2*VDD;
        Vgs_values_psi_ss = V_FB + psi_s_arr + (sqrt(2*eps_si*kB*T*NA)/Cox)*sqrt(q*psi_s_arr/(kB*T) + ((n_i/NA)^2)*exp(q*(psi_s_arr - 0)/(kB*T)));
        Vgs_values_psi_sd = V_FB + psi_s_arr + (sqrt(2*eps_si*kB*T*NA)/Cox)*sqrt(q*psi_s_arr/(kB*T) + ((n_i/NA)^2)*exp(q*(psi_s_arr - Vds)/(kB*T)));
        
        psi_ss = interp1(Vgs_values_psi_ss, psi_s_arr, Vgs);
        psi_sd = interp1(Vgs_values_psi_sd, psi_s_arr, Vgs);

        dpsi = (psi_sd - psi_ss)/100;
        psi_arr = psi_ss : dpsi : psi_sd;

        ints = [];

        for psi = psi_arr
            term = (Cox*(Vgs - V_FB - psi) - sqrt(2*eps_si*q*NA*psi) + (2*kB*T/q)*(Cox^2*(Vgs - V_FB - psi) + eps_si*q*NA)/(Cox*(Vgs - V_FB - psi) + sqrt(2*eps_si*q*NA*psi)));
            term = term * dpsi * u_n_fixed * (W/L);

            ints = [ints term];
        end
        Id = sum(ints);
        Id_vals = [Id_vals Id];
    end
    plot(Vds_arr, Id_vals);
    hold on
    
end

xlabel('$V_{DS}$ (in V)', 'interpreter', 'latex');
ylabel('$I_{DS}$ (in A)', 'interpreter', 'latex');
title(['Brews Model $I_{D}$-$V_{D}$ characteristics for $t_{ox}$ = ', num2str(tox*1e9), ' nm, L = ', num2str(L*1e6),' $\mu$m'], 'interpreter', 'latex');
legend("$V_{G}$ = " + Vgs_arr(1) + " V", "$V_{G}$ = " + Vgs_arr(2) + " V", "$V_{G}$ = " + Vgs_arr(3) + " V", 'interpreter', 'latex',Location='best');
grid on


%% Code for Id-Vg plots at different Vd vals

Vds_arr = [0.5 2 3.5];
Vgs_arr = 0 : 0.01 : VDD;

figure;

for Vds = Vds_arr
    Id_vals = [];
    for Vgs = Vgs_arr
        psi_s_arr = 0: 0.01 : 2*VDD;
        Vgs_values_psi_ss = V_FB + psi_s_arr + (sqrt(2*eps_si*kB*T*NA)/Cox)*sqrt(q*psi_s_arr/(kB*T) + ((n_i/NA)^2)*exp(q*(psi_s_arr - 0)/(kB*T)));
        Vgs_values_psi_sd = V_FB + psi_s_arr + (sqrt(2*eps_si*kB*T*NA)/Cox)*sqrt(q*psi_s_arr/(kB*T) + ((n_i/NA)^2)*exp(q*(psi_s_arr - Vds)/(kB*T)));
        
        psi_ss = interp1(Vgs_values_psi_ss, psi_s_arr, Vgs);
        psi_sd = interp1(Vgs_values_psi_sd, psi_s_arr, Vgs);

        dpsi = (psi_sd - psi_ss)/100;
        psi_arr = psi_ss : dpsi : psi_sd;

        ints = [];

        for psi = psi_arr
            term = (Cox*(Vgs - V_FB - psi) - sqrt(2*eps_si*q*NA*psi) + (2*kB*T/q)*(Cox^2*(Vgs - V_FB - psi) + eps_si*q*NA)/(Cox*(Vgs - V_FB - psi) + sqrt(2*eps_si*q*NA*psi)));
            term = term * dpsi * u_n_fixed * (W/L);

            ints = [ints term];
        end
        Id = sum(ints);
        Id_vals = [Id_vals Id];
    end
    plot(Vgs_arr, Id_vals);
    hold on
    
end

xlabel('$V_{GS}$ (in V)', 'interpreter', 'latex');
ylabel('$I_{DS}$ (in A)', 'interpreter', 'latex');
title(['Brews Model $I_{D}$-$V_{G}$ characteristics for $t_{ox}$ = ', num2str(tox*1e9), ' nm, L = ', num2str(L*1e6),' $\mu$m'], 'interpreter', 'latex');
legend("$V_{D}$ = " + Vds_arr(1) + " V", "$V_{D}$ = " + Vds_arr(2) + " V", "$V_{D}$ = " + Vds_arr(3) + " V", 'interpreter', 'latex',Location='best');
grid on
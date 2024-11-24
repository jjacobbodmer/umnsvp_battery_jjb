function SVP_BATT_TUBE_PIPE_COMPARISON()
% SVP BATT TUBE VS PIPE COOLING JOSHUA BODMER 2022-12-28

% BERGMAN HEAT TRANSFER and GERHART FLUID MECHANICS are recommended to have
% on hand. They are available on libgen.rs and are used in the UMN MechE
% courses ME3332, ME3333.

% The goal of this is to determine whether a flat-plate cooling method will
% be more useful than the tube-bank method used to cool our li-ion cells in
% previous batteries.

% The primary goal of UMNSVP should be to win solar races. Skills desired
% by employers are either incidentally or purposefully honed by members of
% this program during their participation therein; it is necessary for that
% skill-building to occur if the vehicles we make will ever win.

% In the case of ASC and WSC, our goal is to maximize _S_, subject to
% race-provided and team-defined regulations placed upon all systems.
% Specifically with regards to WSC, maximizing _S_ is almost perfectly
% equivalent to minimizng _E_, subject to the aforementioned regs plus some
% maximium value for _l_, our lateness penalty to each checkpoint (I say
% this because there is ...probably... a vehicle design out there that
% makes it late to each checkpoint but uses disproportionately less energy
% to do so, offsetting the penalty). It's reasonable to just set _l_ to
% zero.

% With regards to battery design, there are 3 avenues to minimize _E_: (1)
% minimize battery mass -> rolling and gravity resistance -> _E_ (2)
% minimize aerodynamic resistance -> _E_ (3) minimize heat losses -> _E_

% Heat losses due to electrical resistance are likely an order of magnitude
% less than the previous two sources and are "solved" through cell choice
% more than any other thing.

% Mass and drag, however, are both significant. As of 2022-12-28, there is
% little that modifying battery mass can do at this point to reduce total
% vehicle mass but reducing mass in one system can be a virtuous cycle if
% performed by all others in turn. It can (and did for G1) turn vicious
% instead.

% Currently, aero losses are still able to be optimized. Think of the ducts
% on the vehicle as locations where we can modify the static pressure
% against that surface, possibly reducing it. There are certain locations
% on the front of the pontoons of the vehicle where this is most apparent.
% Minimizing duct-battery-duct aerodynamic resistance is then a viable
% avenue in increasing vehicle competitiveness.

% The optimization problem we have to solve on battery is thus:
%       minimize delta_p_batt_duct s.t. (WSC) max cell temperature < 60 C
%       (WSC) max cell voltage < 4.2 (?) (WSC) min cell voltage > 3.3 (?)
%       (WSC) SF > 1 for box in 20g loads (WSC) must be able to identify
%       cells (UMNSVP) must be able to remove/insert battery (UMNSVP) total
%       battery mass < 200 kg

% One thing about the above formulation is that minimizing pressure drop
% also directly impacts the quantity and type of fans (pumps) we need to
% achieve this flow during vehicle charging.

clear
close all
tic

%% LOOP THROUGH DIFFERENT BATTERY CURRENTS

I_target = 80;      % A, target current
I_off = 30;         % A, offset from for-loop initial
I_max = 90;         % A, max current in loop
for i = 1:I_max - I_off;
air(i) = f_air_init( ...
    1.1025, ...     % kg/m^3, rho
    0.7, ...        % Pr
    0.7, ...        % Pr_s
    17E-6, ...      % m^2/s, nu
    25, ...         % m/s, u_inlet
    0.01, ...       % m^2, A_x_inlet
    40, ...         % C, T_i
    0.027, ...      % W/mK, K
    1000);          % J/kgK, c_p
batt(i) = f_batt_init( ...
    i + I_off, ...  % A, I_batt
    18E-3, ...      % ohm, R_cell
    3.6, ...        % V, V_cell
    18, ...         % Wh, E_cell
    36, ...         % n_ser
    24, ...         % n_par
    24, ...         % n_wid
    20.45E-3, ...   % m, d_cell
    70.55E-3, ...   % m, l_cell
    23E-3, ...      % m, spac_wid
    23E-3, ...      % m, spac_len
    0, ...          % T/F, staggered
    50, ...         % C, T_s
    1E-6, ...       % m, e_pipe
    14E-3);         % m, h_pipe

for u_tube = 0.00:0.01:1E3
    air(i).u_inlet = u_tube;
    tube(i) = f_tube(batt(i), air(i));
    if abs(tube(i).T_batt_o - batt(i).T_s) < 1 && tube(i).T_batt_o < batt(i).T_s
        break;
    end
end
for u_pipe = 0.00:0.01:1E3
    air(i).u_inlet = u_pipe;
    pipe(i) = f_pipe(batt(i), air(i));
    if abs(pipe(i).T_batt_o - batt(i).T_s) < 1 && pipe(i).T_batt_o < batt(i).T_s
        break;
    end
end
end

%% PRINT TO TERMINAL

air(I_target - I_off)
batt(I_target - I_off)
pipe(I_target - I_off)
tube(I_target - I_off)

v_dot_tube_CFM = f_CFM(tube(I_target - I_off).v_dot_tube)
v_dot_pipe_CFM = f_CFM(pipe(I_target - I_off).v_dot_pipe)

%% PLOT

f_plot(batt, air, tube, pipe);

toc
end

%% FUNCTIONS FOR TUBE BANKS (CONSTANT TEMPERATURE HEAT TRANSFER)

function tube = f_tube(batt, air)
tube.A_x_batt = batt.l_cell * batt.width;
tube.u_batt = air.u_inlet .* (air.A_x_inlet ./ tube.A_x_batt);
tube.u_max_tube = f_u_max_tube(tube.u_batt, batt.S_T, batt.S_D, batt.d_cell);
tube.Re_D_max = f_Re(tube.u_max_tube, batt.d_cell, air.nu);
[tube.C_1, tube.m] = f_C_1(tube.Re_D_max, batt.stag, batt.S_T, batt.S_L);
tube.C_2 = f_C_2(batt.n_len, batt.stag);
tube.Nu_D_avg_tube = f_Nu_D_avg_tube(tube.C_1, tube.C_2, tube.Re_D_max, tube.m, air.Pr, air.Pr_s);
tube.h_tube = f_h(tube.Nu_D_avg_tube, air.k, batt.d_cell);

tube.v_dot_tube = tube.u_batt .* tube.A_x_batt;
tube.m_dot_tube = tube.v_dot_tube .* air.rho;
% tube.T_o = f_T_outlet_tube(batt.T_s, air.T_i, batt.d_cell, batt.n_cell, tube.h_tube, air.rho, tube.u_batt, batt.n_wid, batt.S_T, air.c_p);
tube.T_o = air.T_i + batt.q_dot_batt ./ (tube.m_dot_tube .* air.c_p);

tube.A_flux = batt.n_cell .* batt.d_cell .* pi .* batt.l_cell;
tube.q_flux = batt.q_dot_batt ./ tube.A_flux;
tube.T_batt_i = air.T_i + f_delta_T(tube.q_flux, tube.h_tube);
tube.T_batt_o = tube.T_o + f_delta_T(tube.q_flux, tube.h_tube);

% tube.delta_T_lm = f_delta_T_lm(batt.T_s, air.T_i, tube.T_o);
% tube.q_len_tube = f_q_len_tube(batt.n_cell, tube.h_tube, batt.d_cell, tube.delta_T_lm);

% tube.q_dot_tube = tube.q_len_tube * batt.l_cell;
% tube.m_dot_tube = f_m_dot(batt.q_dot_batt, air.c_p, (tube.T_o - air.T_i));
% tube.v_dot_tube = f_v_dot(tube.m_dot_tube, air.rho);

tube.P_L = batt.S_L / batt.d_cell;
tube.P_T = batt.S_T / batt.d_cell;
[tube.ff_tube, tube.xi_tube] = f_ff_xi_tube(batt.stag, tube.Re_D_max, tube.P_L, tube.P_T);
tube.delta_p_tube = f_delta_p_tube(batt.n_len, tube.xi_tube, air.rho, tube.u_max_tube, tube.ff_tube);
end

function Nu_D_avg_tube = f_Nu_D_avg_tube(C_1, C_2, Re_D_max, m, Pr, Pr_s) % BERGMAN 7.58, 7.59
Nu_D_avg_tube = C_1 .* C_2 .* (Re_D_max .^ m) .* (Pr .^ 0.36) .* ((Pr ./ Pr_s) .^ 0.25);
end

function [C_1, m] = f_C_1(Re, stag, S_T, S_L) % BERGMAN TABLE 7.5 AND 7.2 (KEEP 7.2 IN MIND FOR "SINGLE TUBES IN CROSSFLOW")
if stag == 0
    if      Re > 2E5
                            C_1 = 0.021;
                            m = 0.84;
    elseif  Re > 1E3
        if (S_T / S_L) > 0.7
                            C_1 = 0.27;
                            m = 0.63;
        else
                            C_1 = 0;
                            m = 0;
        end
    elseif  Re > 1E2
                            [C_1, m] = f_C_m_single_isolated(Re);
    elseif  Re > 10
                            C_1 = 0.80;
                            m = 0.40;
    else
                            C_1 = 0;
                            m = 0;
    end
else
    if      Re > 2E5
                            C_1 = 0.022;
                            m = 0.84;
    elseif  Re > 1E3
        if (S_T / S_L) < 2
                            C_1 = 0.35 .* ((S_T ./ S_L) .^ 0.2);
                            m = 0.60;
        else
                            C_1 = 0.40;
                            m = 0.60;
        end
    elseif  Re > 1E2
                            [C_1, m] = f_C_m_single_isolated(Re);
    elseif  Re > 10
                            C_1 = 0.90;
                            m = 0.40;
    else
                            C_1 = 0;
                            m = 0;
    end
end
end

function [C, m] = f_C_m_single_isolated(Re) % BERGMAN TABLE 7.2
lnRe = log(Re);
C = -0.0065 .* (lnRe .^ 2) + 0.0048 .* lnRe + 0.9957;
m = 0.0041 .* (lnRe .^ 2) - 0.0197 .* lnRe + 0.3676;
end

function C_2 = f_C_2(N_L, stag) % BERGMAN TABLE 7.6
if N_L < 20
    C_2 = ...
    ((-3E-5 .* N_L .^ 4) + (0.0014 .* N_L .^ 3) + (-0.0211 .* N_L .^ 2) + (0.1446 .* N_L) + 0.5791) - ...
    stag .* ...
    ((7E-6 .* N_L .^ 4) + (-0.0003 .* N_L .^ 3) + (0.0055 .* N_L .^ 2) + (-0.0384 .* N_L) + 0.0946);
else
    C_2 = 1;
end
end

function u_max_tube = f_u_max_tube(u, S_T, S_D, D) % BERGMAN 7.60, 7.61
if S_D <= ((S_T + D) / 2)
    u_max_tube = u .* S_T ./ (2 .* (S_D - D));
else
    u_max_tube = u .* S_T ./ (S_T - D);
end
end

function delta_T_lm = f_delta_T_lm(T_s, T_i, T_o) % BERGMAN 7.62
delta_T_lm = ((T_s - T_i) - (T_s - T_o)) ./ log((T_s - T_i) ./ (T_s - T_o));
end

function T_outlet_tube = f_T_outlet_tube(T_s, T_i, D, N, h_avg, rho, u, N_T, S_T, c_p) % BERGMAN 7.63
T_outlet_tube = -((exp((-pi .* D .* N .* h_avg) ./ (rho .* u .* N_T .* S_T .* c_p)) .* (T_s - T_i)) - T_s);
end

function q_len_tube = f_q_len_tube(N, h_avg, D, delta_T_lm) % BERGMAN 7.64
q_len_tube = N .* h_avg .* D .* delta_T_lm .* pi;
end

function delta_p_tube = f_delta_p_tube(N_L, xi, rho, u_max, ff) % BERGMAN 7.65
delta_p_tube = N_L .* xi .* ((rho .* u_max .^ 2) ./ 2) .* ff;
end

function [ff_tube, xi_tube] = f_ff_xi_tube(stag, Re_D_max, P_L, P_T) % BERGMAN FIGURE 7.14, 7.15
if stag == 0
    xi_tube = 1 .* ((P_T - 1) ./ (P_L - 1));
    ff_tube = 0.21 .* ((P_L - 1) .^ -0.755);
else
    xi_tube = 1;
    ff_tube = 0.43 .* ((P_T - 1) .^ -0.618);
end
end

%% FUNCTIONS FOR RECTANGULAR PIPE (CONSTANT FLUX HEAT TRANSFER)

function pipe = f_pipe(batt, air)
pipe.A_x_pipe = batt.h_pipe .* batt.width;
pipe.u_pipe = air.u_inlet .* (air.A_x_inlet ./ pipe.A_x_pipe);
pipe.Per_pipe = 2 .* batt.width + 2 .* batt.h_pipe;
pipe.D_h = f_D_h(pipe.A_x_pipe, pipe.Per_pipe);
pipe.Re_D_h = f_Re(pipe.u_pipe, pipe.D_h, air.nu);
pipe.ff = f_ff_pipe(batt.e_pipe, pipe.D_h, pipe.Re_D_h);
pipe.Nu_D_h = f_Nu_D_pipe(pipe.ff, pipe.Re_D_h, air.Pr);
pipe.h = f_h(pipe.Nu_D_h, air.k, pipe.D_h);
pipe.A_flux = ((batt.d_cell .^ 2) / 4) .* pi .* batt.n_cell;
pipe.q_flux = batt.q_dot_batt ./ pipe.A_flux;
pipe.delta_T = f_delta_T(pipe.q_flux, pipe.h);
pipe.v_dot_pipe = pipe.u_pipe .* pipe.A_x_pipe;
pipe.m_dot_pipe = pipe.v_dot_pipe .* air.rho;
pipe.T_o = ...
    air.T_i + ...
    batt.length .* (pipe.q_flux .* batt.width) ./ ...
    (pipe.m_dot_pipe .* air.c_p);
pipe.T_batt_i = air.T_i + f_delta_T(pipe.q_flux, pipe.h);
pipe.T_batt_o = pipe.T_o + f_delta_T(pipe.q_flux, pipe.h);
pipe.delta_p_pipe = f_delta_p_pipe(pipe.ff, air.rho, pipe.u_pipe, pipe.D_h, batt.length);
end

function Nu_D_pipe = f_Nu_D_pipe(ff, Re_D, Pr) % BERGMAN 8.60 OR 8.61 OR 8.62
Nu_D_pipe = ... % 8.62
    ((ff ./ 8) .* (Re_D - 1000) .* Pr) ./ ...
    (1 + (12.7 .* ((ff ./ 8) .^ 0.5) .* ((Pr .^ (2 / 3)) - 1)));
if Nu_D_pipe <= 0
    Nu_D_pipe = 0.027 .* (Re_D .^ 0.8) .* (Pr .^ 0.333); % 8.61
end
if Nu_D_pipe <= 0
    Nu_D_pipe = 0.023 .* (Re_D .^ 0.8) .* (Pr .^ 0.3); % 8.60
end
end

function D_h = f_D_h(A_c, P) % BERGMAN 8.66
D_h = 4 .* A_c ./ P;
end

function ff = f_ff_pipe(e, D, Re_D) % GERHART 8.35b
ff = (-1.8 .* log10(((e ./ (D .* 3.7)) .^ 1.11) + (6.9 ./ Re_D))) .^ -2;
end

function delta_p = f_delta_p_pipe(ff, rho, u_avg, D, x) % BERGMAN 8.22a
delta_p = ff .* rho .* u_avg .^ 2 ./ (2 .* D) .* x;
end

%% OTHER HEAT TRANSFER FUNCTIONS

function Re = f_Re(u, L, nu)
Re = u .* L ./ nu;
end

function h = f_h(Nu, k, L)
h = Nu .* k ./ L;
end

function q_flux = f_q_flux(h, delta_T)
% IF CONDUCTIVE, h = k/L...
q_flux_conv = h .* delta_T;
end

function delta_T = f_delta_T(q_flux, h)
% IF CONDUCTIVE, h = k/L...
delta_T = q_flux ./ h;
end

function q_dot = f_q_dot(m_dot, c_p, delta_T)
q_dot = m_dot .* c_p .* delta_T;
end

function m_dot = f_m_dot(q_dot, c_p, delta_T)
m_dot = q_dot ./ c_p ./ delta_T;
end

%% FUNCTIONS FOR PUMPS/FANS

function v_dot = f_v_dot(m_dot, rho)
v_dot = m_dot / rho;
end

function CFM = f_CFM(ms)
CFM = ms .* 2118.8799727597;
end

%% OTHER FUNCTIONS

function batt = f_batt_init(I_batt, R_cell, V_cell, E_cell, n_ser, n_par, n_wid, d_cell, l_cell, spac_wid, spac_len, stag, T_s, e_pipe, h_pipe)
batt.n_ser = n_ser;
batt.n_par = n_par;
batt.n_wid = n_wid;
batt.n_cell = n_ser * n_par;
batt.n_len = batt.n_cell / n_wid;
batt.width = n_wid * spac_wid;
batt.length = batt.n_len * spac_len;
batt.stag = stag;

batt.S_T = spac_wid;
batt.S_L = spac_len;
batt.S_D = (batt.S_T ^ 2 + batt.S_L ^ 2) ^ 0.5;

batt.d_cell = d_cell;
batt.l_cell = l_cell;

batt.I_batt = I_batt;
batt.I_cell = I_batt / n_par;
batt.R_cell = R_cell;
batt.V_cell = V_cell;
batt.V_batt = V_cell * n_ser;
batt.E_batt = E_cell * batt.n_cell;
batt.ESR = batt.n_ser * (batt.R_cell / batt.n_par);

batt.q_dot_cell = (batt.I_cell ^ 2) * batt.R_cell;
batt.q_dot_batt = (batt.I_batt ^ 2) * batt.ESR;

batt.T_s = T_s;
batt.e_pipe = e_pipe;
batt.h_pipe = h_pipe;
end

function air = f_air_init(rho, Pr, Pr_s, nu, u_inlet, A_x_inlet, T_i, k, c_p)
air.rho = rho;
air.Pr = Pr;
air.Pr_s = Pr_s;
air.nu = nu;
air.u_inlet = u_inlet;
air.A_x_inlet = A_x_inlet;
air.m_dot = air.u_inlet * air.A_x_inlet;
air.T_i = T_i;
air.k = k;
air.c_p = c_p;
end

function [] = f_plot(b, a, t, p)
tb = struct2table(b);
ta = struct2table(a);
tt = struct2table(t);
tp = struct2table(p);

figure
plot(tb.I_batt, tt.delta_p_tube ./ tp.delta_p_pipe);
grid on
title("Tube/pipe pressure drop ratio versus battery current")
xlabel("Current through battery (A)")
ylabel("Tube pressure / pipe pressure")

figure
semilogy(tb.I_batt, tt.delta_p_tube)
hold on
semilogy(tb.I_batt, tp.delta_p_pipe)
grid on
legend(["tube", "pipe"])
title("Pressure drop versus battery current")
xlabel("Current through battery (A)")
ylabel("Pressure drop across battery (Pa)")
hold off

figure
loglog(f_CFM(tt.v_dot_tube), tt.delta_p_tube)
hold on
loglog(f_CFM(tp.v_dot_pipe), tp.delta_p_pipe)
grid on
legend(["tube", "pipe"])
title("Pressure drop versus volumetric flow rate")
xlabel("Flow rate (CFM)")
ylabel("Pressure drop across battery (Pa)")
hold off

figure
loglog(tb.I_batt, f_CFM(tt.v_dot_tube))
hold on
loglog(tb.I_batt, f_CFM(tp.v_dot_pipe))
grid on
legend(["tube", "pipe"])
title("Volumetric flow rate versus battery current")
xlabel("Current through battery (A)")
ylabel("Flow rate (CFM)")
hold off

figure
plot(tb.I_batt, tt.C_1)
hold on
plot(tb.I_batt, tt.m)
grid on
legend(["C_1", "m"])
title("Tube bank C_1 and m constants versus battery current")
xlabel("Current through battery (A)")
ylabel("Magnitude")
hold off

figure
plot(tb.I_batt, tt.q_flux)
hold on
plot(tb.I_batt, tp.q_flux)
grid on
legend(["tube", "pipe"])
title("Heat flux versus battery current")
xlabel("Current through battery (A)")
ylabel("Heat flux (W/m^2)")
hold off

figure
plot(tb.I_batt, tt.h_tube)
hold on
plot(tb.I_batt, tp.h)
grid on
legend(["tube", "pipe"])
title("Heat transfer coefficient versus battery current")
xlabel("Current through battery (A)")
ylabel("h (W/m^2K)")
hold off

figure
plot(tb.I_batt, tt.T_batt_o - tt.T_o)
hold on
plot(tb.I_batt, tp.T_batt_o - tp.T_o)
grid on
legend(["tube", "pipe"])
title("Difference between battery and air temperature versus battery current")
xlabel("Current through battery (A)")
ylabel("h (W/m^2K)")
hold off

figure
plot3(tb.I_batt, f_CFM(tt.v_dot_tube), tt.delta_p_tube)
hold on
plot3(tb.I_batt, f_CFM(tp.v_dot_pipe), tp.delta_p_pipe)
plot3(tb.I_batt, f_CFM(tt.v_dot_tube), zeros(size(tt,1),1), "Color", [0 0 0])
plot3(tb.I_batt, f_CFM(tp.v_dot_pipe), zeros(size(tt,1),1), "Color", [0 0 0])
grid on
legend(["tube", "pipe"])
title("Pressure drop and volumetric flow rate versus battery current")
xlabel("Current (A)")
ylabel("Flow rate (CFM)")
zlabel("Pressure drop (Pa)")
hold off

end










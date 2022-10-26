%Paul Adkisson
%4/12/21
%Purpose: ODEs for Rattay Cable Equations
%   Inputs:
%       state_vars = [V(k), n_gate(k), m_gate(k), h_gate(k)] where k = 1, 2, ..., N
%           are spatial indices
%       V_e(k) external voltages at each spatial index k
%       I_int(k) internal current at each spatial index k
%   Outputs:
%       [dV/dt(k), dn/dt(k), dm/dt(k), dh/dt(k)]

function dydt = rattay_odes(state_vars, V_e, I_int)
    %Setup 
    load("rattay_constants.mat");
    V = state_vars(:, 1);
    n_gate = state_vars(:, 2);
    m_gate = state_vars(:, 3);
    h_gate = state_vars(:, 4);
    dVdt = zeros(length(V), 1);
    
    %I_ionic
    I_ionic = pi*d*L*(g_na*m_gate.^3.*h_gate.*(V-E_na) + g_k*n_gate.^4.*(V-E_k) + g_l*(V-E_L)) - I_int;
    %I_ionic = pi*d*L*(g_l*(V-V_rest)) - I_int; %cable lif
    
    %Boundary Conditions
    dVdt(1) = 1/C_m*( Ga*(V(2) - V(1) + V_e(2) - V_e(1)) - I_ionic(1));
    dVdt(end) = 1/C_m*( Ga*(V(end-1) - V(end) + V_e(end-1) - V_e(end)) - I_ionic(end));
    
    %Non-Boundary
    dVdt(2:end-1) = 1/C_m*( Ga*(V(1:end-2) - 2*V(2:end-1) + V(3:end) + V_e(1:end-2) - 2*V_e(2:end-1) + V_e(3:end)) ...
        - I_ionic(2:end-1));
    
    %Gate Eqs
    dndt = alpha_n(V).*(1-n_gate) - beta_n(V).*n_gate;
    dmdt = alpha_m(V).*(1-m_gate) - beta_m(V).*m_gate;
    dhdt = alpha_h(V).*(1-h_gate) - beta_h(V).*h_gate;
    
    %Return
    dydt = [dVdt, dndt, dmdt, dhdt];
end
%%% Paul Adkisson
%%% 7/25/22
%%% Purpose: Define Parameters for Rattay Model Cable Eq for variable
%%% diameter z

function rattay_z_constants(z)
    %Extracellular Electrode Constants
    N = 25; %Number of nodes
    n = (-(N-1)/2:(N-1)/2)'; %node numbers
    D = 10*10^(-4); %Fiber Diameter (cm)
    d = 0.7*D; %Axon Diameter (cm)
    delta_x = 100*D; %Internode Distance (cm)
    x = delta_x*n; %node positions
    %z = 0.1; %electrode height (cm)
    r = sqrt(x.^2 + z^2); %distance from electrode to node (cm)
    L = 2.5*10^(-4); %node length (cm)
    rho_i = 0.055; %intracellular resistivity (kOhm-cm)
    rho_e = 0.3; %extracellular resistivity (kOhm-cm)
    Ga = pi*d^2 / (4*rho_i*delta_x); %axial conductivity (mS)
    c_m = 1; %uF/cm^2
    C_m = c_m*pi*d*L; %uF

    %HH Constants
    R = 8.31451; %J/(mol-K)
    F = 96485.3; %C/mol
    T_C = 37; %Celsius
    T = 273.16 + T_C; %K

    g_na = 120; %mSiemann/cm^2
    g_k = 36; %mSiemann/cm^2
    g_l = 0.3; %mSiemann/cm^2

    Na_in = 12; %mM
    Na_out = 145; %mM
    E_na = -R*T/F * 10^3 * log(Na_in/Na_out); %mV
    K_in = 140; %mM
    K_out = 5; %mM
    E_k = -R*T/F * 10^3 * log(K_in/K_out); %mV
    Cl_in = 4; %mM
    Cl_out = 110; %mM

    P_K = 1;
    P_Na = 0.03;
    P_Cl = 0.1;
    V_rest = R*T/F*10^3 * log( (P_K*K_out + P_Na*Na_out + P_Cl*Cl_in) / (P_K*K_in + P_Na*Na_in + P_Cl*Cl_out) ); %mV

    epsilon = 0.01; %mV (removing singularities with Taylor Expansion)
    alpha_n = @(V) (abs(V+55)>=epsilon) .* 0.01.*(V+55) ./ (1-exp(-(V+55)/10)) + (abs(V+55)<epsilon) .* (0.1 + (V+55)/200 + (V+55).^2/12000);
    beta_n = @(V) 0.125*exp(-(V+65)/80);
    alpha_m = @(V) (abs(V+40)>=epsilon) .* 0.1.*(V+40) ./ (1-exp(-(V+40)/10)) + (abs(V+40)<epsilon) .* (1 + (V+40)/20 + (V+40).^2/1200);
    beta_m = @(V) 4*exp(-(V+65)/18);
    alpha_h = @(V) 0.07*exp(-(V+65)/20);
    beta_h = @(V) 1 ./ (1 + exp(-(V+35)/10));

    tau_n = @(V) 1 / (alpha_n(V) + beta_n(V));
    n_inf = @(V) alpha_n(V) * tau_n(V);
    tau_m = @(V) 1 / (alpha_m(V) + beta_m(V));
    m_inf = @(V) alpha_m(V) * tau_m(V);
    tau_h = @(V) 1 / (alpha_h(V) + beta_h(V));
    h_inf = @(V) alpha_h(V) * tau_h(V);

    I_k_rest = g_k*n_inf(V_rest)^4*(V_rest - E_k);
    I_na_rest = g_na*m_inf(V_rest)^3*h_inf(V_rest)*(V_rest - E_na);
    I_l_rest = -(I_k_rest + I_na_rest);
    E_L = V_rest-I_l_rest/g_l; %mV 

    save rattay_constants.mat
end
function electric_r = mirror_est(z)
    rho_e = 3; %extracellular resistivity (Ohm-m)
    N = 25; %Number of nodes
    n = (-(N-1)/2:(N-1)/2)'; %node numbers
    D = 10*10^(-6); %Fiber Diameter (m)
    delta_x = 100*D; %Internode Distance (m)
    x = delta_x*n; %node positions
    r = sqrt(x.^2 + z.^2); %distance from electrode to node (cm)
    gL = 25*1e-9;
    V_e = rho_e ./ (4*pi*z);
    V_es = rho_e ./ (4*pi*r);
    V_e_bar = mean(V_es, 1);
    electric_r = (V_e_bar - V_e)*gL; %electric_r = Vm_ss/I_el
end
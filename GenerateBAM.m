function GenerateBAM(brains, N_E, N_I, f, p, w_plus, w_minus, w, sim_path)
    N = N_E + N_I;
    adja = w*ones(N, N);
    num_group = floor(f*N_E);
    for i = 0:(p-1)
        group_idx = (i*num_group+1):(i+1)*num_group;
        adja(group_idx, 1:p*num_group) = w_minus;
        adja(group_idx, group_idx) = w_plus;
        adja(p*num_group+1:N_E, 1:p*num_group) = w_minus;
    end
    diag_idx = logical(eye(size(adja)));
    adja(diag_idx) = 0; %Disallow self-connections
    savepath = strcat(sim_path, "/adja.mat");
    save(savepath, "adja")
    
    min_x = sqrt(50*1e-12); %minimum distance of 10um
    max_x = 2*1e-3; %From Levitt et al.
    min_y = sqrt(50*1e-12);
    max_y = 200*1e-6; %From Levitt et al.
    regular_x = (max_x - min_x) / num_group;
    regular_y = (max_y - min_y) / num_group;
    group_idx = 0:num_group-1;
    for brain = brains
        rng(brain);
        ball_x = min_x + regular_x*group_idx+regular_x.*rand(1, num_group);
        ball_y = min_y + regular_y*group_idx+regular_y.*rand(1, num_group);
        ball_r = sqrt(ball_x.^2 + ball_y.^2);
        electric_r = mirror_est(ball_r);
        
        brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
        mkdir(brainpath)
        save(strcat(brainpath, "/r.mat"), "ball_r", "electric_r")
    end
end

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
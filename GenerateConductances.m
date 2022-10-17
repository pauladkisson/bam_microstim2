function GenerateConductances(N_E, N_I, sim_path)
    N = N_E + N_I;
    % Synaptic Conductance = [pyramidal, interneuron]
    G_ampa_rec = [0.05, 0.04]*1e-9 .* 1600 / N_E; %nS
    G_nmda = [0.165, 0.13]*1e-9 .* 1600 / N_E; %nS
    G_gaba = [1.3, 1]*1e-9 .* 400 / N_I; %nS
    AMPA = zeros(N, N);
    AMPA(1:N_E, 1:N_E) = G_ampa_rec(1);
    AMPA(1:N_E, N_E+1:end) = G_ampa_rec(2);
    NMDA = zeros(N, N);
    NMDA(1:N_E, 1:N_E) = G_nmda(1);
    NMDA(1:N_E, N_E+1:end) = G_nmda(2);
    GABA = zeros(N, N);
    GABA(N_E+1:end, 1:N_E) = G_gaba(1);
    GABA(N_E+1:end, N_E+1:end) = G_gaba(2);
    
    %{
    figure;
    imagesc(AMPA)
    colorbar;
    title("AMPA")
    
    figure;
    imagesc(NMDA)
    colorbar;
    title("NMDA")
    
    figure;
    imagesc(GABA)
    colorbar;
    title("GABA")
    %}
    save_path = strcat(sim_path, "/conductances.mat");
    save(save_path, 'AMPA', 'NMDA', 'GABA')
end
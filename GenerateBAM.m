function GenerateBAM(N_E, N_I, f, p, w_plus, w_minus, w, sim_path)
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
end
function GenerateBAM_NoRecX(brains, N_E, N_I, f, p, w_plus, w_minus, w, sim_path)
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
    adja(1:N_E, 1:N_E) = 0; %No Recurrent Excitations
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
        axon_SA = 100*(1e-6^2); %Typical axon surface area = 100um^2
        electric_r = -axon_SA ./ (4*pi*ball_r.^2);
        brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
        mkdir(brainpath)
        save(strcat(brainpath, "/r.mat"), "ball_r", "electric_r")
        
        %{
        %figure;
        %imagesc(adja)
        %colorbar;
        pulse_test_stim = -10*1e-6;
        figure;
        subplot(2, 1, 1)
        scatter(ball_r*1e6, electric_r*pulse_test_stim*1e9)
        title("Pulse Stimulation")
        xlabel("Distance from Electrode to Neuron (um)")
        ylabel("Stimulation Current (nA)")
        subplot(2, 1, 2)
        histogram(electric_r*pulse_test_stim*1e9, 48)
        xlabel("Stimulation Current (nA)")
        ylabel("Number of Neurons")

        galvanic_test_stim = -100*1e-9;
        figure;
        subplot(2, 1, 1)
        scatter(ball_r*1e6, electric_r*galvanic_test_stim*1e12)
        title("Galvanic Stimulation")
        xlabel("Distance from Electrode to Neuron (um)")
        ylabel("Stimulation Current (pA)")
        subplot(2, 1, 2)
        histogram(electric_r*galvanic_test_stim*1e12, 24)
        xlabel("Stimulation Current (pA)")
        ylabel("Number of Neurons")
        %}
    end
end
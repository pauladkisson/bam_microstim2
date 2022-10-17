%%% sandbox for FR spectral analysis

%% Figure 1a
%First Replicate Results from Matzner & Bar-Gad (2015)
f0 = 12; %oscillation frequency (Hz)
r0 = 5; %Average FR (spk/s)
m = 0.5; %Modulation index (strength of oscillation)
T = 1; %Trial duration (s)
dt = 0.001; % (ms)
t = 0:dt:T; %Simulation time (s)
fr = r0 * (1 + m*cos(2*pi*f0*t)); %Instaneous Firing Rate

figure;
plot(t, fr, 'k-', 'Linewidth', 5)
xlabel("Time (s)")
ylabel("Instantaneous Firing Rate (spk/s)")
ylim([0, 10])
title("Figure 1a")
set(gca, 'Fontsize', 24)

%% Figure 1b
T = 60*5;
t = 0:dt:T; %Simulation time (s)
fr = r0 * (1 + m*cos(2*pi*f0*t)); %Instaneous Firing Rate
window = 1000; %size of welch window (# of samples)
n_overlap = 1; %number of overlapping points in the window
nyquist = 1 / (2*dt); %Nyquist Frequency
fr_norm = fr - mean(fr);
[psd, freqs] = pwelch(fr_norm, window, n_overlap, window);
freqs = freqs * nyquist / pi;
disp(max(freqs))
disp(min(freqs))
figure;
plot(freqs, psd, 'k-', 'Linewidth', 5)
xline(f0, 'r--', 'Linewidth', 2)
xlim([0, 100])
title("Figure 1b")
xlabel("Frequency (Hz)")
ylabel("Power ( (spk/s)^2/Hz )")
set(gca, 'Fontsize', 24)

%% Figure 1c
rng(1);
spikes = double(rand(size(fr)) < dt*fr);
figure;
plot(t(t<=1), spikes(t<=1), 'k-', 'Linewidth', 5)
xlabel("Time (s)")
ylabel("Spike")
title("Figure 1c")
yticks([0, 1])
yticklabels([0, 1])
set(gca, 'Fontsize', 24)

%% Figure 1d
norm_spikes = spikes - mean(spikes);
[psd, freqs] = pwelch(spikes, window, n_overlap, window);
freqs = freqs * nyquist / pi;
figure;
plot(freqs, psd, 'k-', 'Linewidth', 5)
xline(f0, 'r--', 'Linewidth', 2)
xlim([0, 100])
xlabel("Frequency (Hz)")
ylabel("Power (A.U.)")
set(gca, 'Fontsize', 24)
title("Figure 1d")

%% Figures 1E, F, G
r0s = [5, 15, 30];
titles = ["Figure 1E", "Figure 1F", "Figure 1G"];
for i = 1:length(r0s)
    r0 = r0s(i);
    fr = r0 * (1 + m*cos(2*pi*f0*t));
    spikes = double(rand(size(fr)) < dt*fr);
    [psd, freqs] = pwelch(spikes, window, n_overlap, window);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    figure;
    plot(freqs, psd_snr, 'k-', 'Linewidth', 5)
    xline(f0, 'r--', 'Linewidth', 2)
    xlim([0, 100])
    ylim([-5, 25])
    xlabel("Frequency (Hz)")
    ylabel("SNR")
    title(titles(i))
    set(gca, 'Fontsize', 24)
end

%% Figure 4B with ground truth f0
T = 1*60;
fs = 1/dt;
window = ones(1000, 1);
win_len = length(window);
t = 0:dt:T; %Simulation time (s)
r0s = 1:5:80;
ms = [0.1, 0.25, 0.4, 0.5, 0.75, 1];
r0_hats = zeros(length(r0s), length(ms));
mhats = zeros(length(r0s), length(ms));
mhat_stds = zeros(length(r0s), length(ms));
num_trials = 10;
for i = 1:length(ms)
    m = ms(i);
    fprintf("m=%0.1f \n", m)
    for j = 1:length(r0s)
        r0 = r0s(j);
        fr = r0 * (1 + m*cos(2*pi*f0*t));
        trial_mhats = zeros(num_trials, 1);
        for k = 1:num_trials
            spikes = double(rand(size(fr)) < dt*fr);
            [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
            freqs = freqs * nyquist / pi;
            psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
            [~, f0_ind] = min(abs(freqs-f0));
            peak_snr = psd_snr(f0_ind);
            r0_hat = sum(spikes) / T;
            mhat = get_mhat(peak_snr, r0_hat, T, fs, win_len);
            trial_mhats(k) = mhat;
        end
        mhats(j, i) = mean(trial_mhats);
        mhat_stds(j, i) = std(trial_mhats);
        r0_hats(j, i) = r0_hat;
    end
end
%% Figure 4B Plot
default_colors = get(gca, 'Colororder');
figure;
hold on
for i = 1:length(ms)
    errorbar(r0s, mhats(:, i), mhat_stds(:, i), 'Linewidth', 5)
    yline(ms(i), '--', 'Color', default_colors(i, :), 'Linewidth', 2)
end
hold off
xlabel("Baseline Firing Rate r0 (spk/s)")
ylabel("Estimated Modulation Index $\hat{m}$", 'Interpreter', 'LaTeX')
title("Figure 4B")
set(gca, 'Fontsize', 24)

%% Figure 4B with NEO-based spike detection of f0
thresh_factor = 15;
T = 1*60;
fs = 1/dt;
window = ones(1000, 1);
win_len = length(window);
t = 0:dt:T; %Simulation time (s)
r0s = 1:5:80;
ms = [0.1, 0.25, 0.4, 0.5, 0.75, 1];
r0_hats = zeros(length(r0s), length(ms));
mhats = zeros(length(r0s), length(ms));
mhat_stds = zeros(length(r0s), length(ms));
num_trials = 10;
for i = 1:length(ms)
    m = ms(i);
    fprintf("m=%0.1f \n", m)
    for j = 1:length(r0s)
        r0 = r0s(j);
        fr = r0 * (1 + m*cos(2*pi*f0*t));
        trial_mhats = zeros(num_trials, 1);
        for k = 1:num_trials
            spikes = double(rand(size(fr)) < dt*fr);
            [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
            freqs = freqs * nyquist / pi;
            psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
            psd_snr = psd_snr(freqs>=1&freqs<=100); %omit baseline shift artifact
            freqs = freqs(freqs>=1&freqs<=100);
            ne_snr = neo(psd_snr);
            sigma = median(ne_snr) / 0.67;
            threshold = sigma*thresh_factor;
            if any(ne_snr>=threshold)
                [peak_snr, f0hat] = max(psd_snr(ne_snr>=threshold));
            else
                peak_snr = 0;
            end
            r0_hat = sum(spikes) / T;
            mhat = get_mhat(peak_snr, r0_hat, T, fs, win_len);
            trial_mhats(k) = mhat;
        end
        mhats(j, i) = mean(trial_mhats);
        mhat_stds(j, i) = std(trial_mhats);
        r0_hats(j, i) = r0_hat;
    end
end
%% Figure 4B NEO Plot
default_colors = get(gca, 'Colororder');
figure;
hold on
for i = 1:length(ms)
    errorbar(r0s, mhats(:, i), mhat_stds(:, i), 'Linewidth', 5)
    yline(ms(i), '--', 'Color', default_colors(i, :), 'Linewidth', 2)
end
hold off
xlabel("Baseline Firing Rate r0 (spk/s)")
ylabel("Estimated Modulation Index $\hat{m}$", 'Interpreter', 'LaTeX')
title("Figure 4B with NEO detection")
set(gca, 'Fontsize', 24)

%% Apply to Poisson input
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC I_b100/bam_constants.mat')
num_trials = 36;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 99);
for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    trialname = sprintf('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC I_b100/spikes/c=0.000/trial%0.0f/input.mat', ...
        trial);
    load(trialname, 'spikes');
    tstart = 0;
    tend = 1;
    T = tend - tstart;
    spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E);
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=1), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(ne_snr) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
%% Plot Poisson
figure;
scatter(f0hats, mhats, 100, 'k')
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
set(gca, 'Fontsize', 24)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1), 'k', 'Linewidth', 2)
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")
set(gca, 'Fontsize', 24)

%% Apply to Disconnected spikes (0-1s)
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC Disconnected/bam_constants.mat')
num_trials = 36;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 99);
for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    trialname = sprintf('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC Disconnected/brain1/data/0.0nA_galvanic/c=0.000/trial%0.0f.mat', ...
        trial);
    load(trialname, 'recspikes');
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    tstart = 1;
    tend = 3;
    T = tend - tstart;
    %spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E); %NS neurons
    %spikes = spikes(t>=tstart&t<=tend, 1:num_group); %P1
    spikes = spikes(t>=tstart&t<=tend, :); %All neurons
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=1), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(ne_snr) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
%% Plot Disconnected
figure;
scatter(f0hats, mhats, 100, 'k')
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
set(gca, 'Fontsize', 24)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1), 'k', 'Linewidth', 2)
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")
set(gca, 'Fontsize', 24)

%% Apply to Connected spikes (0-1s)
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC I_b100/bam_constants.mat')
num_trials = 36;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 99);
for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    trialname = sprintf('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC I_b100/brain1/data/0.0nA_galvanic/c=0.000/trial%0.0f.mat', ...
        trial);
    load(trialname, 'recspikes');
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    tstart = 1;
    tend = 3;
    T = tend - tstart;
    %spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E); %NS neurons
    %spikes = spikes(t>=tstart&t<=tend, 1:num_group); %P1
    spikes = spikes(t>=tstart&t<=tend, :); %All neurons
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=1), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(ne_snr) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
%% Plot Connected
figure;
scatter(f0hats, mhats, 100, 'k')
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
set(gca, 'Fontsize', 24)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1), 'k', 'Linewidth', 2)
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")
set(gca, 'Fontsize', 24)

%% Apply to Connected spikes (Baseline Activity) Long Duration (30s)
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation LongDuration/bam_constants.mat')
num_trials = 10;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 495);
for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    trialname = sprintf('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation LongDuration/brain1/data/0.0nA_galvanic/c=0.000/trial%0.0f.mat', ...
        trial);
    load(trialname, 'recspikes');
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    tstart = 0;
    tend = 30;
    T = tend - tstart;
    %spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E); %NS neurons
    %spikes = spikes(t>=tstart&t<=tend, 1:num_group); %P1
    spikes = spikes(t>=tstart&t<=tend, :); %All neurons
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=5), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(abs(ne_snr)) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
%% Plot Long Duration
figure;
scatter(f0hats, mhats, 100, 'k')
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
set(gca, 'Fontsize', 24)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1), 'k', 'Linewidth', 1)
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")
set(gca, 'Fontsize', 24)

%% Create Oscillatory input for more dramatic oscillations
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation EMBC I_b100/bam_constants.mat')
num_trials = 36;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 99);
r0 = 2400;
m = 0.5;
f0 = 40;
fr = r0 * (1 + m*cos(2*pi*f0*t)); %Instaneous Firing Rate

for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    spikes = rand(length(t), N) < (dt*fr)';
    tstart = 0;
    tend = 1;
    T = tend - tstart;
    %spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E); %NS neurons
    %spikes = spikes(t>=tstart&t<=tend, 1:num_group); %P1
    spikes = spikes(t>=tstart&t<=tend, :); %All neurons
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=1), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(ne_snr) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
%% Plot Poisson Gamma
figure;
hold on
scatter(f0hats, mhats, 100, 'k')
scatter(f0, m, 100, 'r')
xlim([39, 41])
ylim([0.4, 0.6])
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
set(gca, 'Fontsize', 24)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1), 'k', 'Linewidth', 1)
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")
set(gca, 'Fontsize', 24)

%% Manually created oscillations in disconnected neurons
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation Test/bam_constants.mat')
num_trials = 1;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 99);

for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    trialname = sprintf('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation Test/brain1/data/-1000.0nA_galvanic/c=0.000/trial%0.0f.mat', ...
        trial);
    load(trialname, 'recspikes');
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    tstart = 0;
    tend = 1;
    T = tend - tstart;
    %spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E); %NS neurons
    spikes = spikes(t>=tstart&t<=tend, 1:num_group); %P1
    %spikes = spikes(t>=tstart&t<=tend, :); %All neurons
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=1), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(ne_snr) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
figure;
scatter(f0hats, mhats)
ylim([0, 1])
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
median(mhats)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1))
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")

%% Manually created oscillations in Connected neurons
clear;
load('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation FeedforwardGamma/bam_constants.mat')
num_trials = 28;
mhats = zeros(num_trials, 1);
f0hats = zeros(num_trials, 1);
psd_snrs = zeros(num_trials, 99);

for trial = 1:num_trials
    num_group = 120;
    fprintf("Trial: %0.0f \n", trial)
    trialname = sprintf('/Volumes/Paul''s Data/Modeling/bam_microstim/Simulation FeedforwardGamma/brain1/data/-1400.0nA_galvanic/c=0.000/trial%0.0f.mat', ...
        trial);
    load(trialname, 'recspikes');
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    tstart = 3;
    tend = 4;
    T = tend - tstart;
    %spikes = spikes(t>=tstart&t<=tend, 2*num_group+1:N_E); %NS neurons
    spikes = spikes(t>=tstart&t<=tend, 1:num_group); %P1
    %spikes = spikes(t>=tstart&t<=tend, :); %All neurons
    spikes = sum(spikes, 2); %aggregate over all neurons for better resolution
    num_group = 1;
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    thresh_factor = 15;
    window = ones(sum(t<=1), 1);
    win_len = length(window);
    n_overlap = 0;
    [psd, freqs] = pwelch(spikes, window, n_overlap, win_len);
    freqs = freqs * nyquist / pi;
    psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
    psd_snr = psd_snr(freqs>=1&freqs<=100, :); %omit baseline shift artifact
    freqs = freqs(freqs>=1&freqs<=100);
    ne_snr = neo(psd_snr);
    sigma = median(ne_snr) / 0.67; %estimation of std
    threshold = sigma*thresh_factor;
    if any(ne_snr>=threshold)
        peak_snr = max(psd_snr(ne_snr>=threshold));
        f0hat = freqs(psd_snr==peak_snr);
        r0hat = sum(spikes) / T; 
        mhat = get_mhat(peak_snr, r0hat, T, fs, win_len);
        if length(f0hat)~=1 %multiple modulation frequencies shouldn't be possible
            f0hat = 0;
            mhat = 0;
        end
    else
        f0hat = 0;
        mhat = 0;
    end
    f0hats(trial) = f0hat;
    mhats(trial) = mhat;
    psd_snrs(trial, :) = psd_snr;
end
figure;
scatter(f0hats, mhats)
xlabel("Estimated Modulation Frequency $\hat{f_0}$ (Hz)", 'Interpreter', 'LaTeX')
ylabel("Estimated Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
median(mhats)

figure;
errorbar(freqs, mean(psd_snrs, 1), std(psd_snrs, [], 1))
xlabel("Frequency (Hz)")
ylabel("SNR")
title("Power Spectral Density Averaged over trials")

%% Functions
function mhat = get_mhat(peak_snr, r0, T, fs, win_len)
    N_wins = ceil(T*fs / win_len);
    peak_snr_cor = peak_snr*sqrt(N_wins); %welch correction
    mhat = 2*sqrt(peak_snr_cor / (r0*T));
end

function psi = neo(x)
    psi = zeros(size(x));
    psi(1, :) = x(1, :).^2 - x(1, :).*x(2, :);
    psi(end, :) = x(end).^2 - x(end, :).*x(end-1);
    psi(2:end-1, :) = x(2:end-1, :).^2 - x(3:end, :).*x(1:end-2, :);
end
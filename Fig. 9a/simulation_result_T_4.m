clc;
clearvars;
close all;
load("DFT_2.mat")
DFT_SIM = DFT_2; % The response of SIM, which has been optimized by using the gradient descent algorithm
N_x = 2; % Number of antennas in the x-direction
N_y = 2; % Number of antennas in the y-direction
N = N_x * N_y; % Number of antennas
SNR_set = 10 : 1 : 30; % SNR
for n_SNR = 1 : length(SNR_set)
    noise_var = 1 / 10 ^ (SNR_set(n_SNR) / 10); % noise variance
    tic
    rng(1);
    for n_MonteCarlo = 1 : 1000
        psi_x = 1.4 * rand - 0.7; % Elesctrical angle \psi_x in the x direction
        psi_y = 1.4 * rand - 0.7; % Elesctrical angle \psi_y in the y direction
        SV_x = exp(1i * pi * (0 : N_x - 1).' * psi_x);
        SV_y = exp(1i * pi * (0 : N_y - 1).' * psi_y);
        SV = kron(SV_y, SV_x); % Calculate the steering vector
        N_snapshot = 4; % Number of snapshots for each discrete bin
        Spectrum = zeros(N_snapshot * N_x, N_snapshot * N_y); % Fine-resolution angle spectrum
        p = zeros(N_snapshot, N_snapshot); % Output gain of the coarse-resolution angle spectrum
        index_p = zeros(N_snapshot, N_snapshot); % Corresponding index
        %% Find the optimal configuration of the input layer
        for ii = 1 : N_snapshot
            for jj = 1 : N_snapshot
                text_x = (ii - 1) / N_snapshot;
                text_y = (jj - 1) / N_snapshot;
                offset_x = exp(-1i / N_x * 2 * pi * (0 : N_x - 1).' * text_x);
                offset_y = exp(-1i / N_y * 2 * pi * (0 : N_y - 1).' * text_y);
                offset = kron(offset_y, offset_x); % Transmission coefficients of the input layer
                output_with_noise = DFT_SIM * diag(offset) * SV + ...
                    sqrt(noise_var / 2) * (randn(N, 1) + 1i * randn(N, 1)); % DFT in the wave domain
                [p(ii, jj),index_p(ii,jj)] = max(abs(output_with_noise));
                spectrum_coarse = reshape(output_with_noise, N_x, N_y); % Coarse-resolution spectrum in the ii-th slot of the jj-th block
                Spectrum(ii : N_snapshot : N_snapshot * N_x, jj : N_snapshot : N_snapshot * N_y) = spectrum_coarse;
            end
        end
        Spectrum = Spectrum.';
        p_peak = max(max(p)); % The peak of the angle spectrum
        [row, column] = find(p == p_peak); % The index corresponding to the peak
        index = index_p(row, column);
        y_index = ceil(index / N_x);
        x_index = mod(index - 1, N_x) + 1;
        psi_x_est = mod(((x_index - 1) / N_x + (row - 1) / N_snapshot /N_x) * 2 + 1, 2) - 1; % Estimated elesctrical angle \psi_x in the x direction
        psi_y_est = mod(((y_index - 1) / N_y + (column - 1) / N_snapshot / N_y) * 2 + 1, 2) - 1; % Estimated elesctrical angle \psi_y in the y direction
        mse_x_cache(n_MonteCarlo) = real((psi_x - psi_x_est) ^ 2);
        mse_y_cache(n_MonteCarlo) = real((psi_y - psi_y_est) ^ 2);
    end
    mse_x(n_SNR) = mean(mse_x_cache);
    mse_y(n_SNR) = mean(mse_y_cache);
    toc
end
figure;
semilogy(mse_x)
hold on
semilogy(mse_y)
mse_x_4_SIM = mse_x;
mse_y_4_SIM = mse_y;
save mse_x_4_SIM mse_x_4_SIM
save mse_y_4_SIM mse_y_4_SIM
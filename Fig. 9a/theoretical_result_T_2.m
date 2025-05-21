clc;
clearvars;
close all;
SNR_set = 10 : 1 : 30; % SNR
N_snapshot = 2; % Number of snapshots for each discrete bin
N_x = 2; % Number of antennas in the x-direction
N_y = 2; % Number of antennas in the y-direction
N = N_x * N_y; % Number of antennas
DFT_x_real = exp(-1i / N_x * 2 * pi * (0 : N_x - 1).' * (0 : N_x - 1));
DFT_y_real = exp(-1i / N_y * 2 * pi * (0 : N_y - 1).' * (0 : N_y - 1));
DFT = kron(DFT_y_real, DFT_x_real); % Generation of the 2D-DFT matrix
for n_SNR = 1 : length(SNR_set)
    tic
    noise_var = 1 / 10 ^ (SNR_set(n_SNR) / 10); % noise variance
    rng(1);
    for n_MonteCarlo = 1 : 1000
        psi_x = 1.4 * rand - 0.7; % Elesctrical angle \psi_x in the x direction
        psi_y = 1.4 * rand - 0.7; % Elesctrical angle \psi_y in the y direction
        SV_x = exp(1i * pi * (0 : N_x - 1).' * psi_x);
        SV_y = exp(1i * pi * (0 : N_y - 1).' * psi_y);
        SV = kron(SV_y, SV_x); % Calculate the steering vector
        for ii = 1 : N_snapshot
            for jj = 1 : N_snapshot
                text_x = (ii - 1) / N_snapshot;
                text_y = (jj - 1) / N_snapshot;
                offset_x = exp(-1i / N_x * 2 * pi * (0 : N_x - 1).' * text_x);
                offset_y = exp(-1i / N_y * 2 * pi * (0 : N_y - 1).' * text_y);
                offset = kron(offset_y, offset_x); % Transmission coefficients of the input layer
                output = DFT * diag(offset) * SV; % Output of the angle spectrum
                spectrum_coarse = reshape(output, N_x, N_y); % Coarse-resolution spectrum in the ii-th slot of the jj-th block
                Spectrum(ii : N_snapshot : N_snapshot * N_x, jj : N_snapshot : N_snapshot * N_y) = spectrum_coarse;
            end
        end
        Spectrum_peak = max(max(abs(Spectrum))); % Spectrum peak
        [row_peak, column_peak] = find(abs(Spectrum) == Spectrum_peak); % Corresponding index
        for ii = 1 : N * N_snapshot * N_snapshot
            index_x_ii = mod(ii - 1, N_x * N_snapshot) + 1;
            index_y_ii = ceil(ii / N_x / N_snapshot);
            if ii == (column_peak - 1) * N_x * N_snapshot + row_peak
                p = 1;
            else
                x1 = 2 * abs(real(Spectrum(index_x_ii, index_y_ii))) ^ 2 / noise_var;
                x2 = 2 * abs(imag(Spectrum(index_x_ii, index_y_ii))) ^ 2 / noise_var;
                x3 = 2 * abs(real(Spectrum(row_peak, column_peak))) ^ 2 / noise_var;
                x4 = 2 * abs(imag(Spectrum(row_peak, column_peak))) ^ 2 / noise_var;
                a = 1 * (1 + x1) + 1 * (1 + x2) - 1 * (1 + x3) - 1 * (1 + x4);
                b = 1 ^ 2 * (1 + x1 * 2) + 1 ^ 2 * (1 + x2 * 2) + 1 ^ 2 * (1 + x3 * 2) + 1 ^ 2 * (1 + x4 * 2);
                c = 1 ^ 3 * (1 + x1 * 3) + 1 ^ 3 * (1 + x2 * 3) - 1 ^ 3 * (1 + x3 * 3) - 1 ^ 3 * (1 + x4 * 3);
                if c == 0
                    p = 1;
                else
                    h = b ^ 3 / c ^ 2; % Eq. (38)
                    y = - a * sqrt(h / b) + h; % Eq. (39)
                    p = qfunc(((y / h) ^ (1/3) - 1 + 2 / 9 / h) / sqrt( 2 / 9 / h));
                end
            end
            Prob(ii) = p;
            psi_x_est(ii) = mod(((index_x_ii - 1) / N_snapshot / N_x) * 2 + 1, 2) - 1; % Estimated elesctrical angle \psi_x in the x direction
            psi_y_est(ii) = mod(((index_y_ii - 1) / N_snapshot / N_y) * 2 + 1, 2) - 1; % Estimated elesctrical angle \psi_y in the y direction
            error_x(ii) = (psi_x_est(ii) - psi_x) ^ 2 * Prob(ii); % Eq. (36)
            error_y(ii) = (psi_y_est(ii) - psi_y) ^ 2 * Prob(ii); % Eq. (37)
        end
        mse_x_theo_cache(n_MonteCarlo) = sum(error_x);
        mse_y_theo_cache(n_MonteCarlo) = sum(error_y);
    end
    mse_x_theo(n_SNR) = mean(mse_x_theo_cache);
    mse_y_theo(n_SNR) = mean(mse_y_theo_cache);
    toc
end
figure;
semilogy(mse_x_theo, 'o');
hold on
semilogy(mse_y_theo, 'x');
mse_x_2_theo = mse_x_theo;
save mse_x_2_theo mse_x_2_theo
mse_y_2_theo = mse_x_theo;
save mse_y_2_theo mse_y_2_theo
clc;
clearvars;
close all;
M_set = [10 10 10 11 11 11 12 12 12];
Thick_set = [7 8 9 7 8 9 7 8 9];
Spacing_set = [2 2 2 3 3 3 4 4 4];
L_set = [5 6 7 5 6 7 5 6 7];
counter = 0; % Increases from 1 to 81
c = 3 * 10 ^ 8; % Speed of light
f0 = 60 * 10 ^ 9; % Frequency: 60GHz
lambda = c / f0; % Radio wavelength
d_BS = lambda / 2; % Antenna spacing
N_x = 2; % Number of antennas in the x-direction
N_y = 2; % Number of antennas in the y-direction
N = N_x * N_y; % Number of antennas
DFT_x_real = exp(-1i / N_x * 2 * pi * (0 : N_x - 1).' * (0 : N_x - 1));
DFT_y_real = exp(-1i / N_y * 2 * pi * (0 : N_y - 1).' * (0 : N_y - 1));
DFT = kron(DFT_y_real, DFT_x_real); % Generation of the 2D-DFT matrix
DFT_vec = vec(DFT); % Vectorization
NMSE_cache = zeros(9, 9);
for ii_counter = 1:9
    for jj_counter = 1:9
        tic
        rng(1);
        d_x = 2 * lambda / Spacing_set(jj_counter); % Meta-atom spacing in the x-direction
        d_y = 2 * lambda / Spacing_set(jj_counter); % Meta-atom spacing in the y-direction
        L = L_set(jj_counter); % Number of metasurface layers
        M_x = M_set(ii_counter);
        M_y = M_set(ii_counter);
        M = M_x * M_y; % Number of meta-atoms per layer
        d_layer = Thick_set(ii_counter) * lambda / L; % Spacing between two adjacent layers
        W_T = zeros(M, M); % Propagation matrix between two adjacent metasurface layers
        W_T_L = zeros(N, M); % Propagation coefficient matrix from the L-th metasurface layer to the antenna array
        W_T_0 = zeros(M, N); % Propagation coefficient matrix from the input layer to the 1-st metasurface layer
        %% Calculate the inter-layer propagation matrix
        for mm1 = 1 : M
            m1_y = ceil(mm1 / M_x);
            m1_x = mod(mm1 - 1, M_x) + 1;
            for mm2 = 1 : M
                m2_y = ceil(mm2 / M_x);
                m2_x = mod(mm2 - 1, M_x) + 1;
                d_mm1_mm2  = sqrt((m1_x - m2_x) ^ 2 + (m1_y - m2_y) ^ 2) * d_x; % Spacing between the mm1-th meta-atom and the mm2-th meta-atom on the same layer
                d_temp = sqrt(d_layer ^ 2 + d_mm1_mm2 ^ 2);
                W_T(mm2, mm1) = d_x * d_y * ...
                    (d_layer / d_temp / d_temp * (1 / 2 / pi / d_temp - 1i / lambda)) * ...
                    exp(1i * 2 * pi * d_temp / lambda); % Rayleigh-Sommerfeld diffraction equation
            end
            for nn = 1 : N
                n2_y = ceil(nn / N_x);
                n2_x = mod(nn - 1, N_x) + 1;
                d_temp = sqrt(d_layer ^ 2 + ...
                    ((m1_x - (1 + M_x) / 2) * d_x - (n2_x - (1 + N_x) / 2) * d_BS) ^ 2 + ...
                    ((m1_y - (1 + M_y) / 2) * d_y - (n2_y - (1 + N_y) / 2) * d_BS) ^ 2); % Distance from the mm1-th meta-atom on the L-th metasurface layer (nn-th meta-atom on the input layer) to the nn-th antenna (mm1-th meta-atom on the 1-st metasurface layer)
                W_T_L(nn, mm1) = d_x * d_y * ...
                    (d_layer / d_temp / d_temp * (1 / 2 / pi / d_temp - 1i / lambda)) * ...
                    exp(1i * 2 * pi * d_temp / lambda); % Rayleigh-Sommerfeld diffraction equation
                W_T_0(mm1, nn) = d_x * d_y * ...
                    (d_layer / d_temp / d_temp * (1 / 2 / pi / d_temp - 1i / lambda)) * ...
                    exp(1i * 2 * pi * d_temp / lambda); % Rayleigh-Sommerfeld diffraction equation
            end
        end
        %% Initialization
        transmission_coefficient = randn(M, L) + 1i * randn(M, L);
        transmission_coefficient = transmission_coefficient ./ abs(transmission_coefficient); % Initialization of the transmission coefficients of meta-atoms in the SIM
        phase_shift = angle(transmission_coefficient); % Corresponding phase shifts
        W_SIM = diag(transmission_coefficient(:, 1));
        for l = 1 : (L - 1)
            W_SIM = diag(transmission_coefficient(:, l + 1)) * W_T * W_SIM;
        end
        Response = W_T_L * W_SIM * W_T_0; % Calculate the response of the SIM
        Response_vec = vec(Response); % Vectorization
        Scaling_factor = (Response_vec' * Response_vec) \ (Response_vec' * DFT_vec); % Calculate the scaling factor
        nmse_temp = norm(Response_vec * Scaling_factor - DFT_vec) ^ 2 / norm(DFT_vec) ^ 2; % Calculate the NMSE
        learning_rate = 1;
        %% Utilize the gradient descent algorithm to configure the SIM
        for jj = 1 : 100
            learning_rate = learning_rate * 0.8; % Decay rate: 0.8
            for ll = 1 : L
                for mm1 = 1 : M
                    X_left = W_T_0;
                    for ll_left = 1 : (ll - 1)
                        X_left = W_T * diag(transmission_coefficient(:, ll_left)) * X_left;
                    end
                    X_right = W_T_L;
                    for ll_right = 1 : (L - ll)
                        X_right = X_right * diag(transmission_coefficient(:, L + 1 - ll_right)) * W_T;
                    end
                    Gradient_phase_shift(mm1, ll) = sum(sum(-2 * imag(conj(Scaling_factor * Response - DFT) .* ...
                        (Scaling_factor * X_right(:, mm1) * X_left(mm1, :) * transmission_coefficient(mm1, ll))))); % Calculate the gradient based on Eq. (19)
                end
            end
            Gradient_phase_shift = Gradient_phase_shift / max(max(Gradient_phase_shift)) * pi;
            phase_shift = phase_shift - learning_rate * Gradient_phase_shift; % Update the phase shifts of meta-atoms in the SIM
            transmission_coefficient = exp(1i * phase_shift); % Update the transmission coefficients of meta-atoms in the SIM
            W_SIM = diag(transmission_coefficient(:, 1));
            for l = 1 : (L - 1)
                W_SIM = diag(transmission_coefficient(:, l + 1)) * W_T * W_SIM;
            end
            Response = W_T_L * W_SIM * W_T_0; % Update the response of the SIM
            Response_vec = vec(Response); % Vectorization
            Scaling_factor = (Response_vec' * Response_vec) \ (Response_vec' * DFT_vec); % Update the scaling factor
            nmse_temp(jj) = norm(Response_vec * Scaling_factor - DFT_vec) ^ 2 / norm(DFT_vec) ^ 2; % Update the NMSE
        end
        NMSE_cache(ii_counter, jj_counter) = nmse_temp(end);
        counter = counter + 1 % Update the counter
        clear Gradient_phase_shift
        toc
    end
end
NMSE = 10 * log10(NMSE_cache); % NMSE in dB
save NMSE NMSE
min(min(NMSE))
figure;
plot(vec(NMSE));
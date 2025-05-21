clc;
clearvars;
close all;
N_x = 2; % Number of antennas in the x-direction
N_y = 2; % Number of antennas in the y-direction
N = N_x * N_y; % Number of antennas
DFT_x_real = exp(-1i / N_x * 2 * pi * (0 : N_x - 1).' * (0 : N_x - 1));
DFT_y_real = exp(-1i / N_y * 2 * pi * (0 : N_y - 1).' * (0 : N_y - 1));
DFT = kron(DFT_y_real, DFT_x_real); % Generation of the 2D-DFT matrix
%% True electrical angles
rng(14)
psi_x = -rand % Elesctrical angle \psi_x in the x direction
psi_y = -rand % Elesctrical angle \psi_y in the y direction
SV_x = exp(1i * pi * (0 : N_x - 1).' * psi_x);
SV_y = exp(1i * pi * (0 : N_y - 1).' * psi_y);
SV = kron(SV_y, SV_x); % Calculate the steering vector
N_snapshot = 64; % Number of snapshots for each discrete bin, i.e., T_x and T_y in (25)
Spectrum = zeros(N_snapshot * N_x, N_snapshot * N_y); % Fine-resolution angle spectrum
p = zeros(N_snapshot, N_snapshot); % Output gain of the coarse-resolution angle spectrum
%% Find the optimal configuration of the input layer
for ii = 1 : N_snapshot
    for jj = 1 : N_snapshot
        text_x = (ii - 1) / N_snapshot;
        text_y = (jj - 1) / N_snapshot;
        offset_x = exp(-1i / N_x * 2 * pi * (0 : N_x - 1).' * text_x);
        offset_y = exp(-1i / N_y * 2 * pi * (0 : N_y - 1).' * text_y);
        offset = kron(offset_y, offset_x); % Transmission coefficients of the input layer
        output = DFT * diag(offset) * SV; % DFT in the digital domain
        p(ii, jj) = max(abs(output));
        spectrum_coarse = reshape(output, N_x, N_y); % Coarse-resolution spectrum in the ii-th slot of the jj-th block
        Spectrum(ii : N_snapshot : N_snapshot * N_x, jj : N_snapshot : N_snapshot * N_y) = spectrum_coarse;
    end
end
Spectrum = Spectrum.';
p_peak = max(max(p)); % The peak of the angle spectrum
[row, column] = find(p == p_peak); % The index corresponding to the peak
offset_x = exp(-1i * 2 * pi / N_x * (0 : N_x - 1).' * (row - 1) / N_snapshot);
offset_y = exp(-1i * 2 * pi / N_y * (0 : N_y - 1).' * (column - 1) / N_snapshot);
offset = kron(offset_y, offset_x); % Transmission coefficients of the input layer
output = DFT * diag(offset) * SV; % DFT in the digital domain
%% Estimated electrical angles
[~, index] = max(abs(output));
y_index = ceil(index / N_x);
x_index = mod(index - 1, N_x) + 1;
psi_x_est = mod(((x_index - 1) / N_x + (row - 1) / N_snapshot /N_x) * 2 + 1, 2) - 1 % Estimated elesctrical angle \psi_x in the x direction
psi_y_est = mod(((y_index - 1) / N_y + (column - 1) / N_snapshot / N_y) * 2 + 1, 2) - 1 % Estimated elesctrical angle \psi_y in the y direction
%% Plot the angle specturm
xx_set = -1 : 2 / (N_snapshot * N_x) : 1 - 2 / (N_snapshot * N_x);
yy_set = -1 : 2 / (N_snapshot * N_y) : 1 - 2 / (N_snapshot * N_y);
figure;
Spectrum_gain = abs(Spectrum) .^ 2 / N ^ 2;
set(gca, 'Position', [0.2 0.2 0.7 0.7]);
imagesc(xx_set, yy_set, Spectrum_gain);
Spectrum_gain_peak = max(max(Spectrum_gain));
[row, column] = find(Spectrum_gain == Spectrum_gain_peak);
hold on
cross = scatter(xx_set(column), yy_set(row), 'x');
cross.SizeData = 100;
cross.LineWidth = 2.5;
cross.MarkerEdgeColor = 'r';
xlabel('\psi_x')
ylabel('\psi_y')
colorbar
colormap parula
set(gca, 'fontsize', 18);
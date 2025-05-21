clc;
clearvars;
close all;
SNR_set = 10 : 1 : 30; % SNR
load("mse_x_4_SIM.mat")
load("mse_y_4_SIM.mat")
load("mse_x_4_theo.mat")
load("mse_y_4_theo.mat")
load("mse_x_8_SIM.mat")
load("mse_y_8_SIM.mat")
load("mse_x_8_theo.mat")
load("mse_y_8_theo.mat")
figure;
semilogy(SNR_set, mse_x_4_SIM, '--', 'linewidth', 1.7, 'markersize', 12, 'color', [250, 128, 114] / 255)
hold on
semilogy(SNR_set, mse_y_4_SIM, '-.', 'linewidth', 1.7, 'markersize', 10, 'color', [0, 191, 255] / 255)
hold on
semilogy(SNR_set, mse_x_4_theo, ':x', 'linewidth', 1.7, 'markersize', 12, 'color', [60, 179, 113] / 255)
hold on
semilogy(SNR_set, mse_y_4_theo, ':o', 'linewidth', 1.7, 'markersize', 10, 'color', [60, 179, 113] / 255)
hold on
semilogy(SNR_set,mse_x_8_SIM, '--', 'linewidth', 1.7, 'markersize', 12, 'color', [250, 128, 114] / 255)
hold on
semilogy(SNR_set, mse_y_8_SIM, '-.', 'linewidth', 1.7, 'markersize', 10, 'color', [0, 191, 255] / 255)
hold on
semilogy(SNR_set, mse_x_8_theo, ':x', 'linewidth', 1.7, 'markersize', 12, 'color', [60, 179, 113] / 255)
hold on
semilogy(SNR_set, mse_y_8_theo, ':o', 'linewidth', 1.7, 'markersize', 10, 'color', [60, 179, 113] / 255)
hold on
%% Annotation
x = [0.49 0.44];
y = [0.65 0.7];
annotation('textarrow', x, y, 'String', 'T_x = T_y = 4', 'linewidth', 1.2, 'fontsize', 14, 'linestyle', ':')
dim = [.41 .7 .03 .1];
annotation('ellipse', dim, 'linewidth', 1.2, 'linestyle', ':')
x = [0.36 0.41];
y = [0.21 0.26];
annotation('textarrow', x, y, 'String', 'T_x = T_y = 8', 'linewidth', 1.2, 'fontsize', 14, 'linestyle', ':')
dim = [.4 .25 .05 .25];
annotation('ellipse', dim, 'linewidth', 1.2, 'linestyle', ':')
grid on
maxvalue = mse_y_4_theo(1) * 1.2;
minvalue = mse_x_8_SIM(end) / 1.2;
axis([10 30 minvalue maxvalue])
xlabel('SNR [dB]')
ylabel('MSE')
legend('\psi_x, Sim.', '\psi_y, Sim.', '\psi_x, The.', '\psi_y, The.', 'location', 'east')
set(gca, 'fontsize', 14)
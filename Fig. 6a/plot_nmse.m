clc;
clearvars;
close all;
load('NMSE_90.mat');
load('NMSE_95.mat');
load('NMSE_99.mat');
figure;
semilogy(NMSE_90, '--', 'linewidth', 2, 'color', [250,128,114]/255);
hold on
semilogy(NMSE_95, '-', 'linewidth', 2, 'color', [0,191,255]/255);
hold on
semilogy(NMSE_99, '-.', 'linewidth', 2, 'color', [60,179,113]/255);
xlabel('Number of iterations');
yLabel = sprintf("\\mathcal{L}");
ylabel(strcat("$", yLabel, "$"), 'Interpreter', 'latex');
grid on
set(gca, 'fontsize', 18);
legend('\zeta = 0.90', '\zeta = 0.95', '\zeta = 0.99', 'location', 'southwest');
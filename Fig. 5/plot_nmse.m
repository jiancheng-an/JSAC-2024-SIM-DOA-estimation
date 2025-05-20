clc;
clearvars;
close all;
load('NMSE_75.mat');
load('NMSE_80.mat');
load('NMSE_90.mat');
figure;
semilogy(NMSE_75, '--', 'linewidth', 2, 'color', [250,128,114]/255);
hold on
semilogy(NMSE_80, '-', 'linewidth', 2, 'color', [0,191,255]/255);
hold on
semilogy(NMSE_90, '-.', 'linewidth', 2, 'color', [60,179,113]/255);
xlabel('Number of iterations');
yLabel = sprintf("\\mathcal{L}");
ylabel(strcat("$", yLabel, "$"), 'Interpreter', 'latex');
grid on
set(gca, 'fontsize', 18);
legend('\zeta = 0.75', '\zeta = 0.80', '\zeta = 0.90', 'location', 'southwest');
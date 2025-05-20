clc;
clearvars;
close all;
load('NMSE_mean.mat');
load('NMSE_lower_bound.mat');
load('NMSE_upper_bound.mat');
zeta_set = 0.7 : 0.0001 : 1;
figure;
semilogy(zeta_set, NMSE_mean, '--', 'linewidth', 2, 'color', [0, 191, 255] / 255);
hold on
pic = fill([zeta_set, fliplr(zeta_set)], [NMSE_lower_bound, fliplr(NMSE_upper_bound)], 'r');
set(pic,'edgealpha', 0, 'facealpha', 0.4);
xlabel('\zeta')
yLabel = sprintf("\\mathcal{L}");
ylabel(strcat("$", yLabel, "$"), 'Interpreter', 'latex');
grid on
set(gca,'fontsize', 18);
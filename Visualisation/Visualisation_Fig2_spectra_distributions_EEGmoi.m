%% LEAP pipelines: Spectra and distributions for EEG metrics of interest

% This script takes the resulting global power and connectivity from each
% of the pipelines and plots the spectra and distributions of power and
% connectivity. 

% The rest of the code creates raincloud plots of the data to show the
% distributions of values for each pipeline. This section uses the
% Raincloud plot MATLAB code from Micah Allan and colleagues:
% Allen M, Poggiali D, Whitaker K et al. Raincloud plots: a multi-platform 
% tool for robust data visualization [version 2; peer review: 2 approved]. 
% Wellcome Open Res 2021, 4:63 (https://doi.org/10.12688/wellcomeopenres.15191.2)


% Created by Rianne Haartsen, PhD.; 12-2023 
% Birkbeck University of London

%% Add paths and read in data

addpath('xxx')
cd xxx/Figs_v3_Mar23/data_csv 

% Read in data
Pow_manual = readmatrix('Manual_pow_spectra_alltrls.csv');
Pow_MADE = readmatrix('MADE_pow_spectra_alltrls.csv');
Pow_BOND = readmatrix('BOND_pow_spectra_alltrls.csv');
Pow_HAPPE1 = readmatrix('HAPPE_pow_spectra_alltrls.csv');

FC_manual = readmatrix('Manual_FC_spectra_alltrls.csv');
FC_MADE = readmatrix('MADE_FC_spectra_alltrls.csv');
FC_BOND = readmatrix('BOND_FC_spectra_alltrls.csv');
FC_HAPPE1 = readmatrix('HAPPE_FC_spectra_alltrls.csv');


%% Spectra for power and connectivity
spectraFig = figure;
% power
subplot(2,1,1)
freqs = 1:1:32;

mnPow = mean(Pow_manual,1);
sdPow = std(Pow_manual,0,1);
curve1 = mnPow + sdPow;
curve2 = mnPow - sdPow;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.1, 0.4, 0.8],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnPow, 'Color',[0.1, 0.4, 0.8], 'LineWidth', 2);
clear mnPow sdPow curve1 curve2 inBetween freqs2

mnPow = mean(Pow_MADE,1);
sdPow = std(Pow_MADE,0,1);
curve1 = mnPow + sdPow;
curve2 = mnPow - sdPow;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.9, 0.5, 0],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnPow, 'Color',[0.9, 0.5, 0], 'LineWidth', 2);
clear mnPow sdPow curve1 curve2 inBetween freqs2

mnPow = mean(Pow_BOND,1);
sdPow = std(Pow_BOND,0,1);
curve1 = mnPow + sdPow;
curve2 = mnPow - sdPow;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.8, 0, 0.6],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnPow, 'Color',[0.8, 0, 0.6], 'LineWidth', 2);
clear mnPow sdPow curve1 curve2 inBetween freqs2

mnPow = mean(Pow_HAPPE1,1);
sdPow = std(Pow_HAPPE1,0,1);
curve1 = mnPow + sdPow;
curve2 = mnPow - sdPow;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.4, 0.6, 0.1],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnPow, 'Color',[0.4, 0.6, 0.1], 'LineWidth', 2);
clear mnPow sdPow curve1 curve2 inBetween freqs2

xlim([1, 32])
xticks([1 2 3 4 6 7 12 13 30])
xline(1.5); xline(3.5); xline(6.5); xline(12.5); xline(30.5)
xlabel('Frequency (Hz)'); ylabel ('Power (log)')

% FC
subplot(2,1,2)
freqs = 1:1:32;

mnFC = mean(FC_manual,1);
sdFC = std(FC_manual,0,1);
curve1 = mnFC + sdFC;
curve2 = mnFC - sdFC;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.1, 0.4, 0.8],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnFC, 'Color',[0.1, 0.4, 0.8], 'LineWidth', 2);
clear mnFC sdFC curve1 curve2 inBetween freqs2

mnFC = mean(FC_MADE,1);
sdFC = std(FC_MADE,0,1);
curve1 = mnFC + sdFC;
curve2 = mnFC - sdFC;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.9, 0.5, 0],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnFC, 'Color',[0.9, 0.5, 0], 'LineWidth', 2);
clear mnFC sdFC curve1 curve2 inBetween freqs2

mnFC = mean(FC_BOND,1);
sdFC = std(FC_BOND,0,1);
curve1 = mnFC + sdFC;
curve2 = mnFC - sdFC;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.8, 0, 0.6],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnFC, 'Color',[0.8, 0, 0.6], 'LineWidth', 2);
clear mnFC sdFC curve1 curve2 inBetween freqs2

mnFC = mean(FC_HAPPE1,1);
sdFC = std(FC_HAPPE1,0,1);
curve1 = mnFC + sdFC;
curve2 = mnFC - sdFC;
freqs2 = [freqs, fliplr(freqs)];
inBetween = [curve1, fliplr(curve2)];
h = fill(freqs2, inBetween, [0.4, 0.6, 0.1],'EdgeColor','none');
set(h, 'facealpha', 0.2)
hold on;
plot(freqs, mnFC, 'Color',[0.4, 0.6, 0.1], 'LineWidth', 2);
clear mnFC sdFC curve1 curve2 inBetween freqs2

xlim([1, 32])
xticks([1 2 3 4 6 7 12 13 30])
xline(1.5); xline(3.5); xline(6.5); xline(12.5); xline(30.5)
ylim([0 .06])
xlabel('Frequency (Hz)'); ylabel ('Connectivity (dbWPLI)')

%% Distributions frequency bands

addpath /Users/riannehaartsen/Documents/MATLAB/RainCloudPlots-master/tutorial_matlab

figure
% power
subplot(2,4,1) % delta
data_cur{1} = mean(Pow_manual(:,2:3),2);
data_cur{2} = mean(Pow_MADE(:,2:3),2);
data_cur{3} = mean(Pow_BOND(:,2:3),2);
data_cur{4} = mean(Pow_HAPPE1(:,2:3),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-.6 .8]);
box off 
view([-90 90]);
clear data_cur
title('Delta')
xlabel('Global power (log)')

subplot(2,4,2) % theta
data_cur{1} = mean(Pow_manual(:,4:6),2);
data_cur{2} = mean(Pow_MADE(:,4:6),2);
data_cur{3} = mean(Pow_BOND(:,4:6),2);
data_cur{4} = mean(Pow_HAPPE1(:,4:6),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-.6 .8]);
box off 
view([-90 90]);
clear data_cur
title('Theta')

subplot(2,4,3) % alpha
data_cur{1} = mean(Pow_manual(:,7:12),2);
data_cur{2} = mean(Pow_MADE(:,7:12),2);
data_cur{3} = mean(Pow_BOND(:,7:12),2);
data_cur{4} = mean(Pow_HAPPE1(:,7:12),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-.6 .8]);
box off 
view([-90 90]);
clear data_cur
title('Alpha')

subplot(2,4,4) % beta
data_cur{1} = mean(Pow_manual(:,13:30),2);
data_cur{2} = mean(Pow_MADE(:,13:30),2);
data_cur{3} = mean(Pow_BOND(:,13:30),2);
data_cur{4} = mean(Pow_HAPPE1(:,13:30),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-.6 .8]);
box off 
view([-90 90]);
clear data_cur
title('Beta')


% connectivity
subplot(2,4,5) % delta
data_cur{1} = mean(FC_manual(:,2:3),2);
data_cur{2} = mean(FC_MADE(:,2:3),2);
data_cur{3} = mean(FC_BOND(:,2:3),2);
data_cur{4} = mean(FC_HAPPE1(:,2:3),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-100 120]);
box off 
view([-90 90]);
clear data_cur
title('Delta')
xlabel('Global connectivity (dbWPLI)')

subplot(2,4,6) % theta
data_cur{1} = mean(FC_manual(:,4:6),2);
data_cur{2} = mean(FC_MADE(:,4:6),2);
data_cur{3} = mean(FC_BOND(:,4:6),2);
data_cur{4} = mean(FC_HAPPE1(:,4:6),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-80 120]);
box off 
view([-90 90]);
clear data_cur
title('Theta')

subplot(2,4,7) % alpha
data_cur{1} = mean(FC_manual(:,7:12),2);
data_cur{2} = mean(FC_MADE(:,7:12),2);
data_cur{3} = mean(FC_BOND(:,7:12),2);
data_cur{4} = mean(FC_HAPPE1(:,7:12),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-50 80]);
box off 
view([-90 90]);
clear data_cur
title('Alpha')

subplot(2,4,8) % beta
data_cur{1} = mean(FC_manual(:,13:30),2);
data_cur{2} = mean(FC_MADE(:,13:30),2);
data_cur{3} = mean(FC_BOND(:,13:30),2);
data_cur{4} = mean(FC_HAPPE1(:,13:30),2);
h1 = raincloud_plot(data_cur{1}, 'box_on', 1, 'color', [0.1, 0.4, 0.8], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0,'line_width',1);
h2 = raincloud_plot(data_cur{2}, 'box_on', 1, 'color', [0.9, 0.5, 0], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
h3 = raincloud_plot(data_cur{3}, 'box_on', 1, 'color', [0.8, 0, 0.6], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'line_width',1);
h4 = raincloud_plot(data_cur{4}, 'box_on', 1, 'color', [0.4, 0.6, 0.1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
set(gca, 'YLim', [-110 160]);
box off 
view([-90 90]);
clear data_cur
title('Beta')

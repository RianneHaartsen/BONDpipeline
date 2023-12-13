%% LEAP pipelines: Between pipeline comparisons

% For 4 pipelines: Manual, MADE original, BOND-MADE, HAPPE

% 1) Inclusion rates: number of trials for all, soc, toy, and differences
% between conditions
% 2) Power values: for theta and alpha frequencies
% 3) Connectivity values: for theta and alpha frequencies
% 4) Spectrum plots: power and fc across all trls for inspection

% This script also calls to the ICC script written by Arash Salarian:
% Arash Salarian (2023). Intraclass Correlation Coefficient (ICC) 
% (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc), 
% MATLAB Central File Exchange. 

% export to csv files for python


% Created by Rianne Haartsen, PhD.; 06-2023 
% Birkbeck University of London

%% add folder with extra scripts
addpath('xxx/ICC') 

%% Inclusion rates

% Manual
load xxx/Manual_MeasuresOfInterest_report.mat
Incl_Manual = Manual_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_Manual)
    if isempty(Incl_Manual.Ntrls_tot{ii}) || isnan(Incl_Manual.Ntrls_tot{ii})
        Incl_Manual.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_Manual.Ntrls_soc{ii}) || isnan(Incl_Manual.Ntrls_soc{ii})
        Incl_Manual.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_Manual.Ntrls_toy{ii}) || isnan(Incl_Manual.Ntrls_toy{ii})
        Incl_Manual.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
writetable(Incl_Manual,'Manual_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_Manual Manual_MOI_table

% MADE
load xxx/MADE_MeasuresOfInterest_report.mat
Incl_MADE = MADE_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_MADE)
    if isempty(Incl_MADE.Ntrls_tot{ii}) || isnan(Incl_MADE.Ntrls_tot{ii})
        Incl_MADE.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_MADE.Ntrls_soc{ii}) || isnan(Incl_MADE.Ntrls_soc{ii})
        Incl_MADE.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_MADE.Ntrls_toy{ii}) || isnan(Incl_MADE.Ntrls_toy{ii})
        Incl_MADE.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
writetable(Incl_MADE,'MADE_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_MADE MADE_MOI_table

% BOND-MADE
load xxx/BOND_MeasuresOfInterest_report.mat
Incl_BOND = BOND_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_BOND)
    if isempty(Incl_BOND.Ntrls_tot{ii}) || isnan(Incl_BOND.Ntrls_tot{ii})
        Incl_BOND.Ntrls_tot{ii} = 0;
    end
    if isempty(Incl_BOND.Ntrls_soc{ii}) || isnan(Incl_BOND.Ntrls_soc{ii})
        Incl_BOND.Ntrls_soc{ii} = 0;
    end
    if isempty(Incl_BOND.Ntrls_toy{ii}) || isnan(Incl_BOND.Ntrls_toy{ii})
        Incl_BOND.Ntrls_toy{ii} = 0;
    end
end
clear ii
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
writetable(Incl_BOND,'BOND_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_BOND BOND_MOI_table

% HAPPE
load xxx/HAPPE_MeasuresOfInterest_report.mat
Incl_HAPPE = HAPPE_MOI_table(:, [1, 5, 6, 7]);
for ii = 1:height(Incl_HAPPE)
    if isempty(Incl_HAPPE.Neps_tot{ii}) || isnan(Incl_HAPPE.Neps_tot{ii})
        Incl_HAPPE.Neps_tot{ii} = 0;
    end
    if isempty(Incl_HAPPE.Neps_soc{ii}) || isnan(Incl_HAPPE.Neps_soc{ii})
        Incl_HAPPE.Neps_soc{ii} = 0;
    end
    if isempty(Incl_HAPPE.Neps_toy{ii}) || isnan(Incl_HAPPE.Neps_toy{ii})
        Incl_HAPPE.Neps_toy{ii} = 0;
    end
end
clear ii
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
writetable(Incl_HAPPE,'HAPPE_InclRates.csv','Delimiter',',','QuoteStrings',1)
clear Incl_HAPPE HAPPE_MOI_table


%% Spectra
load xxx/Manual_MeasuresOfInterest_report.mat
load xxx/MADE_MeasuresOfInterest_report.mat
load xxx/BOND_MeasuresOfInterest_report.mat
load xxx/HAPPE_MeasuresOfInterest_report.mat

% Power
Manual_Pow_spectra_all = cell2mat(table2array(Manual_MOI_table(:, 10)));
MADE_Pow_spectra_all = cell2mat(table2array(MADE_MOI_table(:, 10)));
BOND_Pow_spectra_all = cell2mat(table2array(BOND_MOI_table(:, 10)));
HAPPE_Pow_spectra_all = cell2mat(table2array(HAPPE_MOI_table(:, 10)));

% set those below threshold to NaNs
Ntrls = [cell2mat(Manual_MOI_table.Ntrls_soc) cell2mat(Manual_MOI_table.Ntrls_toy)...
    cell2mat(MADE_MOI_table.Ntrls_soc) cell2mat(MADE_MOI_table.Ntrls_toy) ...
    cell2mat(BOND_MOI_table.Ntrls_soc) cell2mat(BOND_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPE_MOI_table.Neps_soc) cell2mat(HAPPE_MOI_table.Neps_toy)];
Ind_incl = find(sum(Ntrls >= 20, 2) ==8);
N115_Manual_Pow_spectra = Manual_Pow_spectra_all(Ind_incl,:);
N115_MADE_Pow_spectra = MADE_Pow_spectra_all(Ind_incl,:);
N115_BOND_Pow_spectra = BOND_Pow_spectra_all(Ind_incl,:);
N115_HAPPE_Pow_spectra = HAPPE_Pow_spectra_all(Ind_incl,:);
Ind_incl_pow = Ind_incl;
clear Ind_incl

% FC
Manual_FC_spectra_all = cell2mat(table2array(Manual_MOI_table(:, 14)));
MADE_FC_spectra_all = cell2mat(table2array(MADE_MOI_table(:, 14)));
BOND_FC_spectra_all = cell2mat(table2array(BOND_MOI_table(:, 14)));
HAPPE_FC_spectra_all = cell2mat(table2array(HAPPE_MOI_table(:, 14)));

% set those below threshold to NaNs
Ntrls_FC = [cell2mat(Manual_MOI_table.Ntrls_soc) cell2mat(Manual_MOI_table.Ntrls_toy)...
    cell2mat(MADE_MOI_table.Ntrls_soc) cell2mat(MADE_MOI_table.Ntrls_toy) ...
    cell2mat(BOND_MOI_table.Ntrls_soc) cell2mat(BOND_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPE_MOI_table.Neps_soc) cell2mat(HAPPE_MOI_table.Neps_toy)];
Ind_incl = find(sum(Ntrls_FC >= 90, 2) ==8);
N105_Manual_FC_spectra = Manual_FC_spectra_all(Ind_incl,:);
N105_MADE_FC_spectra = MADE_FC_spectra_all(Ind_incl,:);
N105_BOND_FC_spectra = BOND_FC_spectra_all(Ind_incl,:);
N105_HAPPE_FC_spectra = HAPPE_FC_spectra_all(Ind_incl,:);
Ind_incl_FC = Ind_incl;

% save data
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
% power
writematrix(N115_Manual_Pow_spectra, 'Manual_pow_spectra_alltrls.csv')
writematrix(N115_MADE_Pow_spectra, 'MADE_pow_spectra_alltrls.csv')
writematrix(N115_BOND_Pow_spectra, 'BOND_pow_spectra_alltrls.csv')
writematrix(N115_HAPPE_Pow_spectra, 'HAPPE_pow_spectra_alltrls.csv')
% FC
writematrix(N105_Manual_FC_spectra, 'Manual_FC_spectra_alltrls.csv')
writematrix(N105_MADE_FC_spectra, 'MADE_FC_spectra_alltrls.csv')
writematrix(N105_BOND_FC_spectra, 'BOND_FC_spectra_alltrls.csv')
writematrix(N105_HAPPE_FC_spectra, 'HAPPE_FC_spectra_alltrls.csv')
% other
save('Incl_indices.mat','Ind_incl_pow', 'Ind_incl_FC')

%% ICC between frequency bands and pipelines (across all trials)
Freqs = 1:1:32;
Delta_ind = [find(Freqs == 2),find(Freqs == 3)];
Theta_ind = [find(Freqs == 4),find(Freqs == 6)];
Alpha_ind = [find(Freqs == 7),find(Freqs == 12)];
Beta_ind = [find(Freqs == 13),find(Freqs == 30)];
 
% For power
% Manual - delta, theta, alpha, beta
Manual_pow = [mean(N115_Manual_Pow_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_Manual_Pow_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_Manual_Pow_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_Manual_Pow_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)];
% MADE
MADE_pow = [mean(N115_MADE_Pow_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)];
% BOND
BOND_pow = [mean(N115_BOND_Pow_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)];    
% HAPPE
HAPPE_pow = [mean(N115_HAPPE_Pow_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
% Combined
Power_freqbands_all = [Manual_pow MADE_pow BOND_pow HAPPE_pow];
% run ICC for all combinations
r_vals = zeros(size(Power_freqbands_all,2));
LB_vals = zeros(size(Power_freqbands_all,2));
UB_vals = zeros(size(Power_freqbands_all,2));
p_vals = zeros(size(Power_freqbands_all,2));
for rr = 1:size(Power_freqbands_all,2)
    for cc = 1:size(Power_freqbands_all,2)
        Mat = [Power_freqbands_all(:,rr) Power_freqbands_all(:,cc)];
        [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
        r_vals(rr, cc) = r;
        LB_vals(rr, cc) = LB;
        UB_vals(rr, cc) = UB;
        p_vals(rr,cc) = p;
        clear Mar r p
    end
end

r_vals(r_vals < 0) = 0;

% save data
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
% power
writematrix(r_vals, 'Pow_ICCs_alltrls_rvals_masked.csv')
writematrix(LB_vals, 'Pow_ICCs_alltrls_LBvals.csv')
writematrix(UB_vals, 'Pow_ICCs_alltrls_UBvals.csv')
writematrix(p_vals, 'Pow_ICCs_alltrls_pvals.csv')

clear r_vals LB_vals UB_vals p_vals

% For connectivity values
% Manual - delta, theta, alpha, beta
Manual_FC = [mean(N105_Manual_FC_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_Manual_FC_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_Manual_FC_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_Manual_FC_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)];
% MADE
MADE_FC = [mean(N105_MADE_FC_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_MADE_FC_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_MADE_FC_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_MADE_FC_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)];
% BOND
BOND_FC = [mean(N105_BOND_FC_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_BOND_FC_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_BOND_FC_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_BOND_FC_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)];    
% HAPPE
HAPPE_FC = [mean(N105_HAPPE_FC_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
% Combined
FC_freqbands_all = [Manual_FC MADE_FC BOND_FC HAPPE_FC];
% run ICC for all combinations
r_vals = zeros(size(FC_freqbands_all,2));
LB_vals = zeros(size(FC_freqbands_all,2));
UB_vals = zeros(size(FC_freqbands_all,2));
p_vals = zeros(size(FC_freqbands_all,2));
for rr = 1:size(FC_freqbands_all,2)
    for cc = 1:size(FC_freqbands_all,2)
        Mat = [FC_freqbands_all(:,rr) FC_freqbands_all(:,cc)];
        [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
        r_vals(rr, cc) = r;
        LB_vals(rr, cc) = LB;
        UB_vals(rr, cc) = UB;
        p_vals(rr,cc) = p;
        clear Mar r p
    end
end

r_vals(r_vals < 0) = 0;

% save data
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
% power
writematrix(r_vals, 'FC_ICCs_alltrls_rvals_masked.csv')
writematrix(LB_vals, 'FC_ICCs_alltrls_LBvals.csv')
writematrix(UB_vals, 'FC_ICCs_alltrls_UBvals.csv')
writematrix(p_vals, 'FC_ICCs_alltrls_pvals.csv')

clear r_vals LB_vals UB_vals p_vals



%% Condition differences


Freqs = 1:1:32;
Delta_ind = [find(Freqs == 2),find(Freqs == 3)];
Theta_ind = [find(Freqs == 4),find(Freqs == 6)];
Alpha_ind = [find(Freqs == 7),find(Freqs == 12)];
Beta_ind = [find(Freqs == 13),find(Freqs == 30)];

% Power
load /Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Pipeline_comparisons_subset/PreprocData_Roche1.0manual/Roche1_MeasuresOfInterest_report.mat
load /Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Pipeline_comparisons_subset/Preproc_MADE2_v2_Oct/BOND-MADE_MeasuresOfInterest_report.mat
BOND_MOI_table = MADE_MOI_table; clear MADE_MOI_table
load /Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Pipeline_comparisons_subset/Preproc_MADE_v1_June/MADE_MeasuresOfInterest_report_v1_June22.mat
load /Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Pipeline_comparisons_subset/Preproc_HAPPE_v1_Dec/HAPPE_MeasuresOfInterest_report.mat

% Power
Roche1_Pow_spectra_s = cell2mat(table2array(Roche1_MOI_table(:, 11)));
MADE_Pow_spectra_s = cell2mat(table2array(MADE_MOI_table(:, 11)));
BOND_Pow_spectra_s = cell2mat(table2array(BOND_MOI_table(:, 11)));
HAPPE_Pow_spectra_s = cell2mat(table2array(HAPPE_MOI_table(:, 11)));

Roche1_Pow_spectra_t = cell2mat(table2array(Roche1_MOI_table(:, 12)));
MADE_Pow_spectra_t = cell2mat(table2array(MADE_MOI_table(:, 12)));
BOND_Pow_spectra_t = cell2mat(table2array(BOND_MOI_table(:, 12)));
HAPPE_Pow_spectra_t = cell2mat(table2array(HAPPE_MOI_table(:, 12)));

% set those below threshold to NaNs
Ntrls = [cell2mat(Roche1_MOI_table.Ntrls_soc) cell2mat(Roche1_MOI_table.Ntrls_toy)...
    cell2mat(MADE_MOI_table.Ntrls_soc) cell2mat(MADE_MOI_table.Ntrls_toy) ...
    cell2mat(BOND_MOI_table.Ntrls_soc) cell2mat(BOND_MOI_table.Ntrls_toy) ...
    cell2mat(HAPPE_MOI_table.Neps_soc) cell2mat(HAPPE_MOI_table.Neps_toy)];
Ind_incl = find(sum(Ntrls >= 20, 2) ==8);

N115_Roche1_Pow_spectra_s = Roche1_Pow_spectra_s(Ind_incl,:);
N115_MADE_Pow_spectra_s = MADE_Pow_spectra_s(Ind_incl,:);
N115_BOND_Pow_spectra_s = BOND_Pow_spectra_s(Ind_incl,:);
N115_HAPPE_Pow_spectra_s = HAPPE_Pow_spectra_s(Ind_incl,:);

N115_Roche1_Pow_spectra_t = Roche1_Pow_spectra_t(Ind_incl,:);
N115_MADE_Pow_spectra_t = MADE_Pow_spectra_t(Ind_incl,:);
N115_BOND_Pow_spectra_t = BOND_Pow_spectra_t(Ind_incl,:);
N115_HAPPE_Pow_spectra_t = HAPPE_Pow_spectra_t(Ind_incl,:);

Ind_incl_pow = Ind_incl;
clear Ind_incl


% For power differences
% Roche1 - delta, theta, alpha, beta
Roche1_pow_s = [mean(N115_Roche1_Pow_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_Roche1_Pow_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_Roche1_Pow_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_Roche1_Pow_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)];
Roche1_pow_t = [mean(N115_Roche1_Pow_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_Roche1_Pow_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_Roche1_Pow_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_Roche1_Pow_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)];
Roche1_pow_diff = Roche1_pow_s - Roche1_pow_t;

% MADE
MADE_pow_s = [mean(N115_MADE_Pow_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)];
MADE_pow_t = [mean(N115_MADE_Pow_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_MADE_Pow_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)];
MADE_pow_diff = MADE_pow_s - MADE_pow_t;

% BOND
BOND_pow_s = [mean(N115_BOND_Pow_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)];    
BOND_pow_t = [mean(N115_BOND_Pow_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_BOND_Pow_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
BOND_pow_diff = BOND_pow_s - BOND_pow_t;

% HAPPE
HAPPE_pow_s = [mean(N115_HAPPE_Pow_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
HAPPE_pow_t = [mean(N115_HAPPE_Pow_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N115_HAPPE_Pow_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
HAPPE_pow_diff = HAPPE_pow_s - HAPPE_pow_t;

% Combined
Power_freqbands_diff = [Roche1_pow_diff MADE_pow_diff BOND_pow_diff HAPPE_pow_diff];
% run ICC for all combinations
r_vals = zeros(size(Power_freqbands_diff,2));
LB_vals = zeros(size(Power_freqbands_diff,2));
UB_vals = zeros(size(Power_freqbands_diff,2));
p_vals = zeros(size(Power_freqbands_diff,2));
for rr = 1:size(Power_freqbands_diff,2)
    for cc = 1:size(Power_freqbands_diff,2)
        Mat = [Power_freqbands_diff(:,rr) Power_freqbands_diff(:,cc)];
        [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
        r_vals(rr, cc) = r;
        LB_vals(rr, cc) = LB;
        UB_vals(rr, cc) = UB;
        p_vals(rr,cc) = p;
        clear Mar r p
    end
end

r_vals(r_vals < 0) = 0;

% save data
cd /Users/riannehaartsen/Documents/02a_LEAP_EEG/LEAP_SocNSocVids/01_EEGdata/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
% power
writematrix(r_vals, 'Pow_ICCs_diffs_rvals_masked.csv')
writematrix(LB_vals, 'Pow_ICCs_diffs_LBvals.csv')
writematrix(UB_vals, 'Pow_ICCs_diffs_UBvals.csv')
writematrix(p_vals, 'Pow_ICCs_diffs_pvals.csv')

save('Power_diff_data.mat','Manual_pow_s','Manual_pow_t','Manual_pow_diff',...
    'MADE_pow_s','MADE_pow_t','MADE_pow_diff',...
    'BOND_pow_s','BOND_pow_t','BOND_pow_diff',...
    'HAPPE_pow_s','HAPPE_pow_t','HAPPE_pow_diff')

clear r_vals LB_vals UB_vals p_vals


% For connectivity
Roche1_FC_spectra_s = cell2mat(table2array(Roche1_MOI_table(:, 15)));
MADE_FC_spectra_s = cell2mat(table2array(MADE_MOI_table(:, 15)));
BOND_FC_spectra_s = cell2mat(table2array(BOND_MOI_table(:, 15)));
HAPPE_FC_spectra_s = cell2mat(table2array(HAPPE_MOI_table(:, 15)));

Roche1_FC_spectra_t = cell2mat(table2array(Roche1_MOI_table(:, 16)));
MADE_FC_spectra_t = cell2mat(table2array(MADE_MOI_table(:, 16)));
BOND_FC_spectra_t = cell2mat(table2array(BOND_MOI_table(:, 16)));
HAPPE_FC_spectra_t = cell2mat(table2array(HAPPE_MOI_table(:, 16)));

% set those below threshold to NaNs
Ntrls_FC = Ntrls;
Ind_incl = find(sum(Ntrls_FC >= 90, 2) ==8);
N105_Roche1_FC_spectra_s = Roche1_FC_spectra_s(Ind_incl,:);
N105_MADE_FC_spectra_s = MADE_FC_spectra_s(Ind_incl,:);
N105_BOND_FC_spectra_s = BOND_FC_spectra_s(Ind_incl,:);
N105_HAPPE_FC_spectra_s = HAPPE_FC_spectra_s(Ind_incl,:);

N105_Roche1_FC_spectra_t = Roche1_FC_spectra_t(Ind_incl,:);
N105_MADE_FC_spectra_t = MADE_FC_spectra_t(Ind_incl,:);
N105_BOND_FC_spectra_t = BOND_FC_spectra_t(Ind_incl,:);
N105_HAPPE_FC_spectra_t = HAPPE_FC_spectra_t(Ind_incl,:);

Ind_incl_FC = Ind_incl;

% Roche1
Roche1_FC_s = [mean(N105_Roche1_FC_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_Roche1_FC_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_Roche1_FC_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_Roche1_FC_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)];
Roche1_FC_t = [mean(N105_Roche1_FC_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_Roche1_FC_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_Roche1_FC_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_Roche1_FC_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)];
Roche1_FC_diff = Roche1_FC_s - Roche1_FC_t;
% MADE
MADE_FC_s = [mean(N105_MADE_FC_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_MADE_FC_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_MADE_FC_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_MADE_FC_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)];
MADE_FC_t = [mean(N105_MADE_FC_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_MADE_FC_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_MADE_FC_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_MADE_FC_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)];
MADE_FC_diff = MADE_FC_s - MADE_FC_t;
% BOND
BOND_FC_s = [mean(N105_BOND_FC_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_BOND_FC_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_BOND_FC_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_BOND_FC_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)];   
BOND_FC_t = [mean(N105_BOND_FC_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_BOND_FC_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_BOND_FC_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_BOND_FC_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)];  
BOND_FC_diff = BOND_FC_s - BOND_FC_t;
% HAPPE
HAPPE_FC_s = [mean(N105_HAPPE_FC_spectra_s(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra_s(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra_s(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra_s(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
HAPPE_FC_t = [mean(N105_HAPPE_FC_spectra_t(:,[Delta_ind(1):Delta_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra_t(:,[Theta_ind(1):Theta_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra_t(:,[Alpha_ind(1):Alpha_ind(2)]),2)...
    mean(N105_HAPPE_FC_spectra_t(:,[Beta_ind(1):Beta_ind(2)]),2)]; 
HAPPE_FC_diff = HAPPE_FC_s - HAPPE_FC_t;
% Combined
FC_freqbands_diff = [Roche1_FC_diff MADE_FC_diff BOND_FC_diff HAPPE_FC_diff];

% run ICC for all combinations
r_vals = zeros(size(FC_freqbands_diff,2));
LB_vals = zeros(size(FC_freqbands_diff,2));
UB_vals = zeros(size(FC_freqbands_diff,2));
p_vals = zeros(size(FC_freqbands_diff,2));
for rr = 1:size(FC_freqbands_diff,2)
    for cc = 1:size(FC_freqbands_diff,2)
        Mat = [FC_freqbands_diff(:,rr) FC_freqbands_diff(:,cc)];
        [r, LB, UB, ~, ~, ~, p] = ICC(Mat, 'C-k', .05, 0);
        r_vals(rr, cc) = r;
        LB_vals(rr, cc) = LB;
        UB_vals(rr, cc) = UB;
        p_vals(rr,cc) = p;
        clear Mar r p
    end
end

r_vals(r_vals < 0) = 0;

% save data
cd xxx/Pipeline_comparisons_subset/Figs_v3_Mar23/data_csv
% power
writematrix(r_vals, 'FC_ICCs_diff_rvals_masked.csv')
writematrix(LB_vals, 'FC_ICCs_diff_LBvals.csv')
writematrix(UB_vals, 'FC_ICCs_diff_UBvals.csv')
writematrix(p_vals, 'FC_ICCs_diff_pvals.csv')

save('FC_diff_data.mat','Manual_FC_s','Manual_FC_t','Manual_FC_diff',...
    'MADE_FC_s','MADE_FC_t','MADE_FC_diff',...
    'BOND_FC_s','BOND_FC_t','BOND_FC_diff',...
    'HAPPE_FC_s','HAPPE_FC_t','HAPPE_FC_diff')

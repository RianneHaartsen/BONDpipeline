%% LEAP pipelines: Within pipeline comparisons - split-half reliability for power

% This script reads in the datasets with power from the different
% pipelines. Then it takes the power across all trials and creates dataset
% A and B with alternating trials. This method is repeated for the social
% and toy trials. Next, power is averaged across the trials in the datasets
% A and B and the power spectra are saved.

% In the second step, split half reliability is calculated using correlations 
% for the different frequency bands (delta, theta, alpha, beta) and pipelines. 


% Created by Rianne Haartsen, PhD.; 04-2023 
% Birkbeck University of London


%% Power for subset A and subset B
% Roche 1 manual pipeline
cd '/xxx/PreprocData_Manual/' % folder with the manual measures of interest report
load('Manual_MeasuresOfInterest_report.mat')
Power_folder = '/xxx/PreprocData_Manual/A_Power_data';

for ii = 1:height(Manual_MOI_table)
    ID = Manual_MOI_table.ID{ii};
    load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
    Freqs = 1:1:32;
    % across all trials
    if ~isempty(Power_all)
        Ntrls_tot = size(Power_all.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_all = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_all = nan(1,32);
        PowB_all = nan(1,32);
        Nsplit_all = 0;  
    end
    
    % social trials
    if ~isempty(Power_soc)
        Ntrls_tot = size(Power_soc.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_soc = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_soc = nan(1,32);
        PowB_soc = nan(1,32);
        Nsplit_soc = 0;  
    end
    
    % toy trials
    if ~isempty(Power_toy)
        Ntrls_tot = size(Power_toy.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_toy = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_toy = nan(1,32);
        PowB_toy = nan(1,32);
        Nsplit_toy = 0;  
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            Manual_SplitHalf(ii,:) = data_table_newrow;
            cd '/xxx/PreprocData_Manual/'
            save('Manual_SplitHalf_reliability.mat','Manual_SplitHalf')
            clear Manual_SplitHalf
        else
            cd '/xxx/PreprocData_Manual/'
            load Manual_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            Manual_SplitHalf(ii,:) = data_table_newrow;
            cd '/xxx/PreprocData_Manual/'
            save('Manual_SplitHalf_reliability.mat','Manual_SplitHalf')
            clear Manual_SplitHalf
        end

end
clear Roche1_MOI_table


% MADE pipeline
cd xxx/Preproc_MADE
load('MADE_MeasuresOfInterest_report_v1_June22.mat')
Power_folder = 'xxx/Preproc_MADE/A_Power_data';

for ii = 1:height(MADE_MOI_table)
    ID = MADE_MOI_table.ID{ii};
    load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
    Freqs = 1:1:32;
    % across all trials
    if ~isempty(Power_all)
        Ntrls_tot = size(Power_all.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_all = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_all = nan(1,32);
        PowB_all = nan(1,32);
        Nsplit_all = 0;  
    end
    
    % social trials
    if ~isempty(Power_soc)
        Ntrls_tot = size(Power_soc.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_soc = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_soc = nan(1,32);
        PowB_soc = nan(1,32);
        Nsplit_soc = 0;  
    end
    
    % toy trials
    if ~isempty(Power_toy)
        Ntrls_tot = size(Power_toy.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_toy = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_toy = nan(1,32);
        PowB_toy = nan(1,32);
        Nsplit_toy = 0;  
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADE_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_MADE
            save('MADE_SplitHalf_reliability.mat','MADE_SplitHalf')
            clear MADE_SplitHalf
        else
            cd xxx/Preproc_MADE
            load MADE_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            MADE_SplitHalf(ii,:) = data_table_newrow;
            cd /xxx/Preproc_MADE
            save('MADE_SplitHalf_reliability.mat','MADE_SplitHalf')
            clear MADE_SplitHalf
        end

end
clear MADE_MOI_table


% BOND pipeline
cd xxx/Preproc_BOND
load('BOND_MeasuresOfInterest_report.mat')
Power_folder = 'xxx/Preproc_BOND/A_Power_data';

for ii = 1:height(BOND_MOI_table)
    ID = BOND_MOI_table.ID{ii};
    load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
    Freqs = 1:1:32;
    % across all trials
    if ~isempty(Power_all)
        Ntrls_tot = size(Power_all.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_all = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_all = nan(1,32);
        PowB_all = nan(1,32);
        Nsplit_all = 0;  
    end
    
    % social trials
    if ~isempty(Power_soc)
        Ntrls_tot = size(Power_soc.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_soc = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_soc = nan(1,32);
        PowB_soc = nan(1,32);
        Nsplit_soc = 0;  
    end
    
    % toy trials
    if ~isempty(Power_toy)
        Ntrls_tot = size(Power_toy.logpowspctrm,1);
        if mod(Ntrls_tot,2) == 1
            Nmax = Ntrls_tot - 1;
        else
            Nmax = Ntrls_tot;
        end
        Ind_A = 1:2:(Nmax-1);
        Ind_B = 2:2:Nmax;
        PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
        PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
        Nsplit_toy = Nmax/2;
        % clean up 
        clear Ind_A Ind_B Nmax Ntrls_tot
    else
        PowA_toy = nan(1,32);
        PowB_toy = nan(1,32);
        Nsplit_toy = 0;  
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            BOND_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_BOND
            save('BOND_SplitHalf_reliability.mat','BOND_SplitHalf')
            clear BOND_SplitHalf
        else
            cd xxx/Preproc_BOND
            load BOND_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            BOND_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_BOND
            save('BOND_SplitHalf_reliability.mat','BOND_SplitHalf')
            clear BOND_SplitHalf
        end

end
clear BOND_MOI_table


% HAPPE pipeline
cd xxx/Preproc_HAPPE
load('HAPPE_MeasuresOfInterest_report.mat')
Power_folder = 'xxx/Preproc_HAPPE/A_Power_data';

for ii = 1:height(HAPPE_MOI_table)
    ID = HAPPE_MOI_table.ID{ii};
    if exist(strcat(Power_folder, '/',ID,'_Power_data.mat'),'file')
        load(strcat(Power_folder, '/',ID,'_Power_data.mat'))
    
        Freqs = 1:1:32;
        % across all trials
        if ~isempty(Power_all)
            Ntrls_tot = size(Power_all.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_all = mean(squeeze(mean(Power_all.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_all = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_all = nan(1,32);
            PowB_all = nan(1,32);
            Nsplit_all = 0;  
        end

        % social trials
        if ~isempty(Power_soc)
            Ntrls_tot = size(Power_soc.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_soc = mean(squeeze(mean(Power_soc.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_soc = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_soc = nan(1,32);
            PowB_soc = nan(1,32);
            Nsplit_soc = 0;  
        end

        % toy trials
        if ~isempty(Power_toy)
            Ntrls_tot = size(Power_toy.logpowspctrm,1);
            if mod(Ntrls_tot,2) == 1
                Nmax = Ntrls_tot - 1;
            else
                Nmax = Ntrls_tot;
            end
            Ind_A = 1:2:(Nmax-1);
            Ind_B = 2:2:Nmax;
            PowA_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_A,:,:),1)),1);
            PowB_toy = mean(squeeze(mean(Power_toy.logpowspctrm(Ind_B,:,:),1)),1);
            Nsplit_toy = Nmax/2;
            % clean up 
            clear Ind_A Ind_B Nmax Ntrls_tot
        else
            PowA_toy = nan(1,32);
            PowB_toy = nan(1,32);
            Nsplit_toy = 0;  
        end
        
    else
        Freqs = nan(1,32);
        Nsplit_all = 0; PowA_all = nan(1,32); PowB_all = nan(1,32);
        Nsplit_soc = 0; PowA_soc = nan(1,32); PowB_soc = nan(1,32);
        Nsplit_toy = 0; PowA_toy = nan(1,32); PowB_toy = nan(1,32);
    end
    
    % add to tracking table
        if ii == 1
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPE_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_HAPPE
            save('HAPPE_SplitHalf_reliability.mat','HAPPE_SplitHalf')
            clear HAPPE_SplitHalf
        else
            cd xxx/Preproc_HAPPE
            load HAPPE_SplitHalf_reliability.mat
            data_table_newrow = table({ID}, {Freqs}, {Nsplit_all}, {PowA_all}, {PowB_all},...
                {Nsplit_soc}, {PowA_soc}, {PowB_soc}, {Nsplit_toy}, {PowA_toy}, {PowB_toy});
            data_tabale_newrow.Properties.VariableNames = {'ID', 'Freqs', 'Nsplit_all', 'PowA_all', 'PowB_all',...
                'Nsplit_soc', 'PowA_soc', 'PowB_soc', 'Nsplit_toy', 'PowA_toy', 'PowB_toy'};
            HAPPE_SplitHalf(ii,:) = data_table_newrow;
            cd xxx/Preproc_HAPPE
            save('HAPPE_SplitHalf_reliability.mat','HAPPE_SplitHalf')
            clear HAPPE_SplitHalf
        end

end
clear HAPPE_MOI_table

%% Correlations for split half reliability

Freqs = 1:1:32;
Delta_ind = [find(Freqs == 2),find(Freqs == 3)];
Theta_ind = [find(Freqs == 4),find(Freqs == 6)];
Alpha_ind = [find(Freqs == 7),find(Freqs == 12)];
Beta_ind = [find(Freqs == 13),find(Freqs == 30)];

% Check Ntrls
cd xxx/PreprocData_Manual
load Manual_SplitHalf_reliability.mat
cd xxx/Preproc_MADE
load MADE_SplitHalf_reliability.mat
cd xxx/Preproc_BOND
load BOND_SplitHalf_reliability.mat
cd xxx/Preproc_HAPPE
load HAPPE_SplitHalf_reliability.mat
    
% set those below threshold to NaNs
Ntrls = [cell2mat(Manual_SplitHalf.Var3) cell2mat(Manual_SplitHalf.Var6) cell2mat(Manual_SplitHalf.Var9)...
    cell2mat(MADE_SplitHalf.Var3) cell2mat(MADE_SplitHalf.Var6) cell2mat(MADE_SplitHalf.Var9) ...
    cell2mat(BOND_SplitHalf.Var3) cell2mat(BOND_SplitHalf.Var6) cell2mat(BOND_SplitHalf.Var9) ...
    cell2mat(HAPPE_SplitHalf.Var3) cell2mat(HAPPE_SplitHalf.Var6) cell2mat(HAPPE_SplitHalf.Var9)];
Ind_incl_SH = find(sum(Ntrls >= 20, 2) ==12);

Ind_incl_SH_fc = find(sum(Ntrls >= 90, 2) ==12);

% all conditions
N113_Manual_PowA_spectra = cell2mat(table2array(Manual_SplitHalf(Ind_incl_SH, 4)));
N113_Manual_PowB_spectra = cell2mat(table2array(Manual_SplitHalf(Ind_incl_SH, 5)));
N113_MADE_PowA_spectra = cell2mat(table2array(MADE_SplitHalf(Ind_incl_SH, 4)));
N113_MADE_PowB_spectra = cell2mat(table2array(MADE_SplitHalf(Ind_incl_SH, 5)));
N113_BOND_PowA_spectra = cell2mat(table2array(BOND_SplitHalf(Ind_incl_SH, 4)));
N113_BOND_PowB_spectra = cell2mat(table2array(BOND_SplitHalf(Ind_incl_SH, 5)));
N113_HAPPE_PowA_spectra = cell2mat(table2array(HAPPE_SplitHalf(Ind_incl_SH, 4)));
N113_HAPPE_PowB_spectra = cell2mat(table2array(HAPPE_SplitHalf(Ind_incl_SH, 5)));


% Split-Half Reliability: A vs B - all trials
SH_all_r_vals = zeros(4,4);
SH_all_p_vals = zeros(4,4);

% Roche1 
% delta
[r,p] = corrcoef([mean(N113_Manual_PowA_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) mean(N113_Manual_PowB_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)]); 
SH_all_r_vals(1,1) = r(1,2); SH_all_p_vals(1,1) = p(1,2); clear r p
% theta
[r,p] = corrcoef([mean(N113_Manual_PowA_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) mean(N113_Manual_PowB_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)]); 
SH_all_r_vals(1,2) = r(1,2); SH_all_p_vals(1,2) = p(1,2); clear r p
% alpha
[r,p] = corrcoef([mean(N113_Manual_PowA_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) mean(N113_Manual_PowB_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)]); 
SH_all_r_vals(1,3) = r(1,2); SH_all_p_vals(1,3) = p(1,2); clear r p
% beta
[r,p] = corrcoef([mean(N113_Manual_PowA_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) mean(N113_Manual_PowB_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)]); 
SH_all_r_vals(1,4) = r(1,2); SH_all_p_vals(1,4) = p(1,2); clear r p

% MADE 
% delta
[r,p] = corrcoef([mean(N113_MADE_PowA_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) mean(N113_MADE_PowB_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)]); 
SH_all_r_vals(2,1) = r(1,2); SH_all_p_vals(2,1) = p(1,2); clear r p
% theta
[r,p] = corrcoef([mean(N113_MADE_PowA_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) mean(N113_MADE_PowB_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)]); 
SH_all_r_vals(2,2) = r(1,2); SH_all_p_vals(2,2) = p(1,2); clear r p
% alpha
[r,p] = corrcoef([mean(N113_MADE_PowA_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) mean(N113_MADE_PowB_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)]); 
SH_all_r_vals(2,3) = r(1,2); SH_all_p_vals(2,3) = p(1,2); clear r p
% beta
[r,p] = corrcoef([mean(N113_MADE_PowA_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) mean(N113_MADE_PowB_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)]); 
SH_all_r_vals(2,4) = r(1,2); SH_all_p_vals(2,4) = p(1,2); clear r p

% BOND
% delta
[r,p] = corrcoef([mean(N113_BOND_PowA_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) mean(N113_BOND_PowB_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)]); 
SH_all_r_vals(3,1) = r(1,2); SH_all_p_vals(3,1) = p(1,2); clear r p
% theta
[r,p] = corrcoef([mean(N113_BOND_PowA_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) mean(N113_BOND_PowB_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)]); 
SH_all_r_vals(3,2) = r(1,2); SH_all_p_vals(3,2) = p(1,2); clear r p
% alpha
[r,p] = corrcoef([mean(N113_BOND_PowA_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) mean(N113_BOND_PowB_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)]); 
SH_all_r_vals(3,3) = r(1,2); SH_all_p_vals(3,3) = p(1,2); clear r p
% beta
[r,p] = corrcoef([mean(N113_BOND_PowA_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) mean(N113_BOND_PowB_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)]); 
SH_all_r_vals(3,4) = r(1,2); SH_all_p_vals(3,4) = p(1,2); clear r p

% HAPPE
% delta
[r,p] = corrcoef([mean(N113_HAPPE_PowA_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) mean(N113_HAPPE_PowB_spectra(:,[Delta_ind(1):Delta_ind(2)]),2)]); 
SH_all_r_vals(4,1) = r(1,2); SH_all_p_vals(4,1) = p(1,2); clear r p
% theta
[r,p] = corrcoef([mean(N113_HAPPE_PowA_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) mean(N113_HAPPE_PowB_spectra(:,[Theta_ind(1):Theta_ind(2)]),2)]); 
SH_all_r_vals(4,2) = r(1,2); SH_all_p_vals(4,2) = p(1,2); clear r p
% alpha
[r,p] = corrcoef([mean(N113_HAPPE_PowA_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) mean(N113_HAPPE_PowB_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2)]); 
SH_all_r_vals(4,3) = r(1,2); SH_all_p_vals(4,3) = p(1,2); clear r p
% beta
[r,p] = corrcoef([mean(N113_HAPPE_PowA_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) mean(N113_HAPPE_PowB_spectra(:,[Beta_ind(1):Beta_ind(2)]),2)]); 
SH_all_r_vals(4,4) = r(1,2); SH_all_p_vals(4,4) = p(1,2); clear r p




%% Split-half reliability: A vs B condition differences
N113_Roche1_PowAs_spectra = cell2mat(table2array(Manual_SplitHalf(Ind_incl_SH, 7)));
N113_Roche1_PowAt_spectra = cell2mat(table2array(Manual_SplitHalf(Ind_incl_SH, 10)));
N113_Roche1_PowBs_spectra = cell2mat(table2array(Manual_SplitHalf(Ind_incl_SH, 8)));
N113_Roche1_PowBt_spectra = cell2mat(table2array(Manual_SplitHalf(Ind_incl_SH, 11)));

N113_MADE_PowAs_spectra = cell2mat(table2array(MADE_SplitHalf(Ind_incl_SH, 7)));
N113_MADE_PowAt_spectra = cell2mat(table2array(MADE_SplitHalf(Ind_incl_SH, 10)));
N113_MADE_PowBs_spectra = cell2mat(table2array(MADE_SplitHalf(Ind_incl_SH, 8)));
N113_MADE_PowBt_spectra = cell2mat(table2array(MADE_SplitHalf(Ind_incl_SH, 11)));

N113_BOND_PowAs_spectra = cell2mat(table2array(BOND_SplitHalf(Ind_incl_SH, 7)));
N113_BOND_PowAt_spectra = cell2mat(table2array(BOND_SplitHalf(Ind_incl_SH, 10)));
N113_BOND_PowBs_spectra = cell2mat(table2array(BOND_SplitHalf(Ind_incl_SH, 8)));
N113_BOND_PowBt_spectra = cell2mat(table2array(BOND_SplitHalf(Ind_incl_SH, 11)));

N113_HAPPE_PowAs_spectra = cell2mat(table2array(HAPPE_SplitHalf(Ind_incl_SH, 7)));
N113_HAPPE_PowAt_spectra = cell2mat(table2array(HAPPE_SplitHalf(Ind_incl_SH, 10)));
N113_HAPPE_PowBs_spectra = cell2mat(table2array(HAPPE_SplitHalf(Ind_incl_SH, 8)));
N113_HAPPE_PowBt_spectra = cell2mat(table2array(HAPPE_SplitHalf(Ind_incl_SH, 11)));


% Split-Half Reliability: A vs B - condition differences
SH_conddiff_r_vals = zeros(4,4);
SH_conddiff_p_vals = zeros(4,4);

% Roche1 
% delta 
ConDiff_A = mean(N113_Roche1_PowAs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_Roche1_PowAt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
ConDiff_B = mean(N113_Roche1_PowBs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_Roche1_PowBt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(1,1) = r(1,2); SH_conddiff_p_vals(1,1) = p(1,2); clear r p ConDiff_A ConDiff_B
% theta
ConDiff_A = mean(N113_Roche1_PowAs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_Roche1_PowAt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
ConDiff_B = mean(N113_Roche1_PowBs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_Roche1_PowBt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(1,2) = r(1,2); SH_conddiff_p_vals(1,2) = p(1,2); clear r p ConDiff_A ConDiff_B
% alpha
ConDiff_A = mean(N113_Roche1_PowAs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_Roche1_PowAt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
ConDiff_B = mean(N113_Roche1_PowBs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_Roche1_PowBt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(1,3) = r(1,2); SH_conddiff_p_vals(1,3) = p(1,2); clear r p ConDiff_A ConDiff_B
% beta
ConDiff_A = mean(N113_Roche1_PowAs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_Roche1_PowAt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
ConDiff_B = mean(N113_Roche1_PowBs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_Roche1_PowBt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(1,4) = r(1,2); SH_conddiff_p_vals(1,4) = p(1,2); clear r p ConDiff_A ConDiff_B

% MADE 
% delta 
ConDiff_A = mean(N113_MADE_PowAs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_MADE_PowAt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
ConDiff_B = mean(N113_MADE_PowBs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_MADE_PowBt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(2,1) = r(1,2); SH_conddiff_p_vals(2,1) = p(1,2); clear r p ConDiff_A ConDiff_B
% theta
ConDiff_A = mean(N113_MADE_PowAs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_MADE_PowAt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
ConDiff_B = mean(N113_MADE_PowBs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_MADE_PowBt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(2,2) = r(1,2); SH_conddiff_p_vals(2,2) = p(1,2); clear r p ConDiff_A ConDiff_B
% alpha
ConDiff_A = mean(N113_MADE_PowAs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_MADE_PowAt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
ConDiff_B = mean(N113_MADE_PowBs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_MADE_PowBt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(2,3) = r(1,2); SH_conddiff_p_vals(2,3) = p(1,2); clear r p ConDiff_A ConDiff_B
% beta
ConDiff_A = mean(N113_MADE_PowAs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_MADE_PowAt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
ConDiff_B = mean(N113_MADE_PowBs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_MADE_PowBt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(2,4) = r(1,2); SH_conddiff_p_vals(2,4) = p(1,2); clear r p ConDiff_A ConDiff_B

% BOND
% delta 
ConDiff_A = mean(N113_BOND_PowAs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_BOND_PowAt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
ConDiff_B = mean(N113_BOND_PowBs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_BOND_PowBt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(3,1) = r(1,2); SH_conddiff_p_vals(3,1) = p(1,2); clear r p ConDiff_A ConDiff_B
% theta
ConDiff_A = mean(N113_BOND_PowAs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_BOND_PowAt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
ConDiff_B = mean(N113_BOND_PowBs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_BOND_PowBt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(3,2) = r(1,2); SH_conddiff_p_vals(3,2) = p(1,2); clear r p ConDiff_A ConDiff_B
% alpha
ConDiff_A = mean(N113_BOND_PowAs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_BOND_PowAt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
ConDiff_B = mean(N113_BOND_PowBs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_BOND_PowBt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(3,3) = r(1,2); SH_conddiff_p_vals(3,3) = p(1,2); clear r p ConDiff_A ConDiff_B
% beta
ConDiff_A = mean(N113_BOND_PowAs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_BOND_PowAt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
ConDiff_B = mean(N113_BOND_PowBs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_BOND_PowBt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(3,4) = r(1,2); SH_conddiff_p_vals(3,4) = p(1,2); clear r p ConDiff_A ConDiff_B

% HAPPE
% delta 
ConDiff_A = mean(N113_HAPPE_PowAs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_HAPPE_PowAt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
ConDiff_B = mean(N113_HAPPE_PowBs_spectra(:,[Delta_ind(1):Delta_ind(2)]),2) - mean(N113_HAPPE_PowBt_spectra(:,[Delta_ind(1):Delta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(4,1) = r(1,2); SH_conddiff_p_vals(4,1) = p(1,2); clear r p ConDiff_A ConDiff_B
% theta
ConDiff_A = mean(N113_HAPPE_PowAs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_HAPPE_PowAt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
ConDiff_B = mean(N113_HAPPE_PowBs_spectra(:,[Theta_ind(1):Theta_ind(2)]),2) - mean(N113_HAPPE_PowBt_spectra(:,[Theta_ind(1):Theta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(4,2) = r(1,2); SH_conddiff_p_vals(4,2) = p(1,2); clear r p ConDiff_A ConDiff_B
% alpha
ConDiff_A = mean(N113_HAPPE_PowAs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_HAPPE_PowAt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
ConDiff_B = mean(N113_HAPPE_PowBs_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2) - mean(N113_HAPPE_PowBt_spectra(:,[Alpha_ind(1):Alpha_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(4,3) = r(1,2); SH_conddiff_p_vals(4,3) = p(1,2); clear r p ConDiff_A ConDiff_B
% beta
ConDiff_A = mean(N113_HAPPE_PowAs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_HAPPE_PowAt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
ConDiff_B = mean(N113_HAPPE_PowBs_spectra(:,[Beta_ind(1):Beta_ind(2)]),2) - mean(N113_HAPPE_PowBt_spectra(:,[Beta_ind(1):Beta_ind(2)]),2);
[r,p] = corrcoef([ConDiff_A ConDiff_B]); 
SH_conddiff_r_vals(4,4) = r(1,2); SH_conddiff_p_vals(4,4) = p(1,2); clear r p ConDiff_A ConDiff_B



%% save data
cd xxx/Figs_v3_Mar23/data_csv
% all trials
writematrix(SH_all_r_vals, 'SplitHalf_all_rvals.csv')
writematrix(SH_all_p_vals, 'SplitHalf_all_pvals.csv')
% condition differences
writematrix(SH_conddiff_r_vals, 'SplitHalf_diff_rvals.csv')
writematrix(SH_conddiff_p_vals, 'SplitHalf_diff_pvals.csv')

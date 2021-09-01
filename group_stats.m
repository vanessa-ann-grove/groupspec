%% VIPA rsEEG Static Spectral Analysis- Group Comparison
%
% Vanessa Grove, MSc
% 
% The following script is designed for the basic analysis of resting state
% electrophysiogical data between groups.
%
% To use code: For UEA analysis- Edit section 1a of code according to 
% analysis being run. For UEA:UoL analysis- specify 'UoL' as analysis_type
% variable and run code as is. If just looking for result visualisation,
% set visualisation_only variable to '1'. 
%
%% Specify Type of Analysis

analysis_type = 1; %1 = UEA, 2 = UoL; %EDIT
visualisation_only = 1; %0 if running statistics %EDIT
    
%% Section 1a: Data Parameters- UEA Analysis 
if analysis_type == 1

% Title of analysis (to be displayed on figures)
fig_title = 'FMS (n = 23) v. HC (n = 14) Baseline';

% Data to be compared
dataset1 = 'APPENDED_Pain'; %EDIT
dataset2 = 'APPENDED_ctl'; %EDIT
existing_var = 'FMSvHC_baseline'; %EDIT
Exp = [3,3]; %EDIT
stats_title = {'FMS:HC baseline stats'};

Comparison_Type = 1; %EDIT: (1 = baseline, 2 = post-VR) 


% Statistics Parameters
d1 = 1; % Baseline dataset (default = 1)
d2 = 2; 
s1 = 3; % Number of row in gp_compare for statistics
test_stat = 'Mann Whitney U/ Rank Sum';
uncorrected_pval = .05;
critical_pval = uncorrected_pval; %/ 7; 
N = 5000;

%Locate appropriate MATLAB Path
testpath = 'S:\\';
if exist([testpath 'VIPA Study\EEG Data\MATLAB\\experiment_ids_order.mat'], 'file')
    path1 = testpath;
else
    path1 = '\\ueahome\CfE-Research\\';
end

%Specify output location or reload existing gp_compare variable
if visualisation_only == 0
gp_compare = struct([]);
gp_compare(d1).id = dataset1;
gp_compare(d2).id = dataset2;  
gp_compare(s1).id = stats_title;
gp_compare(s1).Test_Stat = test_stat;
gp_compare(s1).Critical_Pval = critical_pval;
gp_compare(s1).Npermutes = N;
else if visualisation_only == 1
    load([path1 'VIPA Study\EEG Data\MATLAB\gp_compare\gp_compare_', existing_var, '.mat']);
    end
end


%% Section 1b: Data Setup

%Load desired data and store power values in dB
if visualisation_only == 0
all_data = {dataset1, dataset2};

for i = 1:length(all_data)
    path2 = [path1 'VIPA Study\EEG Data\MATLAB\E' num2str(Exp(i)) '_fft_data\'];
    dataload = [path2 'fft_data_' cell2mat(all_data(i)) '.mat'];
    load(dataload)
    
    freq_bounds = dsearchn(fft_data(1).frex, [2 45]');
    gp_compare(1).frex = fft_data(1).frex(freq_bounds(1):freq_bounds(2));
    
    if Comparison_Type == 1 
        a = 1;
    else if Comparison_Type == 2 && Exp(i) == 3
            a = 6;
        else if Comparison_Type == 2 && Exp(i) == 2
                a = 4;
            end
        end
    end
    
    powvals = fft_data(a).allpower; 
    gp_compare(i).powvals = powvals(:,freq_bounds(1):freq_bounds(2));
    gp_compare(i).scalpavg = fft_data(a).scalpavg;
    gp_compare(i).scalppow = mean(gp_compare(i).powvals,1);
end
end

%% Section 2: Data Parameters (UoL analysis)

%{
%% DO NOT UNCOMMENT! %%
data_title = {'UEA FMS', 'UEA HC', 'UoL FMS', 'UoL HC'}; 
stats_title = {'UEA FMS:HC', 'UoL FMS:HC', 'FMS UEA:UoL' 'HC UEA:UoL'};
test_stat = 'Mann Whitney U/ Rank Sum';
uncorrected_pval = .05;
critical_pval = uncorrected_pval; %/ 7; %Bonferonni correction for 7 freq bands
N = 5000;

gp_compare = struct([]);
for i = 1:length(data_title)
    gp_compare(i).id = data_title(i);
    if i == 1
        load('fft_data_APPENDED_FMS.mat')
        freq_bounds = dsearchn(fft_data(1).frex, [2 45]');
        gp_compare(1).frex = fft_data(1).frex(freq_bounds(1):freq_bounds(2));
    else if i == 2
            load('fft_data_APPENDED_ctl.mat')
        else if i == 3 || 4
             load('fft_data_UoL.mat')
            end
        end
    end
  
    if i == 1 || 2 || 3
       powvals = fft_data(1).allpower; 
       gp_compare(i).powvals = powvals(:,freq_bounds(1):freq_bounds(2));
       gp_compare(i).scalpavg = fft_data(1).scalpavg;
       gp_compare(i).scalppow = mean(gp_compare(i).powvals,1);
    else if i == 4
       powvals = fft_data(2).allpower; 
       gp_compare(i).powvals = powvals(:,freq_bounds(1):freq_bounds(2));
       gp_compare(i).scalpavg = fft_data(2).scalpavg;
       gp_compare(i).scalppow = mean(gp_compare(i).powvals,1);
        end
    end
end
%}
else if analysis_type == 2

testpath = 'S:\\';
if exist([testpath 'VIPA Study\EEG Data\MATLAB\\experiment_ids_order.mat'], 'file')
    path1 = testpath;
else
    path1 = '\\ueahome\CfE-Research\\';
end

load([path1 'VIPA Study\EEG Data\MATLAB\gp_compare\gp_compare_UoL.mat'])
baseline_data = [1 3 1 2 1 3];
compare_data = [2 4 3 4 4 2];
test_stat = 'Mann Whitney U/ Rank Sum';
uncorrected_pval = .05;
critical_pval = uncorrected_pval; %/ 7;
N = 5000;
stats_title = {'UEA FMS:HC', 'UoL FMS:HC', 'FMS UEA:UoL' 'HC UEA:UoL' 'UEA FMS:UoL HC' 'UoL FMS:UEA HC'};
    end
end

%% Section 3: Statistics

for i = 1:length(stats_title)
%Frequency band boundries
hz = gp_compare(1).frex;
freq_bounds = dsearchn(hz, [2 4 8 12 16 23 30 40]');
 if exist('baseline_data', 'var') 
    d1 = baseline_data(i);
    d2 = compare_data(i);
    s1 = i + 4;
 end
 
%Data
if visualisation_only == 0
    disp('Running scalp statistics........')
preInt = gp_compare(d1).powvals; 
postInt = gp_compare(d2).powvals;
ncmp = size(preInt,1); %number of channels
gp_compare(s1).id = stats_title(i);
gp_compare(s1).Test_Stat = test_stat;
gp_compare(s1).Critical_Pval = critical_pval;
gp_compare(s1).Npermutes = N;

for k=2:length(freq_bounds)

%% Section 4: Statistical Analysis (Whole Scalp, Frequency Bands)
        %Step 1: Select data from freq band
        preData  = preInt(:,freq_bounds(k-1):freq_bounds(k));
            n1 = size(preData,2);
        postData = postInt(:,freq_bounds(k-1):freq_bounds(k));
            n2 = size(postData,2);
        allData = horzcat(preData,postData)';
        truelabels = cat(1,ones(n1,1),2*ones(n2,1));
    
        % Step 2: Compute true test statistic 
            % Wilcoxon signed rank test
                gp_compare(d1).bandpow(k-1,:) = mean(allData(truelabels==1,:),1);
                gp_compare(d2).bandpow(k-1,:) = mean(allData(truelabels==2,:),1);
                
                [p h stats] = ranksum(mean(allData(truelabels==1,:),1),...
                    mean(allData(truelabels==2,:),1));
                true_stat = p;
                gp_compare(s1).scalpstat(k-1,:) = true_stat;
            
        % Step 3: Permute
        for permi=1:N
            % Step 4: Shuffle data labels
            shuflabels = truelabels(randperm(n1+n2));
            %Step 5: Compute test stat with shuffled labels and store
                [p h stats] = ranksum(mean(allData(shuflabels==1,:),1),...
                    mean(allData(shuflabels==2,:),1));
                gp_compare(s1).scalpperm(k-1,permi) = p;
        end
        
        %Step 6: Calculate p value and create new variable to only store significant test stats
        scalp_stat = gp_compare(s1).scalpstat(k-1);
        abs_xtreme = gp_compare(s1).scalpperm(k-1,:);
        pval = sum(abs_xtreme < scalp_stat,2) / N;
        gp_compare(s1).scalppvals(k-1) = pval;
        trialsig = pval  < critical_pval;
        gp_compare(s1).scalpsig(k-1) = trialsig; 

%% Section 5: Statistical Analysis (Single channel: concatenated data)
% Only to be performed if freq band is significant at scalp level
trialsig = gp_compare(s1).scalpsig(k-1);
if trialsig == 0
        gp_compare(s1).truestat(k-1,:) = 1+zeros(1,ncmp);
        gp_compare(s1).chanperm(k-1,:) = zeros(1,N);
        gp_compare(s1).chan_pvals(k-1,:) = 1+zeros(1,ncmp);
        gp_compare(s1).sigdif(k-1,:) = zeros(1,ncmp);
        
else if trialsig == 1
        disp('Significant difference found at scalp level! Running channel statistics......')
        % Step 7: Compute true test statistic 
            %Sign Rank Test for independent samples
            for chani = 1:ncmp
                [p h stats] = ranksum(allData(truelabels==1,chani),allData(truelabels==2,chani));
               gp_compare(s1).truestat(k-1,chani) = p;  %
            end
    
        % Step 8: Permute and store min value
        for permi=1:N
            % Step 9: Shuffle data labels
            shuflabels = truelabels(randperm(n1+n2));
            %Step 10: Compute test stat with shuffled labels and store
                perm_stat = zeros(1,ncmp);
                for chani = 1:ncmp
                    [p h stats] = ranksum(allData(shuflabels==1,chani),allData(shuflabels==2,chani));
                    perm_stat(chani) = p;
                end
                 gp_compare(s1).chanperm(k-1,permi) = min(perm_stat);
        end

        %Step 11: Calculate p value and create new variable to demonstrate significant results
        chan_stat = gp_compare(s1).truestat(k-1,:);
        abs_xtreme =  gp_compare(s1).chanperm(k-1,:)';
        pval = sum(abs_xtreme < chan_stat) / N;
        gp_compare(s1).chan_pvals(k-1,:) = pval;
        trialsig = pval < critical_pval;
        gp_compare(s1).sigdif(k-1,:) = trialsig;

        end
    end
end

gp_compare(d1).scalpstat = mean(gp_compare(d1).bandpow,2);
gp_compare(d2).scalpstat = mean(gp_compare(d2).bandpow,2);
end

%% Section 6: Visualisation- Spectrum and Bar Graph of Scalp Power, Significant Channels

figure(i),clf
labs = {'Delta (2-4 Hz)' 'Theta (4-8 Hz)' 'Alpha (8-12 Hz)'...
    'Beta 1 (12-16 Hz)' 'Beta 2 (16-23 Hz)' 'Beta 3 (23-30 Hz)' 'Gamma (30-40 Hz)'};
leg = {gp_compare(d1).id, gp_compare(d2).id};

subplot(312)
combined_dat5 = [gp_compare(d1).scalpavg;gp_compare(d2).scalpavg]';
b = bar(combined_dat5);
n = 2;
hold on
datai = s1;
sigData2 = gp_compare(datai).scalpsig';
    for freqi=1:7
        maxval = combined_dat5(freqi,1);
        if sigData2(freqi) == 1
           plot(freqi-0.15, 1.1*maxval ,'*k')  
        end
    end
xticklabels(labs)
ylabel('Average Power (\muV^2)')
legend(leg)
title('Average Scalp Power')

% Full Spectra
hold on
freq_bounds = dsearchn(gp_compare(1).frex, [2 45]');
frex = gp_compare(1).frex(freq_bounds(1):freq_bounds(2));
subplot(313)
for l = [d1,d2]
    spectra = mean(gp_compare(l).powvals,1);
    p = plot(frex,spectra);
       p.LineWidth = 1.5;
    hold on
end
ax = gca;
ax.XLim = [2 45];
ax.XTick = [2 4 8 12 16 23 30 45];
for o = 1:length(ax.XTick)
    xline(ax.XTick(o))
end
xlabel('Frequency (Hz)')
ylabel('Log Power Spectral Density 10*log_{10}(\muV^{2}/Hz)')
title(stats_title(i))
legend(leg)
hold on

%% Significant Channel Visualisation

figtitle = 'Group Differences By Channel';
% Define Colour Limits
min_max = zeros(2,7);
for freqi = 1:7
    freq_xtrems = zeros(2,1);    
        %dat = gp_compare(datai).truestat(freqi,:); %tstat
        dat = gp_compare(d2).bandpow(freqi,:)-gp_compare(d1).bandpow(freqi,:); 
        freq_xtrems(1,(datai-1-1)) = min(dat);
        freq_xtrems(2,(datai-1-1)) = max(dat);
        min_max(1,freqi) = min(freq_xtrems(1,:));
        min_max(2,freqi) = max(freq_xtrems(2,:));
end 

sigData = gp_compare(datai).sigdif';
sigData2 = gp_compare(datai).scalpsig';

     % tstat as heatmap data
        %data2plot = gp_compare(datai).truestat;
       
     % channel power difference (HC-FMS) as heatmap data   
        preData = gp_compare(d1).bandpow;
        postData = gp_compare(d2).bandpow;
        data2plot = postData - preData;
        load([path1 'VIPA Study\EEG Data\MATLAB\E' num2str(Exp(i)) '_fft_data\fft_data_V2.mat'])
        chaninfo = fft_data(1).chaninfo;
        clear fft_data
        
    for freqi = 1:size(data2plot,1)
        mi = min_max(1,freqi);
        ma = min_max(2,freqi);
        clear sigchanlocs
        for chani = 1:64
            if sigData(chani,freqi)== 1
                sigchanlocs(chani) = chaninfo(chani);
            end
        end
   
        if exist('sigchanlocs', 'var') && sigData2(freqi)==1
            subplot(3,7,(freqi))
            topoplot(data2plot(freqi,:),chaninfo,'electrodes','off','style','map','maplimits',[mi ma]);
            hold on
            topoplot([],sigchanlocs,  'electrodes', 'on');
            colorbar
            title(labs(freqi))
            hold off
        else
            subplot(3,7,(freqi))
            topoplot(data2plot(freqi,:),chaninfo, 'electrodes','off','style','map','maplimits',[mi ma]);
            colorbar
            title(labs(freqi))
        end
    end
set(gcf, 'name', cell2mat(stats_title(i)))

%% Print Results for User View

freq_bands = {'Delta band' 'Theta band' 'Alpha band' 'Beta1 band'...
    'Beta2 band' 'Beta3 band' 'Gamma band'};

disp([char(leg(1)) ' v. ' char(leg(2)) ' analysis complete!'])
all_sig = gp_compare(s1).scalpsig;

if any(all_sig,'all')
    disp('Significant Results Observed:')
    for f = 1:length(freq_bands)
       if all_sig(f) == 1 
           disp([char(freq_bands(f)) ' is significantly different!'])
       end
    end
else 
   disp('No significant results observed')
    end
end
disp('All Analyses Complete. Dont forget to save gp_compare variable!!!!')

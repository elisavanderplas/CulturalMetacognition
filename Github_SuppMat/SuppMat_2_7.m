%% Correlation plot between social and non-social postxacc impact & traits (Supplementary Material 2.8)

close all;clear all;
culture = {'PKU';'UCL'};
fs = filesep;

surveydat = []; 

%% now load the beta coefficients for each individuals impact of pre, post and post x accuracy interaction in Exp 1 (made with EXP2_ExtractSumStats.R). 
baseDir2 = ['~' fs 'Dropbox' fs 'PKU_collaboration' fs 'Github' fs 'DATA' fs 'EXP2'] ;
filename2 = 'regression_betas_IDs_EXP2.csv';    % use for fMRI
cd(baseDir2);
temp= readtable(filename2,'TreatAsEmpty',{'.','NA'});

% specify traits
b_s_post= table2array(temp(:,5));
b_ns_post= table2array(temp(:,7));
b_s_accpost = table2array(temp(:,8));
b_ns_accpost = table2array(temp(:,9));

%% now load the questionnaire data
for n = 1:length(culture)
    nat = culture{n};
    
    dirData = [baseDir2  fs nat '_data' fs nat '_data' fs nat '_betas' fs];
    filename = '_subject_log';
    
    datafile = [nat filename '.csv']; %betas
    cd(dirData);
    surveydat = [surveydat; readtable(datafile)];
    
end

%% Correlation plot between social and non-social postxacc impact & traits (Supplementary Material 2.8)
%holism
[RHO,PVAL] = corr(zscore(b_ns_accpost),zscore(surveydat.constr2))
[RHO,PVAL] = corr(zscore(b_s_accpost),zscore(surveydat.constr2))
%IQ
[RHO,PVAL] = corr(zscore(b_ns_accpost),zscore(surveydat.IQ))
[RHO,PVAL] = corr(zscore(b_s_accpost),zscore(surveydat.IQ))
%BCIS
[RHO,PVAL] = corr(zscore(b_ns_accpost),zscore(surveydat.BCIS))
[RHO,PVAL] = corr(zscore(b_s_accpost),zscore(surveydat.BCIS))
%PKU separately
[RHO,PVAL] = corr(zscore(b_ns_accpost(1:53)),zscore(surveydat.BCIS(1:53)))
[RHO,PVAL] = corr(zscore(b_s_accpost(1:53)),zscore(surveydat.BCIS(1:53)))
%UCL separately
[RHO,PVAL] = corr(zscore(b_ns_accpost(54:end)),zscore(surveydat.BCIS(54:end)))
[RHO,PVAL] = corr(zscore(b_s_accpost(54:end)),zscore(surveydat.BCIS(54:end)))
%BCIS_SC
[RHO,PVAL] = corr(zscore(b_ns_accpost),zscore(surveydat.BCIS_SC))
[RHO,PVAL] = corr(zscore(b_s_accpost),zscore(surveydat.BCIS_SC))
%BCIS_SR
[RHO,PVAL] = corr(zscore(b_ns_accpost),zscore(surveydat.BCIS_SR))
[RHO,PVAL] = corr(zscore(b_s_accpost),zscore(surveydat.BCIS_SR))


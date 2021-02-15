%% Makes Fig 4 of PKU-UCL paper
% EVDP 2019 elisa.plas.18@ucl.ac.uk

culture = {'PKU';'UCL'};
fs = filesep;

%set fig properties
c.corr =  [0.082, 0.615, 0.835];
c.err = [0.835, 0.250, 0.082];
c.grey = [0.5, 0.5, 0.5];
size = 90;
ms = 20; 
axis_text = 24; 
axis_nr = 18; 
h03 = figure;
set(h03,'units','points','position',[10,10,1600,600])

%set path to used toolboxes
addpath(['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'tools' fs 'Corr_toolbox_v2']); 
addpath(['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'tools' fs 'Corr_toolbox_v2' fs 'LIBRA']);
   
%load betas Exp 1, made in R with 'EXP1_ExtractSumStats.r 
baseDir1 = ['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'DATA' fs 'EXP1'] ;
cd(baseDir1);
fn = 'regression_betas_IDs_EXP1.csv'; 
temp = readtable(fn,'TreatAsEmpty',{'.','NA'});
b_ns_accpost_exp1= table2array(temp(:,6));

%load betas Exp 2, made in R with 'EXP2_ExtractSumStats.r'
baseDir2 = ['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'DATA' fs 'EXP2'] ;
fn = 'regression_betas_IDs_EXP2.csv';  
cd(baseDir2);
temp= readtable(fn,'TreatAsEmpty',{'.','NA'});
b_s_post= table2array(temp(:,5));
b_ns_post= table2array(temp(:,7));
b_s_accpost = table2array(temp(:,8));
b_ns_accpost = table2array(temp(:,9));
    
subplot(3,2,[1,3,5])
%load betas Exp 1, made in R with 'EXP1_Analyses.R' line 170-203
for n = 1:length(culture)
    nat = culture{n};
    baseDir1 = ['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'DATA' fs 'EXP1'] ;
    dirData = [baseDir1 fs nat '_data' fs nat '_data' fs nat '_betas' fs]; 
    filename = 'regression_betas_';  
    
    for acc = 1:2
        correct = {'corr_'; 'err_'};
        suffix = correct{acc}; 
        datafile = [filename suffix nat '.csv']; 
        cd(dirData);
        dat1{n, acc} = readtable(datafile);
    end
    
    baseDir2 = ['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'DATA' fs 'EXP2'] ;
    dirData = [baseDir2  fs nat '_data' fs nat '_data' fs nat '_betas' fs];   
    %load betas Exp 2, made with 'EXP2_NS_Analyses.r' line 132-165 and 'EXP2_S_Analyses.r' line 164-210
    for acc = 1:2
        correct = {'corr_'; 'err_'};
        suffix = correct{acc}; 
        datafile = [filename 'EXP2_NS_' suffix nat '.csv']; 
        cd(dirData);
        dat2{n, acc} = readtable(datafile);
        datafile = [filename 'EXP2_S_' suffix nat '.csv']; 
        cd(dirData);
        dat3{n, acc} = readtable(datafile);
    end
    
end
grid off; box off; hold all
xpos = [1 4 7];
aa = plot(xpos,  [dat1{1,1}{3,2} dat2{1,1}{3,2} dat3{1,1}{3,2}] , 'd', 'MarkerSize',ms,  'MarkerFaceColor', c.corr, 'Linewidth', 3, 'MarkerEdgeColor', c.corr);
ab = plot(xpos,  [dat1{1,2}{3,2} dat2{1,2}{3,2} dat3{1,2}{3,2}] , 'd', 'MarkerSize', ms,'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
xpos = [2 5 8];
aa = plot(xpos,  [dat1{2,1}{3,2} dat2{2,1}{3,2} dat3{2,1}{3,2}] , 'o', 'MarkerSize',ms,  'MarkerFaceColor', c.corr, 'Linewidth', 3, 'MarkerEdgeColor', c.corr);
ab = plot(xpos,  [dat1{2,2}{3,2} dat2{2,2}{3,2} dat3{2,2}{3,2}] , 'o', 'MarkerSize', ms,'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
xpos = [1 2 4 5 7 8];
errorbar(xpos,   [dat1{1,1}{3,2} dat1{2,1}{3,2} dat2{1,1}{3,2} dat2{2,1}{3,2} dat3{1,1}{3,2} dat3{2,1}{3,2}] , [dat1{1,1}{8,2} dat1{2,1}{8,2} dat2{1,1}{8,2} dat2{2,1}{8,2} dat3{1,1}{8,2} dat3{2,1}{8,2}] , '.', 'Color', c.corr*0.2, 'LineWidth', 2);
errorbar(xpos,   [dat1{1,2}{3,2} dat1{2,2}{3,2} dat2{1,2}{3,2} dat2{2,2}{3,2} dat3{1,2}{3,2} dat3{2,2}{3,2}] , [dat1{1,2}{8,2} dat1{2,2}{8,2} dat2{1,2}{8,2} dat2{2,2}{8,2} dat3{1,2}{8,2} dat3{2,2}{8,2}], '.', 'Color', c.err*0.2, 'LineWidth', 2);

hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7);%zero line
xlabels = {'PKU', 'UCL'};

set(gca, 'XLim', [0 9], 'XTick', xpos, 'XTickLabel',['PKU  '; '  UCL'], 'YLim', [-0.5 0.2],'YTick', [-0.5:0.1:0.2], 'FontSize',axis_nr);
ylabel([{'Post-decision evidence'};{'impact on confidence (a.u.)'}], 'FontSize', axis_text);

[lgd, handle] = legend([aa, ab],{'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',1)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',axis_text);

text(0.4, -0.555, 'Experiment 1' , 'FontSize', axis_text-2)
text(0.66, -0.58, 'perceptual', 'FontSize', axis_text-4)
text(3.4, -0.555, 'Experiment 2', 'FontSize', axis_text-2)
text(3.66, -0.58,  'perceptual', 'FontSize', axis_text-4)
text(6.37 , -0.555, ' Experiment 2', 'FontSize', axis_text-2)
text(6.9, -0.58, ' social', 'FontSize', axis_text-4)
   
% Correlation plot questionnaire data and non-social post processing
subplot(3,2,[2,4,6])

scatter(zscore(b_ns_post(54:106)), zscore(b_s_post(54:106)),size,'Marker', 'o', 'MarkerEdgeColor','k', ...
    'MarkerFaceColor',[1, 1, 1],...
    'LineWidth',3)
hold all
scatter(zscore(b_ns_post(1:53)), zscore(b_s_post(1:53)), size,'Marker', 'd','MarkerEdgeColor','k', ...
    'MarkerFaceColor',[0.6 0.6 0.6],...
    'LineWidth',3)
refline; 
set(gca, 'YLim', [-1 4], 'XLim', [-1.5 4], 'YTick', [-1 0, 1, 2, 3, 4], 'FontSize', axis_nr);
ylabel([{'Impact social'};{'post-decision evidence (a.u.)'}], 'FontSize', axis_text);
xlabel(['Impact perceptual post-decision evidence (a.u.)'], 'FontSize', axis_text);
%stats
corr_UCL = robust_correlation(zscore(b_ns_post(54:106)),zscore(b_s_post(54:106)));
corr_PKU = robust_correlation(zscore(b_ns_post(1:53)),zscore(b_s_post(1:53)));


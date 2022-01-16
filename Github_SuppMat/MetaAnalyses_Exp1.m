%% Meta-analyses
% EVDP 2019 elisa.plas.18@ucl.ac.uk
fs = filesep;

culture = {'PKU', 'UCL'};
sj_mat = {[101:115, 117:141],[201:204, 206:234, 236:242]};
c.corr =  [0.082, 0.615, 0.835];
c.err = [0.835, 0.250, 0.082];

for n = 1:length(culture)
    nat = culture{n};
    baseDir =  ['~' fs 'Dropbox' fs 'Github' fs 'CulturalMetacognition' fs];
    dirData = [baseDir 'DATA' fs 'EXP1' fs nat '_data' fs nat '_data' fs];
    addpath([baseDir fs 'tools' fs 'HMeta-d' fs 'Matlab' fs]);
    
    filename = 'Data_sub_';
    subjects = sj_mat{n};
    
    for s = 1:length(subjects)
        % Load data for this subject
        datafile = [filename num2str(subjects(s)) '_2.mat'];
        cd(dirData);
        load(datafile);
        
        precoh = locDATA.dots_coherence';
        postcoh = locDATA.post_coherence';
        dir = locDATA.dots_direction/360;
        dir(dir==0.5) = -1;
        RT = log(locDATA.reaction_time_button);
        action = locDATA.button_response - 1;
        conf = locDATA.mouse_response;
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        coherence = unique(precoh);
       
        %extra compute var for meta-d'
        %%computing meta-d' across sites
        edges = [0:0.2:1];
        binned_conf = discretize(conf, edges);
        % calculate meta-d' variables
        for e = 1:length(edges) %%for all possible ratings
            nR_S1_corr{s}(e) = sum(binned_conf==e & dir==-1 & acc==1); %how often reported confrating r when dir==left & acc==1
            nR_S1_err{s}(e) = sum(binned_conf==e & dir==-1 & acc==0); %how often reported confrating r when dir==left & acc==0
            nR_S2_corr{s}(e) = sum(binned_conf==e & dir==1 & acc==1);%idem for dir==right
            nR_S2_err{s}(e) = sum(binned_conf==e & dir==1 & acc==0);
        end
        if n == 1
            metaData{1}.nR_S1{s} = [fliplr(nR_S1_corr{s}), nR_S1_err{s}];
            metaData{1}.nR_S2{s} = [fliplr(nR_S2_err{s}), nR_S2_corr{s}];
        else
            metaData{2}.nR_S1{s} = [fliplr(nR_S1_corr{s}), nR_S1_err{s}];
            metaData{2}.nR_S2{s} = [fliplr(nR_S2_err{s}), nR_S2_corr{s}];
        end
    end
end
%%compute meta-d' difference
mcmc_params.nsamples = 30000;
PKU_FIT = fit_meta_d_mcmc_group(metaData{1}.nR_S1, metaData{1}.nR_S2);
UCL_FIT = fit_meta_d_mcmc_group(metaData{2}.nR_S1, metaData{2}.nR_S2);

sampleDiff = PKU_FIT.mcmc.samples.mu_logMratio - UCL_FIT.mcmc.samples.mu_logMratio;
hdi = calc_HDI(sampleDiff(:));
fprintf(('\n Mratio group values = %.2f and %.2f'), exp(PKU_FIT.mu_logMratio), exp(UCL_FIT.mu_logMratio));
fprintf(['\n Estimated difference in Mratio between groups: ', num2str(exp(PKU_FIT.mu_logMratio) - exp(UCL_FIT.mu_logMratio))])
fprintf(['\n HDI on difference in log(Mratio): ', num2str(hdi) '\n\n'])
color_TOM = {[0.427, 0.298, 0.803], [0.960, 0.325, 0.019]};%http://doc.instantreality.org/tools/color_calculator/

% compute probability of difference
temp = sampleDiff < 0;
p_theta = (sum(temp(:) == 1))/30000;
1-p_theta

axis_text = 21
% subplot(1,2,2)
hold all
HDI = calc_HDI(exp(PKU_FIT.mcmc.samples.mu_logMratio(:)));
leg1=xline(HDI(1), '--', 'color', color_TOM{1}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{1}, 'linewidth', 2)
histogram(exp(PKU_FIT.mcmc.samples.mu_logMratio(:)),500, 'facecolor', color_TOM{1}, 'facealpha', 0.4, 'edgecolor', color_TOM{1}, 'edgealpha', 0.4);
histogram(exp(UCL_FIT.mcmc.samples.mu_logMratio(:)),500,'facecolor', color_TOM{2}, 'facealpha', 0.4, 'edgecolor', color_TOM{2}, 'edgealpha', 0.4);
HDI = calc_HDI(exp(UCL_FIT.mcmc.samples.mu_logMratio(:)));
leg2 = xline(HDI(1), '--', 'color', color_TOM{2}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{2}, 'linewidth', 2)
xline(0, '-', 'color', 'k', 'linewidth', 1)
ylabel('No. of samples', 'FontSize', axis_text);
set(gca, 'XLim', [1 2.6], 'YLim', [0 250], 'FontSize',22);
xlabel('Posterior distribution', 'FontSize',25);
set(gcf, 'color', 'w')
fakeline1 = plot(2,2000, 'color', color_TOM{1}, 'linewidth', 2);
fakeline2 = plot(2,2000, 'color', color_TOM{2}, 'linewidth', 2);

[lgd, handle] = legend([fakeline1, fakeline2], {'PKU', 'UCL'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',7)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',15);
hold off

subplot(1,2,1)  
nat = {'PKU'; 'UCL'};
correct = {'corr_'; 'err_'};
    
for n = 1:2
    for acc = 1:2
    dirData = [baseDir 'DATA' fs 'EXP1' fs nat{n} '_data' fs nat{n} '_data' fs nat{n} '_betas' fs];
    suffix = ['RT_' correct{acc} nat{n}];
        %% Load betas
        cd(dirData)
        datafile = [suffix '.csv'];
        dat{n, acc} = readtable(datafile);
    end
end
    
    ms = 22;
    axis_text = 26;
    axis_nr = 18;
    hold all
    plot(1,  [dat{1,1}{2,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', c.corr, 'Linewidth', 3, 'MarkerEdgeColor', c.corr);
    plot(1,  [dat{1,2}{2,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
    x1 = plot(2,  [dat{2,1}{2,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', c.corr,'Linewidth', 3, 'MarkerEdgeColor', c.corr);
    x2= plot(2,  [dat{2,2}{2,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
    errorbar(1:2, [dat{1,1}{2,2} dat{2,1}{2,2}], [dat{1,1}{4,2} dat{2,1}{4,2}], '.', 'Color', c.corr*0.2, 'LineWidth', 2);
    errorbar(1:2, [dat{1,2}{2,2} dat{2,2}{2,2}], [dat{1,2}{4,2} dat{2,2}{4,2}], '.', 'Color', c.err*0.2 , 'LineWidth', 2);
    
    hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7); %zeroline
    set(gca, 'XLim', [0 3], 'XTick', 1:2, 'YLim', [-0.06 0.02],'YTick', [-0.06:0.02:0.02], 'FontSize',axis_nr, 'XTickLabels', {'PKU', 'UCL'}, 'Fontsize', 18);
    ylabel([{'logRT impact'};{'on confidence (a.u.)'}], 'FontSize', axis_text);
    set(gcf, 'color', 'w')
    
    
[lgd, handle] = legend([x1, x2],{'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',1)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);
    


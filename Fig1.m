%% Plot Fig1 of PKU-UCL paper
% EVDP 2019 elisa.plas.18@ucl.ac.uk
fs = filesep;

culture = {'PKU', 'UCL'};
sj_mat = {[101:109, 111:115, 117:141],[201:204, 206:227, 229:234, 236:242]};

h01 = figure;
%set fig properties
set(h01,'units','points','position',[10,10,500,400])
c.green = {[0.380, 0.596, 0.243],[1, 1, 1]};
linestyle = {'d', 'o'};

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
        
        %index pre-decision evidence levels
        for i = 1:3
            precoh_index(locDATA.dots_coherence==coherence(i))=i;
        end
        
        %%get performance-levels across conditions
        for pre = 1:3
            accuracy(s,pre) = nanmean(acc(precoh_index == pre));
        end
        mean_acc(n,s) = nanmean(acc);
        mean_conf(n,s) = nanmean(conf);
        mean_rt(n,s) = nanmean(RT);
        
        %extra compute for meta-d'
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
    %%overlay subjects as dots
    hold all; box off
    xpos = 1:3;
    for i = 1:length(subjects)
        r = -0.1 + (0.1+0.1).*rand(1,1);
        scat = plot(xpos+r, accuracy(i,:),linestyle{n}, 'MarkerSize', 5, 'MarkerFaceColor', c.green{n}, 'MarkerEdgeColor', c.green{1}, 'LineWidth', 0.5);
        scat.Color(4) = 0.01; %set facealpha transparency
    end
    nso{n} = errorbar(xpos, nanmean(accuracy),nanstd(accuracy)./sqrt(length(subjects)), '-', 'LineWidth', 4, 'Color',[0.207, 0.388, 0.090]);
end
%set scattered linestyle for UCL
nso{2}.LineStyle = '--';

fakeline_PKU = errorbar(xpos, nanmean(accuracy)+2000, nanstd(accuracy)./sqrt(length(subjects)), '-d', 'LineWidth', 4, 'MarkerFaceColor','k', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'Color', 'k'); %for the legend
fakeline_UCL = errorbar(xpos, nanmean(accuracy)+2000, nanstd(accuracy)./sqrt(length(subjects)), '--o', 'LineWidth', 4, 'MarkerFaceColor', 'w', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'Color', 'k');

[lgd, handle] = legend([fakeline_PKU, fakeline_UCL], {' PKU', ' UCL'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'dot');
set(linehandle, 'LineWidth',8)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);

set(gca, 'YLim', [0.5 1], 'XLim', [0.5  3.5], 'XTick', xpos, 'YTick', [.5:0.1:1], 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', 20);
ylabel('Choice accuracy', 'FontSize', 30);
xlabel('Pre-decision evidence', 'FontSize', 30);

%line 221 of manuscript, basic accuracy differences
nanmean(mean_acc(1,:))%get mean/sem, 1 = PKU, 2 = UCL
nanstd(mean_acc(1,:)./sqrt(length(mean_acc(1,:))))
[H,P,CI] = ttest2(mean_acc(1,:), mean_acc(2,:));

%line 159 of manuscript, basic confidence differences
nanmean(mean_conf(1,:))%get mean/sem, 1 = PKU, 2 = UCL
nanstd(mean_conf(1,:)./sqrt(length(mean_conf(1,:))))
[H,P,CI] = ttest2(mean_conf(1,:), mean_conf(2,:));

%line 353 of manuscript, basic RT differences
nanmean(mean_rt(1,:))%get mean/sem, 1 = PKU, 2 = UCL
nanstd(mean_rt(1,:)./sqrt(length(mean_rt(1,:))))
[H,P,CI] = ttest2(mean_rt(1,:), mean_rt(2,:));

fig1b = figure; 
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
ylabel('No. of samples', 'FontSize', 18);
set(gca, 'XLim', [1 2.6], 'YLim', [0 250], 'FontSize',22);
xlabel('Posterior distribution', 'FontSize',18);
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

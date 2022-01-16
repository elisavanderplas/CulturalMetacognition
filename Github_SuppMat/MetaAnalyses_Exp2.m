%% Plot first-order performance in the social and non-social condition of Exp 2
fs = filesep;

culture = {'PKU';'UCL'};
sj_mat = {[403:418, 420:432, 434:443, 445:459]; [25:79]}; 

h02 = figure; 
set(h02,'units','points','position',[10,10,600,550])

%Colours
c.green_ns = {[0.250, 0.415, 0.145],[1, 1, 1]};
c.green_s = {[0.647, 0.835, 0.521],[1, 1, 1]};

linestyle = {'d', 'o'};

for n = 1:length(culture)
    nat = culture{n};

    baseDir =  ['~' fs 'Dropbox' fs 'Github' fs 'CulturalMetacognition' fs];
    dirData = [baseDir 'DATA' fs 'EXP2' fs nat '_data' fs nat '_data' fs]; 

    filename = 'Data_sub_'; 
    subjects = sj_mat{n};
    
    for s = 1:length(subjects)
        
        %% Load data for this subject
        datafile = [filename num2str(subjects(s)) '_2.mat'];
        cd(dirData);
        load(datafile);
        
        precoh_index = [];
        precoh = locDATA.dots_coherence';
        dir = locDATA.dots_direction/360;
        dir(dir==0.5) = -1;
        conf = locDATA.mouse_response;
        action = locDATA.button_response - 1;
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        task = locDATA.condition;
        acc_t = acc +1; 
        coherence = unique(precoh);
        RT = log(locDATA.reaction_time_button);
        
        for i = 1:3
            precoh_index(locDATA.dots_coherence==coherence(i))=i;
        end
        
        for pre = 1:3
            accuracy_social(s,pre) = nanmean(acc(precoh_index == pre & task == 1));
            accuracy_nonsocial(s,pre) = nanmean(acc(precoh_index == pre & task == 0));
        end
        
        action_adv = locDATA.a_adv -1;
        transformed_action_adv = action_adv;
        transformed_action_adv(action_adv == 0) = -1;
        acc_adv = dir == transformed_action_adv;
        acc_advt = acc_adv +1; 
        
        task = locDATA.condition; %1 = social, 0 = nonsocial
        conf_adv = locDATA.conf_adv; 
        conf_adv(conf_adv == 99) = NaN;
        agree = transformed_action_adv == transformed_action;
        
        mean_conf(n,s) = nanmean(conf);
        mean_rt(n,s) = nanmean(RT);
        
        %extra compute var for meta-d'
        %%computing meta-d' across sites
        edges = [0:0.2:1];
        binned_conf = discretize(conf, edges);
        % calculate meta-d' variables
        for e = 1:length(edges) %%for all possible ratings
            nR_S1_corr_s{s}(e) = sum(binned_conf==e & dir==-1 & acc==1 & task == 1); %how often reported confrating r when dir==left & acc==1
            nR_S1_err_s{s}(e) = sum(binned_conf==e & dir==-1 & acc==0 & task == 1); %how often reported confrating r when dir==left & acc==0
            nR_S2_corr_s{s}(e) = sum(binned_conf==e & dir==1 & acc==1 & task == 1);%idem for dir==right
            nR_S2_err_s{s}(e) = sum(binned_conf==e & dir==1 & acc==0 & task == 1);
        end
        if n == 1
            metaData_s{1}.nR_S1{s} = [fliplr(nR_S1_corr_s{s}), nR_S1_err_s{s}];
            metaData_s{1}.nR_S2{s} = [fliplr(nR_S2_err_s{s}), nR_S2_corr_s{s}];
        else
            metaData_s{2}.nR_S1{s} = [fliplr(nR_S1_corr_s{s}), nR_S1_err_s{s}];
            metaData_s{2}.nR_S2{s} = [fliplr(nR_S2_err_s{s}), nR_S2_corr_s{s}];
        end
        
           % calculate meta-d' variables
        for e = 1:length(edges) %%for all possible ratings
            nR_S1_corr_ns{s}(e) = sum(binned_conf==e & dir==-1 & acc==1 & task == 0); %how often reported confrating r when dir==left & acc==1
            nR_S1_err_ns{s}(e) = sum(binned_conf==e & dir==-1 & acc==0 & task == 0); %how often reported confrating r when dir==left & acc==0
            nR_S2_corr_ns{s}(e) = sum(binned_conf==e & dir==1 & acc==1 & task == 0);%idem for dir==right
            nR_S2_err_ns{s}(e) = sum(binned_conf==e & dir==1 & acc==0 & task == 0);
        end
        if n == 1
            metaData_ns{1}.nR_S1{s} = [fliplr(nR_S1_corr_ns{s}), nR_S1_err_ns{s}];
            metaData_ns{1}.nR_S2{s} = [fliplr(nR_S2_err_ns{s}), nR_S2_corr_ns{s}];
        else
            metaData_ns{2}.nR_S1{s} = [fliplr(nR_S1_corr_ns{s}), nR_S1_err_ns{s}];
            metaData_ns{2}.nR_S2{s} = [fliplr(nR_S2_err_ns{s}), nR_S2_corr_ns{s}];
        end
    end 
end

%%compute meta-d' difference
mcmc_params.nsamples = 30000;
PKU_s_FIT = fit_meta_d_mcmc_group(metaData_s{1}.nR_S1, metaData_s{1}.nR_S2);
UCL_s_FIT = fit_meta_d_mcmc_group(metaData_s{2}.nR_S1, metaData_s{2}.nR_S2);
PKU_ns_FIT = fit_meta_d_mcmc_group(metaData_ns{1}.nR_S1, metaData_ns{1}.nR_S2);
UCL_ns_FIT = fit_meta_d_mcmc_group(metaData_ns{2}.nR_S1, metaData_ns{2}.nR_S2);

sampleDiff_s = PKU_s_FIT.mcmc.samples.mu_logMratio - UCL_s_FIT.mcmc.samples.mu_logMratio;
hdi_s = calc_HDI(sampleDiff_s(:));
fprintf(('\n Mratio group values = %.2f and %.2f'), exp(PKU_s_FIT.mu_logMratio), exp(UCL_s_FIT.mu_logMratio));
fprintf(['\n Estimated difference in Mratio between groups: ', num2str(exp(PKU_s_FIT.mu_logMratio) - exp(UCL_s_FIT.mu_logMratio))])
fprintf(['\n HDI on difference in log(Mratio): ', num2str(hdi_s) '\n\n'])
color_TOM = {[0.427, 0.298, 0.803], [0.960, 0.325, 0.019]};%http://doc.instantreality.org/tools/color_calculator/

% compute probability of difference
temp_s = sampleDiff_s < 0;
p_theta_s = (sum(temp_s(:) == 1))/30000
subplot(1,2,2)
axis_text = 21
hold all
HDI = calc_HDI(exp(PKU_s_FIT.mcmc.samples.mu_logMratio(:)));
leg1=xline(HDI(1), '--', 'color', color_TOM{1}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{1}, 'linewidth', 2)
histogram(exp(PKU_s_FIT.mcmc.samples.mu_logMratio(:)),500, 'facecolor', color_TOM{1}, 'facealpha', 0.4, 'edgecolor', color_TOM{1}, 'edgealpha', 0.4);
histogram(exp(UCL_s_FIT.mcmc.samples.mu_logMratio(:)),500,'facecolor', color_TOM{2}, 'facealpha', 0.4, 'edgecolor', color_TOM{2}, 'edgealpha', 0.4);
HDI = calc_HDI(exp(UCL_s_FIT.mcmc.samples.mu_logMratio(:)));
leg2 = xline(HDI(1), '--', 'color', color_TOM{2}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{2}, 'linewidth', 2)
xline(0, '-', 'color', 'k', 'linewidth', 1)
ylabel('No. of samples', 'FontSize', axis_text);
set(gca, 'XLim', [1 1.6], 'YLim', [0 250], 'FontSize',22);
xlabel('Posterior distribution', 'FontSize',25);
set(gcf, 'color', 'w')
fakeline1 = plot(2,2000, 'color', color_TOM{1}, 'linewidth', 2);
fakeline2 = plot(2,2000, 'color', color_TOM{2}, 'linewidth', 2);
title('Social condition', 'Fontsize', 22)
% [lgd, handle] = legend([fakeline1, fakeline2], {'PKU', 'UCL'},'location', 'SouthEast');
% linehandle = findobj(handle, 'type', 'line');
% set(linehandle, 'LineWidth',7)
% legend boxoff
% texthandle = findobj(handle, 'type', 'text');
% set(texthandle,'FontSize',15);
% hold off

%%compute meta-d' difference
sampleDiff_ns = PKU_ns_FIT.mcmc.samples.mu_logMratio - UCL_ns_FIT.mcmc.samples.mu_logMratio;
hdi_ns = calc_HDI(sampleDiff_ns(:));
fprintf(('\n Mratio group values = %.2f and %.2f'), exp(PKU_ns_FIT.mu_logMratio), exp(UCL_ns_FIT.mu_logMratio));
fprintf(['\n Estimated difference in Mratio between groups: ', num2str(exp(PKU_ns_FIT.mu_logMratio) - exp(UCL_ns_FIT.mu_logMratio))])
fprintf(['\n HDI on difference in log(Mratio): ', num2str(hdi_ns) '\n\n'])
color_TOM = {[0.427, 0.298, 0.803], [0.960, 0.325, 0.019]};%http://doc.instantreality.org/tools/color_calculator/

% compute probability of difference
temp_ns = sampleDiff_ns < 0;
p_theta_ns = (sum(temp_ns(:) == 1))/30000

hold all
subplot(1,2,1)
HDI = calc_HDI(exp(PKU_ns_FIT.mcmc.samples.mu_logMratio(:)));
leg1=xline(HDI(1), '--', 'color', color_TOM{1}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{1}, 'linewidth', 2)
histogram(exp(PKU_ns_FIT.mcmc.samples.mu_logMratio(:)),500, 'facecolor', color_TOM{1}, 'facealpha', 0.4, 'edgecolor', color_TOM{1}, 'edgealpha', 0.4);
histogram(exp(UCL_ns_FIT.mcmc.samples.mu_logMratio(:)),500,'facecolor', color_TOM{2}, 'facealpha', 0.4, 'edgecolor', color_TOM{2}, 'edgealpha', 0.4);
HDI = calc_HDI(exp(UCL_ns_FIT.mcmc.samples.mu_logMratio(:)));
leg2 = xline(HDI(1), '--', 'color', color_TOM{2}, 'linewidth', 2);
xline(HDI(2), '--', 'color', color_TOM{2}, 'linewidth', 2)
xline(0, '-', 'color', 'k', 'linewidth', 1)
ylabel('No. of samples', 'FontSize', axis_text);
set(gca, 'XLim', [1.3 2], 'YLim', [0 250], 'FontSize',22);
xlabel('Posterior distribution', 'FontSize',25);
set(gcf, 'color', 'w')
fakeline1 = plot(2,2000, 'color', color_TOM{1}, 'linewidth', 2);
fakeline2 = plot(2,2000, 'color', color_TOM{2}, 'linewidth', 2);
title('Perceptual condition', 'FontSize', 22)
[lgd, handle] = legend([fakeline1, fakeline2], {'PKU', 'UCL'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',7)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',15);
hold off

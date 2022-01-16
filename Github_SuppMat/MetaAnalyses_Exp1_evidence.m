%% Meta-analyses
% EVDP 2019 elisa.plas.18@ucl.ac.uk
fs = filesep;

culture = {'PKU', 'UCL'};
sj_mat = {[101:109, 111:115, 117:141],[201:204, 206:227, 229:234, 236:242]};
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
        
        %%levels of evidence
        level = unique(postcoh);
       
        % calculate meta-d' variables
        for e = 1:length(edges) %%for all possible ratings
            for l = 1:3
            nR_S1_corr_1{s}(e) = sum(binned_conf==e & dir==-1 & acc==1 & (postcoh==level(1))'); %how often reported confrating r when dir==left & acc==1
            nR_S1_err_1{s}(e) = sum(binned_conf==e & dir==-1 & acc==0 & (postcoh==level(1))'); %how often reported confrating r when dir==left & acc==0
            nR_S2_corr_1{s}(e) = sum(binned_conf==e & dir==1 & acc==1 & (postcoh==level(1))');%idem for dir==right
            nR_S2_err_1{s}(e) = sum(binned_conf==e & dir==1 & acc==0 & (postcoh==level(1))');
            
            nR_S1_corr_2{s}(e) = sum(binned_conf==e & dir==-1 & acc==1 & (postcoh==level(2))'); %how often reported confrating r when dir==left & acc==1
            nR_S1_err_2{s}(e) = sum(binned_conf==e & dir==-1 & acc==0 & (postcoh==level(2))'); %how often reported confrating r when dir==left & acc==0
            nR_S2_corr_2{s}(e) = sum(binned_conf==e & dir==1 & acc==1 & (postcoh==level(2))');%idem for dir==right
            nR_S2_err_2{s}(e) = sum(binned_conf==e & dir==1 & acc==0 & (postcoh==level(2))');
            
            nR_S1_corr_3{s}(e) = sum(binned_conf==e & dir==-1 & acc==1 & (postcoh==level(3))'); %how often reported confrating r when dir==left & acc==1
            nR_S1_err_3{s}(e) = sum(binned_conf==e & dir==-1 & acc==0 & (postcoh==level(3))'); %how often reported confrating r when dir==left & acc==0
            nR_S2_corr_3{s}(e) = sum(binned_conf==e & dir==1 & acc==1 & (postcoh==level(3))');%idem for dir==right
            nR_S2_err_3{s}(e) = sum(binned_conf==e & dir==1 & acc==0 & (postcoh==level(3))');
            end
        end

            metaData{1}.nR_S1{s} = [fliplr(nR_S1_corr_1{s}), nR_S1_err_1{s}];
            metaData{1}.nR_S2{s} = [fliplr(nR_S2_err_1{s}), nR_S2_corr_1{s}];
            metaData{2}.nR_S1{s} = [fliplr(nR_S1_corr_2{s}), nR_S1_err_2{s}];
            metaData{2}.nR_S2{s} = [fliplr(nR_S2_err_2{s}), nR_S2_corr_2{s}];
            metaData{3}.nR_S1{s} = [fliplr(nR_S1_corr_3{s}), nR_S1_err_3{s}];
            metaData{3}.nR_S2{s} = [fliplr(nR_S2_err_3{s}), nR_S2_corr_3{s}];
    end
end
%%compute meta-d' difference
mcmc_params.nsamples = 30000;
level1_FIT = fit_meta_d_mcmc_group(metaData{1}.nR_S1, metaData{1}.nR_S2);
level2_FIT = fit_meta_d_mcmc_group(metaData{2}.nR_S1, metaData{2}.nR_S2);
level3_FIT = fit_meta_d_mcmc_group(metaData{3}.nR_S1, metaData{3}.nR_S2);

subplot(1,2,2)
hold all
HDI = calc_HDI(exp(level1_FIT.mcmc.samples.mu_logMratio(:)));
leg1=xline(HDI(1), '--', 'color', 'b', 'linewidth', 2);
xline(HDI(2), '--', 'color', 'b', 'linewidth', 2)
HDI = calc_HDI(exp(level2_FIT.mcmc.samples.mu_logMratio(:)))
leg2 = xline(HDI(1), '--', 'color', 'g', 'linewidth', 2);
xline(HDI(2), '--', 'color', 'g', 'linewidth', 2)
HDI = calc_HDI(exp(level3_FIT.mcmc.samples.mu_logMratio(:)));
leg3 = xline(HDI(1), '--', 'color', 'r', 'linewidth', 2);

fprintf(('\n Mratio group values = %.2f, %.2f and %.2f'), exp(level1_FIT.mu_logMratio), exp(level2_FIT.mu_logMratio), exp(level3_FIT.mu_logMratio))

sampleDiff = level1_FIT.mcmc.samples.mu_logMratio - level2_FIT.mcmc.samples.mu_logMratio;
p_theta = (sum(temp(:) == 1))/30000;

histogram(exp(level1_FIT.mcmc.samples.mu_logMratio(:)),500, 'facecolor', 'b', 'facealpha', 0.4, 'edgecolor','b', 'edgealpha', 0.4);
histogram(exp(level2_FIT.mcmc.samples.mu_logMratio(:)),500,'facecolor', 'g', 'facealpha', 0.4, 'edgecolor', 'g', 'edgealpha', 0.4);
histogram(exp(level3_FIT.mcmc.samples.mu_logMratio(:)),500,'facecolor', 'r', 'facealpha', 0.4, 'edgecolor','r', 'edgealpha', 0.4);

xline(0, '-', 'color', 'k', 'linewidth', 1)
ylabel('No. of samples', 'FontSize', 25);
set(gca, 'XLim', [0 4], 'YLim', [0 250], 'FontSize',22);
xlabel('Posterior distribution', 'FontSize',25);
set(gcf, 'color', 'w')
fakeline1 = plot(2,2000, 'color','b', 'linewidth', 2);
fakeline2 = plot(2,2000, 'color', 'g', 'linewidth', 2);
fakeline3 = plot(2,2000, 'color', 'r', 'linewidth', 2);
[lgd, handle] = legend([fakeline1, fakeline2 fakeline3], {'low PDE', 'med. PDE', 'str. PDE'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',7)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',23);
hold off


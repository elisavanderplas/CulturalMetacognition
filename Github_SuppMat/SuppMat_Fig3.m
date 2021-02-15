%% Plots confidence of advisers in PKU and UCL (supplementary material 2.3)
% EVDP 2019 elisa.plas.18@ucl.ac.uk

fs = filesep;
c.corr = {[0.647, 0.835, 0.521],[0.368, 0.713, 0.145],[0.250, 0.415, 0.145]};
c.err = {[0.909, 0.784, 0.772],[0.850, 0.070, 0.070], [0.709, 0.090, 0.031]};

%% discret. 33% quant. and save in three bins
allConf_adviser = cell(1,3);
allAcc_adviser = cell(1,3);
allConf_adviser_corr = cell(1,3);
allConf_adviser_err = cell(1,3);

culture = {'PKU', 'UCL'};
sj_mat = {[403:418, 420:432, 434:435, 437:443, 445:459];[25:76,79]};

for n = 1:length(culture)
    nat = culture{n};
    
    baseDir =  ['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs 'DATA' fs 'EXP2' fs];
    dirData = [baseDir nat '_data' fs nat '_data' fs];
   
    filename = 'Data_sub_';   
    subjects = sj_mat{n};
    
    for s = 1:length(subjects)

        datafile = [filename num2str(subjects(s)) '_2.mat'];
        cd(dirData);
        load(datafile);
        
        precoh = locDATA.dots_coherence';
        postcoh = locDATA.post_coherence';
        task = locDATA.condition;
        acc_adv = locDATA.acc_adv;
        conf_adv = locDATA.conf_adv;
        coherence = unique(precoh);
        for i = 1:3
            postcoh_index(locDATA.post_coherence == coherence(i)) =i;
        end
        
        for i = 1:length(conf_adv) %one dir conf adv
            if conf_adv(i) < 0.5
                unsigned_conf_adv(i) = 1 - conf_adv(i);
            else
                unsigned_conf_adv(i) = conf_adv(i);
            end
        end
        
        unsigned_conf_adv(unsigned_conf_adv == 99.000) = NaN; 
        for post = 1:3
            for acc = 1:2
                allConf_adviser{post} = [allConf_adviser{post}, unsigned_conf_adv(postcoh_index == post & task == 1)];
                allAcc_adviser{post} = [allAcc_adviser{post}, acc_adv(postcoh_index == post & task == 1)];
                allConf_adviser_corr{post} = [allConf_adviser_corr{post}, unsigned_conf_adv(postcoh_index == post & task == 1 & acc_adv== 1)];
                allConf_adviser_err{post} = [allConf_adviser_err{post}, unsigned_conf_adv(postcoh_index == post & task == 1 & acc_adv== 0)];
            end
        end
    end
end

h05 = figure;
set(h05,'units','points','position',[10,10,900,500])
xpos = [0.6, 0.75, 0.9];

for post = 1:3
    subplot(1,2,1)
    hold all; box off
    yyaxis right
    box{post} = errorbar(xpos(post),mean(allAcc_adviser{post}), std(allAcc_adviser{post})./sqrt(length(subjects)),'o-', 'markersize', 10, 'markeredgecolor', 'k', 'markerfacecolor',[0.2,0.2,0.2]*(4-post) ,'linewidth', 2, 'color',[0.2,0.2,0.2]*(4-post));
    ylabel('Choice Accuracy', 'FontSize', 20)
    set(gca, 'ycolor', 'k', 'YLim', [0 1], 'YTick', [0:0.2:1], 'FontSize', 20);
    yyaxis left
    [f, xi]  = ksdensity(allConf_adviser_corr{post});
    corr{post} = plot(xi, f, '-','color', c.corr{post}, 'linewidth', 5);
    ylabel('Confidence Density', 'FontSize', 20)
    xlabel('Confidence Adviser', 'FontSize', 20);
    set(gca,'ycolor', 'k','YLim', [0 10], 'YTick', [0:2:10], 'XLim', [0.4 1.1],'XTick', [0.5:0.1:1], 'FontSize', 20);
    
    subplot(1,2,2)
    hold all; box off
    [f, xi]  = ksdensity(allConf_adviser_err{post});
    err{post} = plot(xi, f, '-','color', c.err{post}, 'linewidth', 5);
    ylabel('Confidence Density', 'FontSize', 20)
    xlabel('Confidence Adviser', 'FontSize', 20);
    set(gca, 'ycolor', 'k', 'YLim', [0 10], 'YTick', [0:2:10], 'XLim', [0.4 1.1],'XTick', [0.5:0.1:1], 'FontSize', 20);
end
    
    [lgd, handle] = legend([box{1} box{2} box{3}], {'weak', 'medium', 'strong'},'location', 'NorthEast', 'FontSize', 20);
    legend boxoff
    linehandle = findobj(handle, 'type', 'line');
    set(linehandle, 'LineWidth',3)
  

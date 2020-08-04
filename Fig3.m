%% Plot first-order performance in the social and non-social condition of Exp 2
fs = filesep;

culture = {'PKU';'UCL'};
sj_mat = {[403:418, 420:432, 434:435, 437:443, 445:459]; [25:76, 79]}; 

h02 = figure; 
set(h02,'units','points','position',[10,10,600,550])

%Colours
c.green_ns = {[0.250, 0.415, 0.145],[1, 1, 1]};
c.green_s = {[0.647, 0.835, 0.521],[1, 1, 1]};

linestyle = {'d', 'o'};

for n = 1:length(culture)
    nat = culture{n};

    baseDir =  ['~' fs 'Dropbox' fs 'CulturalMetacognition_2020' fs];
    dirData = [baseDir 'DATA' fs 'EXP2' fs nat '_data' fs nat '_data' fs]; 

    filename = 'fMRI_pilotData_sub_'; 
    suffix = '';
    subjects = sj_mat{n};
    
    for s = 1:length(subjects)
        
        %% Load data for this subject
        datafile = [filename num2str(subjects(s)) suffix '_2.mat'];
        cd(dirData);
        load(datafile);
        cd(cwd);
        
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
        
       for t = 1:length(agree)
            if conf(t) < 0.5 && transformed_action(t) == -1
                transformed_action_post(t) = 1;
            elseif conf(t) > 0.5 && transformed_action(t) == 1
                transformed_action_post(t) = 1;
            else
                transformed_action_post(t) = -1;
            end
       end
        %how often do sj and adviser disagree?
        agree_pr(n,s) = (sum(task == 1 & agree == 0));      
        % how often disagreement about final confidence judgement?
        agree_post =  transformed_action_post == transformed_action_adv;
        % how often do adviser and participant resolve disagreement?
        switch2adv(n,s) = (sum(task == 1 & agree == 0 & agree_post == 1));
        
        % compute how often adviser and participant resolve disagreement,
        % on trials on which the participant was initially wrong or correct
        for acy = 1:2
            switch2adv_acc(n,s,acy) = (sum(task ==1 & agree == 0 & agree_post == 1 & acc_advt == acy));
            agree_pr_acc(n,s,acy) = (sum(task == 1 & agree == 0 & acc_advt == acy)); 
        end
        mean_conf(n,s) = nanmean(conf);
    end
    mean_acc_social = nanmean(accuracy_social); 
    sem_acc_social = nanstd(accuracy_social)./sqrt(length(subjects));
    mean_acc_nonsocial = nanmean(accuracy_nonsocial);
    sem_acc_nonsocial = nanstd(accuracy_nonsocial)./sqrt(length(subjects));
    
    hold all
    xpos =1:3;
    %%overlay subjects as dots - nonsocial
    for i = 1:length(subjects)
     r = -0.1 + (0.1+0.1).*rand(1,1);
    scat = plot(xpos+r, accuracy_nonsocial(i,:),  linestyle{n}, 'MarkerSize', 5, 'MarkerFaceColor',c.green_ns{n}, 'MarkerEdgeColor', c.green_ns{1}, 'LineWidth', 0.5);
    scat.Color(4) = 0.1; 
    end
    %%overlay subjects as dots - social
    for i = 1:length(subjects)
     r = -0.1 + (0.1+0.1).*rand(1,1);
    scat = plot(xpos+r, accuracy_social(i,:),  linestyle{n}, 'MarkerSize', 5, 'MarkerFaceColor', c.green_s{n}, 'MarkerEdgeColor', c.green_s{1}, 'LineWidth', 0.5);
    scat.Color(4) = 0.1; 
    end
    nso{n} = errorbar(xpos, mean_acc_nonsocial,sem_acc_nonsocial, '-', 'LineWidth', 5, 'Color', c.green_ns{1}); 
    so{n} = errorbar(xpos, mean_acc_social,sem_acc_social, '-', 'LineWidth', 5, 'Color', [0.588, 0.780, 0.458]); 
end
fakeline_PKU = errorbar(xpos, mean_acc_social+2000, sem_acc_social, '-d', 'LineWidth', 5, 'MarkerFaceColor','k', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'Color', 'k'); %for the legend
fakeline_UCL = errorbar(xpos, mean_acc_social+2000, sem_acc_social, '--o', 'LineWidth', 5, 'MarkerFaceColor', 'w', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'Color', 'k'); 
nso{2}.LineStyle = '--';
so{2}.LineStyle = '--';

[lgd, handle] = legend([fakeline_PKU, fakeline_UCL, nso{1}, so{1}], {'PKU', 'UCL', 'perceptual', 'social'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'dot');
set(linehandle, 'LineWidth',8)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);

set(gca, 'YLim', [0.5 1], 'XLim', [0.5  3.5], 'XTick', xpos, 'YTick', .5:0.1:1, 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', 20);
ylabel('Choice accuracy', 'FontSize', 30);
xlabel('Pre-decision evidence', 'FontSize', 30);

%line 230 of manuscript, basic confidence differences
nanmean(mean_conf(1,:));%mean/sem, 1 = PKU, 2 = UCL
nanstd(mean_conf(1,:)./sqrt(length(mean_conf(1,:))));
[H,P,CI,STATS] = ttest2(mean_conf(1,:), mean_conf(2,:)); 

%line 254-264, basic advice-taking differences 
p_switch = switch2adv./agree_pr;%proportion of trials on which participant changed their mind towards advised dir
p_switch_acc = switch2adv_acc./agree_pr_acc; %same as p_switch, but separates out initially erroneous vs correct trials 

disp(['average compliance with adviser in PKU : ' num2str(mean(p_switch(1,:)))]);
disp(['average compliance with adviser in UCL : ' num2str(mean(p_switch(2,:)))]);

[H,P,CI,STATS] = ttest2(p_switch(1,:), p_switch(2,:)); 
disp(['p-level difference : ' num2str(P)]);

disp(['average compliance after initial error in PKU : ' num2str(mean(p_switch_acc(1,:,2)))]);
disp(['average compliance after initial error in UCL : ' num2str(mean(p_switch_acc(2,:,2)))]);
[H,P,CI,STATS] = ttest2(p_switch_acc(1,:,2), p_switch_acc(2,:,2)); 
disp(['compliance with correct advisers P-value: ' num2str(P)]);

disp(['average compliance after having been initially correct in PKU : ' num2str(mean(p_switch_acc(1,:,1)))]);
disp(['average compliance after having been initially correct in UCL : ' num2str(mean(p_switch_acc(2,:,1)))]);
[H,P,CI,STATS] = ttest2(p_switch_acc(1,:,1), p_switch_acc(2,:,1)); 
disp(['compliance with correct advisers P-value: ' num2str(P)]);
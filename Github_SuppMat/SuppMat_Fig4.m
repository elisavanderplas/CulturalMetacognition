% elisa.plas.18@ucl.ac.uk

%% v-plot also for social and non-social condition Exp 2

clear all; close all
fs = filesep;

c.corr =  {[0.082, 0.615, 0.835],[0.082, 0.615, 0.835]};
c.err= {[0.835, 0.250, 0.082], [0.835, 0.250, 0.082]};

h0 = figure;
set(h0,'units','points','position',[10,10,1200,2400])

culture = {'PKU';'UCL'};
sj_mat = {[403:418, 420:432, 434:435, 437:443, 445:459]; [25:76, 79]}; 

for n = 1:length(culture)
    nat = culture{n};
    subjects = sj_mat{n};
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
        postcoh_index = [];
        
        precoh = locDATA.dots_coherence';
        postcoh = locDATA.post_coherence';
        dir = locDATA.dots_direction/360;
        dir(dir==0.5) = -1;
        action = locDATA.button_response - 1;
        conf = locDATA.mouse_response;
        logrt = log(locDATA.reaction_time_button);
        conf_rt = log(locDATA.reaction_time_mouse);
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        conf_adv = locDATA.conf_adv; 
        task = locDATA.condition; %1 = social, 0 = nonsocial
        
        %% group conf adviser into three bins

            conf_adv(conf_adv == 99.000) = NaN;
  
    binned_conf_adv = discretize(conf_adv,3);
  
        %% split id vars into social & nonsocial
        
        % Add indicator variables for pre/post confidence
        coherence = unique(precoh);
        
        for i = 1:3
            precoh_index(locDATA.dots_coherence==coherence(i))=i;
        end
        for i = 1:3
            postcoh_index(locDATA.post_coherence==coherence(i))=i;
        end
        
        conf(conf < 0) = 0;
        
             
          %% Confidence/performance for different tasks (social/nonsocial)
          for post = 1:3 
              for pre = 1:3
            conf_cor_social{n}(s,pre, post) = nanmean(conf(binned_conf_adv == post & acc == 1 & task == 1 & precoh_index == pre));
            conf_cor_nonsocial{n}(s, pre, post) = nanmean(conf(postcoh_index == post & acc == 1 & task == 0 & precoh_index == pre));
            conf_err_nonsocial{n}(s,pre, post) = nanmean(conf(postcoh_index == post & acc == 0 & task == 0 & precoh_index == pre));     
            conf_err_social{n}(s,pre, post) = nanmean(conf(binned_conf_adv == post & acc == 0 & task == 1 & precoh_index == pre)); 
              end
          end

    end
    
end

%% plot non-social condition
marker = {'-d', '--o'};
ms = 11; 
axis_text = 17; 
axis_nr = 15; 
for pre = 1:3
    pre_level = {'Weak', 'Medium', 'Strong'};
    subplot(2,3,pre);
    
    for n = 1:2
        hold all; box off; grid off; 
        spacing = 0.5;
        xpos = [1 2 3];
        
        len_subj{n} = length(conf_cor_nonsocial{n}(:,pre,:));
        mean_conf_cor{n} = nanmean(conf_cor_nonsocial{n}(:,pre,:));
        sem_conf_cor{n} = nanstd(conf_cor_nonsocial{n}(:,pre,:))./sqrt(len_subj{n});
        mean_conf_err{n} = nanmean(conf_err_nonsocial{n}(:,pre,:));
        sem_conf_err{n} = nanstd(conf_err_nonsocial{n}(:,pre,:))./sqrt(len_subj{n});
        
        % Corrects
        for group = 1
            x = xpos((group*3)-2:group*3);
            mu = mean_conf_cor{n}((group*3)-2:group*3);
            se = sem_conf_cor{n}((group*3)-2:group*3);
            for p = 1:length(x)-1
                pX = [x(p) x(p+1) x(p+1) x(p)];
                pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
                h = fill(pX,pY,c.corr{n});
                set(h,'edgecolor',c.corr{n},'facealpha',.1,'LineStyle','none');
            end
            hLine(1,n) = plot(x, squeeze(mu), marker{n}, 'markersize', ms, 'Color', c.corr{n},'markerfacecolor', 'w', 'MarkerEdgeColor', c.corr{n},'LineWidth', 3);
        end
        
        % Errors
        for group = 1
            x = xpos((group*3)-2:group*3);
            mu = mean_conf_err{n}((group*3)-2:group*3);
            se = sem_conf_err{n}((group*3)-2:group*3);
            for p = 1:length(x)-1
                pX = [x(p) x(p+1) x(p+1) x(p)];
                pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
                h = fill(pX,pY, c.err{n});
                set(h,'edgecolor',c.err{n},'facealpha',.1,'LineStyle','none');
            end
            hLine(2,n) = plot(x, squeeze(mu),marker{n}, 'markersize', ms, 'Color', c.err{n},'markerfacecolor', 'w','MarkerEdgeColor', c.err{n}, 'LineWidth', 3);
        end
    end
    %get lines to display in the legend
    fakeline_PKU = plot(x+300,squeeze(mu)+ 200,'d-', 'markersize', ms, 'Color', 'k', 'markerfacecolor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    fakeline_UCL = plot(x+300, squeeze(mu)+200,'o--', 'markersize', ms, 'Color', 'k', 'markerfacecolor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    fakeline_corr = plot(x+300,squeeze(mu)+ 200,'-', 'markersize', ms, 'Color', c.corr{2}, 'LineWidth', 3);
    fakeline_err = plot(x+300, squeeze(mu)+200,'-', 'markersize', ms, 'Color',c.err{2},  'LineWidth', 3);

    set(gca, 'YLim', [0 1], 'XLim', [0.65 3.15], 'XTick', xpos, 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', axis_nr);
    ylabel('Confidence', 'FontSize', axis_text);
   
    text(xpos(2)-0.5, 1.03, pre_level{pre}, 'FontSize', axis_nr)
    text(xpos(2)-1.4, -0.12,'Post-decision evidence', 'FontSize', axis_text)
end

text(xpos(2)-5.5, 1.09,'Pre-decision evidence', 'FontSize', axis_text) 

%% plot social condition
marker = {'--s', '-o'};
ms = 11; 
axis_text = 17; 
axis_nr = 15; 
for pre = 1:3
    pre_level = {'Weak', 'Medium', 'Strong'};
    subplot(2,3,3+pre);
    
    for n = 1:2
        hold all; box off; grid off; 
        spacing = 0.5;
        xpos = [1 2 3];
        
        len_subj{n} = length(conf_cor_social{n}(:,pre,:));
        mean_conf_cor{n} = nanmean(conf_cor_social{n}(:,pre,:));
        sem_conf_cor{n} = nanstd(conf_cor_social{n}(:,pre,:))./sqrt(len_subj{n});
        mean_conf_err{n} = nanmean(conf_err_social{n}(:,pre,:));
        sem_conf_err{n} = nanstd(conf_err_social{n}(:,pre,:))./sqrt(len_subj{n});
        
        % Corrects
        for group = 1
            x = xpos((group*3)-2:group*3);
            mu = mean_conf_cor{n}((group*3)-2:group*3);
            se = sem_conf_cor{n}((group*3)-2:group*3);
            for p = 1:length(x)-1
                pX = [x(p) x(p+1) x(p+1) x(p)];
                pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
                h = fill(pX,pY,c.corr{n});
                set(h,'edgecolor',c.corr{n},'facealpha',.1,'LineStyle','none');
            end
            hLine(1,n) = plot(x, squeeze(mu), marker{n}, 'markersize', ms, 'Color', c.corr{n},'markerfacecolor', 'w', 'MarkerEdgeColor', c.corr{n},'LineWidth', 3);
        end
        
        % Errors
        for group = 1
            x = xpos((group*3)-2:group*3);
            mu = mean_conf_err{n}((group*3)-2:group*3);
            se = sem_conf_err{n}((group*3)-2:group*3);
            for p = 1:length(x)-1
                pX = [x(p) x(p+1) x(p+1) x(p)];
                pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
                h = fill(pX,pY, c.err{n});
                set(h,'edgecolor',c.err{n},'facealpha',.1,'LineStyle','none');
            end
            hLine(2,n) = plot(x, squeeze(mu),marker{n}, 'markersize', ms, 'Color', c.err{n},'markerfacecolor', 'w','MarkerEdgeColor', c.err{n}, 'LineWidth', 3);
        end
    end
    %get lines to display in the legend
    fakeline_PKU = plot(x+300,squeeze(mu)+ 200,'d-', 'markersize', ms, 'Color', 'k', 'markerfacecolor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    fakeline_UCL = plot(x+300, squeeze(mu)+200,'o--', 'markersize', ms, 'Color', 'k', 'markerfacecolor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    fakeline_corr = plot(x+300,squeeze(mu)+ 200,'-', 'markersize', ms, 'Color', c.corr{2}, 'LineWidth', 3);
    fakeline_err = plot(x+300, squeeze(mu)+200,'-', 'markersize', ms, 'Color',c.err{2},  'LineWidth', 3);

    set(gca, 'YLim', [0 1], 'XLim', [0.65 3.15], 'XTick', xpos, 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', axis_nr);
    ylabel('Confidence', 'FontSize', axis_text);
   
    text(xpos(2)-0.5, 1.03, pre_level{pre}, 'FontSize', axis_nr)
    text(xpos(2)-1.4, -0.12,'Post-decision evidence', 'FontSize', axis_text)
end

text(xpos(2)-5.5, 1.09,'Pre-decision evidence', 'FontSize', axis_text) 


%%set legend
[lgd, handle] = legend([fakeline_PKU, fakeline_UCL,  fakeline_corr, fakeline_err], {'PKU', 'UCL', 'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'dot');
set(linehandle, 'LineWidth',6)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',axis_nr);
set(gcf, 'color', 'w')

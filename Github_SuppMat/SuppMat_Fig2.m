% Steve Fleming stephen.fleming@ucl.ac.uk 2017
% Adapted by Elisa van der Plas elisavanderplas@gmail.com

clear all; close all
fs = filesep;

c.corr =  {[0.6, 0.6, 0.6],[0.082, 0.615, 0.835],[0.082, 0.615, 0.835] };
c.err= {[0.6, 0.6, 0.6],[0.835, 0.250, 0.082], [0.835, 0.250, 0.082]};

h0 = figure;
set(h0,'units','points','position',[10,10,1200,2400])

culture = {'NYU';'UCL';'PKU';};
sj_mat = {[12:28,30:37]; [201:204, 206:227, 229:234, 236:242];[101:109, 111:115, 117:141]};

for n = 1:length(culture)
    nat = culture{n};
    subjects = sj_mat{n};
    baseDir = '~/Dropbox/PKU_collaboration/Github/DATA/EXP1/';
    dirData = [baseDir  nat '_data/' nat '_data/' ];
    cwd = pwd;
    filename = 'fMRI_pilotData_sub_';
    suffix = '';
    
    for s = 1:length(subjects)
        
        %% Load data for this subject
        datafile = [filename num2str(subjects(s)) suffix '_2.mat'];
        cd(dirData);
        load(datafile);
        cd(cwd);
        
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
        
        % Add indicator variables for pre/post confidence
        coherence = unique(precoh);
        
        for i = 1:3
            precoh_index(locDATA.dots_coherence==coherence(i))=i;
        end
        for i = 1:3
            postcoh_index(locDATA.post_coherence==coherence(i))=i;
        end
        
        conf(conf < 0) = 0;
        
        %% Confidence in 3x3 conditions
        for post = 1:3
            for pre = 1:3
                conf_cor{n}(s,pre, post) = nanmean(conf(precoh_index == pre & postcoh_index == post & acc == 1));
                conf_err{n}(s,pre, post) = nanmean(conf(precoh_index == pre & postcoh_index == post & acc == 0));
            end
        end
    end
    
    dirBeta = [dirData nat '_betas' fs];   
    filename = 'regression_betas_';
    %% Load data for this dataset
    datafile = [filename nat '.csv']; %betas
    cd(dirBeta);
    dat{n}= readtable(datafile);
    
    for cor = 1:2
        correct = {'corr_'; 'err_'}; 
        
    suffix = correct{cor}; % use for fMRI data out-of-scanner
    %% Load data for this dataset
    datafile = [filename suffix nat '.csv']; %betas
    dat_acc{n, cor} = readtable(datafile);
    end
    % specify betas
    b_acc(n)= dat{n}{2,2};
    b_pre(n)= dat{n}{3,2};
    b_pos(n)= dat{n}{4,2};
    b_accXpre(n)= dat{n}{6,2};
    b_accXpos(n)=dat{n}{7,2};
    b_preXpos(n) =dat{n}{8,2};
    b_accXpreXpos(n)= dat{n}{9,2};
end

%% 3 x 3 confidence
marker = {'--s', '--o', '-d'};
ms = 11; 
axis_text = 17; 
axis_nr = 15; 
for pre = 1:3
    pre_level = {'  Weak', 'Medium', ' Strong'};
    subplot(2,4,pre);
    
    for n = 1:3
        hold all; box off; grid off; 
        spacing = 0.5;
        xpos = [1 2 3];
        
        len_subj{n} = length(conf_cor{n}(:,pre,:));
        mean_conf_cor{n} = nanmean(conf_cor{n}(:,pre,:));
        sem_conf_cor{n} = nanstd(conf_cor{n}(:,pre,:))./sqrt(len_subj{n});
        mean_conf_err{n} = nanmean(conf_err{n}(:,pre,:));
        sem_conf_err{n} = nanstd(conf_err{n}(:,pre,:))./sqrt(len_subj{n});
        
        % Corrects - model
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
        
        % Corrects - subjects
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
    fakeline_NYU = plot(x+300, squeeze(mu)+200,'s--', 'markersize', ms, 'Color', [0.6, 0.6, 0.6], 'markerfacecolor', 'w', 'MarkerEdgeColor', [0.6, 0.6, 0.6], 'LineWidth', 3);
    fakeline_corr = plot(x+300,squeeze(mu)+ 200,'-', 'markersize', ms, 'Color', c.corr{3}, 'LineWidth', 3);
    fakeline_err = plot(x+300, squeeze(mu)+200,'-', 'markersize', ms, 'Color',c.err{3},  'LineWidth', 3);

    set(gca, 'YLim', [0 1], 'XLim', [0.65 3.15], 'XTick', xpos, 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', axis_nr);
    ylabel('Confidence', 'FontSize', axis_text);
   
    text(xpos(2)-0.5, 1.03, pre_level{pre}, 'FontSize', axis_nr)
    text(xpos(2)-1.4, -0.12,'Post-decision evidence', 'FontSize', axis_text)
end

text(xpos(2)-5.5, 1.09,' Pre-decision evidence', 'FontSize', axis_text) 

%%set legend
[lgd, handle] = legend([fakeline_PKU, fakeline_UCL, fakeline_NYU, fakeline_corr, fakeline_err], {'PKU', 'UCL', 'NYU', 'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'dot');
set(linehandle, 'LineWidth',6)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',axis_nr);


% specify variables
acc_v= [-1 1];
pre_v= [-0.5 0 0.5];
pos_v= [-0.5 0 0.5];
con = [];% perform posterior predictives

for n = 1:3
for i_acc= 1:2
    for i_pre= 1:3
        for i_pos= 1:3
           con(n, i_acc,i_pre,i_pos)= acc_v(i_acc)* b_acc(n) + ...
                                   pre_v(i_pre)* b_pre(n) + ...
                                   pos_v(i_pos)* b_pos(n) + ...
                                   acc_v(i_acc)* pre_v(i_pre)* b_accXpre(n) + ...
                                   acc_v(i_acc)* pos_v(i_pos)* b_accXpos(n) + ...
                                   pre_v(i_pre)* pos_v(i_pos)* b_preXpos(n) + ...
                                   acc_v(i_acc)* pre_v(i_pre)* pos_v(i_pos)* b_accXpreXpos(n);
        end
    end
end 

for i = 1:3
% plot posterior predictives
subplot(2,4,4+i)
pre_level = {'Weak', 'Medium', 'Strong'}; 
hold on; grid off; box off; 

pa(n)= plot(squeeze(con(n,1,i,:)), marker{n},'markersize', ms, 'Color', c.err{n}, 'markerfacecolor', 'w', 'MarkerEdgeColor', c.err{n}, 'LineWidth', 3);
pb(n) = plot(squeeze(con(n,2,i,:)), marker{n},'markersize', ms, 'Color', c.corr{n}, 'markerfacecolor', 'w', 'MarkerEdgeColor', c.corr{n}, 'LineWidth', 3);

    
set(gca, 'YLim', [-.6 0.4], 'XLim', [0.65 3.15], 'XTick', xpos, 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', axis_nr);
ylabel([{'Post-decision evidence'};{'impact on confidence (a.u)'}], 'FontSize', axis_text);
    
text(xpos(2)-0.5, 0.45, pre_level{i}, 'FontSize', axis_nr)
text(xpos(2)-1.5, -1.52,'Post-decision evidence', 'FontSize', axis_text)
    
end
end
 
text(xpos(2)-5.5, 0.53,' Pre-decision evidence', 'FontSize', axis_text) 

subplot(2,4,4);
hold all; grid off; box off
s = 2; %predecision evidence
xpos = 1; 
aa = plot(xpos,  [dat_acc{3,1}{s,2}  ], 'd', 'MarkerSize',ms,  'MarkerFaceColor', c.corr{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{3});
ab = plot(xpos,  [dat_acc{3,2}{s,2}  ] , 'd', 'MarkerSize', ms,'MarkerFaceColor',c.err{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{3});
xpos = 2;
aa = plot(xpos,  [dat_acc{2,1}{s,2} ] , 'o', 'MarkerSize',ms,  'MarkerFaceColor', c.corr{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{3});
ab = plot(xpos,  [dat_acc{2,2}{s,2}] , 'o', 'MarkerSize', ms,'MarkerFaceColor',c.err{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{3});
xpos = 3; 
aa = plot(xpos,  [dat_acc{1,1}{s,2} ] , 's', 'MarkerSize',ms,  'MarkerFaceColor', c.corr{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{3});
ab = plot(xpos,  [dat_acc{1,2}{s,2}] , 's', 'MarkerSize', ms,'MarkerFaceColor',c.err{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{3});
xpos = 1:3;
errorbar(xpos,  [dat_acc{3,1}{s,2} dat_acc{2,1}{s,2} dat_acc{1,1}{s,2}] , [dat_acc{3,1}{s+5,2}  dat_acc{2,1}{s+5,2} dat_acc{1,1}{s+5,2}] , '.', 'Color', c.corr{3}*0.2, 'LineWidth', 2);
errorbar(xpos, [dat_acc{3,2}{s,2} dat_acc{2,2}{s,2} dat_acc{1,2}{s,2}] , [dat_acc{3,2}{s+5,2}  dat_acc{2,2}{s+5,2} dat_acc{1,2}{s+5,2}], '.', 'Color', c.err{3}*0.2, 'LineWidth', 2);

hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7); %zeroline
set(gca, 'XLim', [0 4], 'XTick', xpos, 'XTickLabel',['PKU  '; '  UCL' ; '  NYU'], 'YLim', [-0.5 0.2],'YTick', [-0.5:0.1:0.2], 'FontSize',axis_nr);
ylabel([{'Pre-decision evidence'};{'impact on confidence (a.u.)'}], 'FontSize', axis_text);

subplot(2,4,8);
hold all;grid off; box off
s = 3; %post decision evidence
xpos = 1; 
aa = plot(xpos,  [dat_acc{3,1}{s,2}  ], 'd', 'MarkerSize',ms+2,  'MarkerFaceColor', c.corr{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{3});
ab = plot(xpos,  [dat_acc{3,2}{s,2}  ] , 'd', 'MarkerSize', ms+2,'MarkerFaceColor',c.err{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{3});
xpos = 2;
aa = plot(xpos,  [dat_acc{2,1}{s,2} ] , 'o', 'MarkerSize',ms+2,  'MarkerFaceColor', c.corr{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{3});
ab = plot(xpos,  [dat_acc{2,2}{s,2}] , 'o', 'MarkerSize', ms+2,'MarkerFaceColor',c.err{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{3});
xpos = 3; 
aa = plot(xpos,  [dat_acc{1,1}{s,2} ] , 's', 'MarkerSize',ms+2,  'MarkerFaceColor', c.corr{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{3});
ab = plot(xpos,  [dat_acc{1,2}{s,2}] , 's', 'MarkerSize', ms+2,'MarkerFaceColor',c.err{3}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{3});
xpos = 1:3;
errorbar(xpos,  [dat_acc{3,1}{s,2} dat_acc{2,1}{s,2} dat_acc{1,1}{s,2}] , [dat_acc{3,1}{s+5,2}  dat_acc{2,1}{s+5,2} dat_acc{1,1}{s+5,2}] , '.', 'Color', c.corr{3}*0.3, 'LineWidth', 2);
errorbar(xpos, [dat_acc{3,2}{s,2} dat_acc{2,2}{s,2} dat_acc{1,2}{s,2}] , [dat_acc{3,2}{s+5,2}  dat_acc{2,2}{s+5,2} dat_acc{1,2}{s+5,2}], '.', 'Color', c.err{3}*0.3, 'LineWidth', 2);

hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7); %zeroline

set(gca, 'XLim', [0 4], 'XTick', xpos, 'XTickLabel',['PKU  '; '  UCL' ; '  NYU'], 'YLim', [-0.5 0.2],'YTick', [-0.5:0.1:0.2], 'FontSize',axis_nr);
ylabel([{'Post-decision evidence'};{'impact on confidence (a.u.)'}], 'FontSize', axis_text);


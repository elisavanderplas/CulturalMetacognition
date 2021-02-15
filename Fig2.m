%% Plot Fig2 of PKU-UCL paper
% EVDP 2019 elisa.plas.18@ucl.ac.uk
fs = filesep;

culture = {'PKU', 'UCL'};
sj_mat = {[101:109, 111:115, 117:141],[201:204, 206:227, 229:234, 236:242]};

h01 = figure;
set(h01,'units','points','position',[10,10,1400,600])
c.corr =  [0.082, 0.615, 0.835];
c.err = [0.835, 0.250, 0.082];

for n = 1:length(culture)
    nat = culture{n};
    baseDir =  ['~' fs 'Dropbox' fs 'CulturalMetacognition-master' fs];
    dirData = [baseDir 'DATA' fs 'EXP1' fs nat '_data' fs nat '_data' fs];
    
    filename = 'regression_betas_'; %made in R with 'Regression_EXP1.r'
    
    for acc = 1:2
        correct = {'corr_'; 'err_'};
        suffix = correct{acc};
        
        %% Load betas
        datafile = [filename suffix nat '.csv'];
        cd([dirData fs nat '_betas' fs]);
        dat{n, acc} = readtable(datafile);
    end
    
    filename = 'Data_sub_';
    subjects = sj_mat{n};
    
    for s = 1:length(subjects)
        %% Load data for this subject
        datafile = [filename num2str(subjects(s)) '_2.mat'];
        cd(dirData);
        load(datafile);
        
        precoh = locDATA.dots_coherence';
        postcoh = locDATA.post_coherence';
        dir = locDATA.dots_direction/360;
        dir(dir==0.5) = -1;
        action = locDATA.button_response - 1;
        conf = locDATA.mouse_response;
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        coherence = unique(precoh);
        
        %index pre/post-decision evidence levels
        for i = 1:3
            precoh_index(locDATA.dots_coherence==coherence(i))=i;
            postcoh_index(locDATA.post_coherence==coherence(i))=i;
        end
        
        % load confidence across 3x3x2 conditions
        for post = 1:3
            for pre = 1:3
                conf_cor{n}(s,pre, post) = nanmean(conf(precoh_index == pre & postcoh_index == post & acc == 1));
                conf_err{n}(s,pre, post) = nanmean(conf(precoh_index == pre & postcoh_index == post & acc == 0));
            end
        end
    end
end

marker = {'-d', '--o'};
ms = 20;
axis_text = 26;
axis_nr = 18;

for pre = 1:3
    pre_level = {'  Weak', 'Medium', ' Strong'};
 
    subplot(1,4,pre);
    hold all; box off; grid off 
    for n = 1:2
        x = 1:3;
        len_subj = length(conf_cor{n}(:,pre,:));
        
        % Corrects
        mu = nanmean(conf_cor{n}(:,pre,:));
        se = nanstd(conf_cor{n}(:,pre,:))./sqrt(len_subj);
        for p = 1:length(x)-1
            pX = [x(p) x(p+1) x(p+1) x(p)];
            pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
            h = fill(pX,pY,c.corr);
            set(h,'edgecolor',c.corr,'facealpha',.1,'LineStyle','none');
        end
        plot(x, squeeze(mu),marker{n}, 'markersize', ms-6, 'Color', c.corr, 'markerfacecolor', 'w', 'MarkerEdgeColor', c.corr, 'LineWidth', 4);
        
        %Error
        mu = nanmean(conf_err{n}(:,pre,:));
        se = nanstd(conf_err{n}(:,pre,:))./sqrt(len_subj);
        for p = 1:length(x)-1
            pX = [x(p) x(p+1) x(p+1) x(p)];
            pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
            h = fill(pX,pY, c.err);
            set(h,'edgecolor',c.err,'facealpha',.1,'LineStyle','none');
        end
        plot(x, squeeze(mu),marker{n}, 'markersize', ms-6, 'Color', c.err, 'markerfacecolor', 'w', 'MarkerEdgeColor', c.err, 'LineWidth', 4);
    end
    
    %get lines to display in the legend
    fakeline_PKU = plot(x+300,squeeze(mu)+200,'d-' , 'markersize', ms, 'Color', 'k', 'markerfacecolor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    fakeline_UCL = plot(x+300,squeeze(mu)+200,'o--', 'markersize', ms, 'Color', 'k', 'markerfacecolor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    fakeline_err = plot(x+300,squeeze(mu)+200,'-'  , 'markersize', ms, 'Color', c.err, 'LineWidth', 3);
    fakeline_corr= plot(x+300,squeeze(mu)+200,'-'  , 'markersize', ms, 'Color',c.corr,  'LineWidth', 3);
    
    set(gca, 'YLim', [0 1], 'XLim', [0.65 3.15], 'XTick', 1:3, 'XTickLabel', {'Weak', 'Medium', 'Strong'}, 'FontSize', axis_nr);
    ylabel('Confidence', 'FontSize', axis_text);
    
    text(1.5, 1.03, pre_level{pre}, 'FontSize', axis_nr)
    text(0.2, -0.09,'Post-decision evidence', 'FontSize', axis_text-2)
end

text(-4.2, 1.07,' Pre-decision evidence', 'FontSize', axis_text)

[lgd, handle] = legend([fakeline_PKU, fakeline_UCL, fakeline_corr, fakeline_err], {'PKU', 'UCL', 'correct', 'error'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'dot');
set(linehandle, 'LineWidth',7)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);

%% Second Fig
subplot(1,4,4)

hold all; grid off; box off
plot(1,  [dat{1,1}{3,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', c.corr,'Linewidth', 3, 'MarkerEdgeColor', c.corr);
plot(1,  [dat{1,2}{3,2}] , 'd', 'MarkerSize', ms, 'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
plot(2,  [dat{2,1}{3,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', c.corr,'Linewidth', 3, 'MarkerEdgeColor', c.corr);
plot(2,  [dat{2,2}{3,2}] , 'o', 'MarkerSize', ms, 'MarkerFaceColor', c.err, 'Linewidth', 3, 'MarkerEdgeColor', c.err);
errorbar(1:2, [dat{1,1}{3,2} dat{2,1}{3,2}], [dat{1,1}{8,2} dat{2,1}{8,2}], '.', 'Color', c.corr*0.2, 'LineWidth', 2);
errorbar(1:2, [dat{1,2}{3,2} dat{2,2}{3,2}], [dat{1,2}{8,2} dat{2,2}{8,2}], '.', 'Color', c.err*0.2 , 'LineWidth', 2);

hline1 = line([0 22], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7); %zeroline
set(gca, 'XLim', [0 3], 'XTick', 1:2, 'XTickLabel',['PKU  '; '  UCL'], 'YLim', [-0.5 0.2],'YTick', -0.5:0.1:0.2, 'FontSize',axis_nr);
ylabel([{'Post-decision evidence'};{'impact on confidence (a.u.)'}], 'FontSize', axis_text);


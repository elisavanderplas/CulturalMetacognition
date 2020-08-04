%% Plot Fig1 of PKU-UCL paper
% EVDP 2019 elisa.plas.18@ucl.ac.uk
% Adapted from Steve Fleming 2016 stephen.fleming@ucl.ac.uk 

close all
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
    baseDir =  ['~' fs 'Dropbox' fs 'PKU_collaboration' fs 'Github' fs];
    dirData = [baseDir 'DATA' fs 'EXP1' fs nat '_data' fs nat '_data' fs];
       
    filename = 'fMRI_pilotData_sub_';    
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
        action = locDATA.button_response - 1;
        conf = locDATA.mouse_response;
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        coherence = unique(precoh);
        
        %% extra EvdP 20 July: 
        sigma(n,s) = var(conf); 
        payment(n,s) = locDATA.QSR_score(end); 
       
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


%% extra EvdP 20 July: 
 [H,P,CI,STATS] = ttest2(sigma(1,:), sigma(2,:)); 
[H,P,CI] = ttest2(payment(1,:),payment(2,:)); 
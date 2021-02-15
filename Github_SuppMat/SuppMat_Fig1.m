% Plot some characteristics of the data from all three cultures next to each other
% EVDP 2019 elisa.plas.18@ucl.ac.uk

clear all;close all
fs = filesep;

culture = {'PKU'; 'UCL'; 'NYU'};
sj_mat = {[101:109, 111:115, 117:141];[201:204, 206:227, 229:234, 236:242];[12:28,30:37]};

s01 = figure;
set(s01,'units','points','position',[10,10,1000,800])

%set path to PsychFit function
addpath('~/Dropbox/CulturalMetacognition-master/tools/psychFit-master')

%Colours
c.corr =  {[0.019, 0.211, 0.501],[0.4, 0.650, 0.976], [0 0.45 0.74]};
c.err = {[0.372, 0.164, 0.007], [0.976, 0.580, 0.4], [0.768, 0.337, 0.011]};

avrg_mean_coh = cell(1,3);
acc1 = repmat(NaN,length(sj_mat), length(sj_mat{1}));
for n = 1:length(culture)
    nat = culture{n};

    baseDir = '~/Dropbox/CulturalMetacognition-master/DATA/EXP1/';
    dirData = [baseDir nat '_data/' nat '_data/' ];           

    cwd = pwd;
    accuracy{n} = [];  
    coherence{n} = []; 
    
    filename = 'Data_sub_';   
    subjects = sj_mat{n};
    
    %% Sort subject data
    allConf_cor = cell(1,9);   % stores aggregate confidence data, 9 x (Ntrials*Nsubjects) matrix; slowest changing factor is post
    allConf_err = cell(1,9);
 
    for s = 1:length(subjects)
        
        %% Load calibration data for this subject
        datafile = [filename num2str(subjects(s)) '_1.mat'];
        cd(dirData);
        load(datafile);
        
        ntrials_calib = length(locDATA.dots_direction); 
        
        acc1(n,s) = mean(locDATA.accuracy); 
        
        percCoherence = [] ;
         % Get the Psychometric function for the calibration
         for i = 1:ntrials_calib
            if locDATA.dots_direction(i) ==180
                percCoherence(i) = locDATA.dots_coherence(i)*-1;
            else
                percCoherence(i) = locDATA.dots_coherence(i);
            end
        end
        
        interval1 = [-1, -0.48, -0.24, -0.12, -0.08, -0.03, 0.03, 0.08, 0.12, 0.24, 0.48, 1];
        
        for j = 1:12
            Pright1(s,j) = sum(locDATA.button_response ==2 & percCoherence == interval1(j));
        end
        
        %get nr of trials per interval-condition
        N = repmat((ntrials_calib/length(interval1)), 1, length(interval1));
        pArray = [0 0.5];
        
        %load the model
        fitparams1 = psychFit(interval1, Pright1(s,:), N, pArray, 'normal');
        base1 = linspace(min(interval1), max(interval1), ntrials_calib/length(interval1));
        %store in 3D matrix w/ 1st dimension=sj
        pred1(s,:) = cumNormPred(base1, fitparams1(1), fitparams1(2));
         
        datafile = [filename num2str(subjects(s)) '_2.mat'];
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
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        coh = unique(precoh);
        ntrials = length(precoh);
        
        %save for each country the coherence and accuracy levels 
        coherence{n} = [coherence{n}; coh'];
        accuracy{n} = [accuracy{n}; acc']; 
        acc2(n,s) = mean(acc); 
        RT_action(n,s) = mean(locDATA.reaction_time_mouse);
        RT_conf(n,s) = mean(locDATA.reaction_time_button); 
        
        % Add indicator variables for pre/post confidence
        for i = 1:3
            precoh_index(precoh == coh(i)) =i; 
        end
        
        for i = 1:ntrials
            if locDATA.dots_direction(i) ==180
                percCoherence(i) = precoh_index(i)*-1;
            else
                percCoherence(i) = precoh_index(i); 
            end
        end
        
        interval2 = [-3 -2 -1 1 2 3]; 
        for j = 1:length(interval2)
            Pright2(s, j) = sum(locDATA.button_response==2 & percCoherence == interval2(j)); 
        end
        
        N = repmat(ntrials/length(interval2),1, length(interval2)); 
        pArray = [0 0.5]; 
        
        %load the model
        fitparams = psychFit(interval2, Pright2(s,:), N, pArray, 'normal'); 
        base2 = linspace(min(interval2), max(interval2), ntrials/length(interval2)); 
        %store in 3D matrix w/ 1st dimension=sj
        pred2(s,:,:) = cumNormPred(base2, fitparams(1), fitparams(2));   
             
    end
    mean_coh{n} = nanmean(coherence{n});
    for i = 1:length(coherence{n})
    avrg_mean_coh{n} = [avrg_mean_coh{n} ; nanmean(coherence{n}(i,:))];
    sem_coh{n} = nanstd(coherence{n})./sqrt(length(coherence{n}));
    avrg_sem_coh{n} = nanstd(avrg_mean_coh{n}./sqrt(length(coherence{n})));
    end
    
    subplot(2,3,n)
    mu_pred1 = mean(pred1,1);
    mu_pred2 = mean(pred2,1);
    sem_pred1 = std(pred1,1)./(sqrt(size(pred1,1))); 
    sem_pred2 = std(pred2,1)./(sqrt(size(pred2,1))); 
    mu_pred2 = squeeze(mu_pred2);
    sem_pred2 = squeeze(sem_pred2);
    
    mu_Pright1 = mean(Pright1); 
    sem_Pright1 = std(Pright1)./(sqrt(length(Pright1)));
    mu_Pright2 = mean(Pright2); 
    sem_Pright2 = std(Pright2)./(sqrt(length(Pright2)));
    
    hold all; grid on
    x = linspace(min(interval1),max(interval1), ntrials_calib/length(interval1));
    lines = plot(x, mu_pred1, 'LineWidth', 3, 'Color', 'r');

        x2 = [x, fliplr(x)];
        inBetween = [mu_pred1 + sem_pred1, fliplr(mu_pred1 - sem_pred1)]; 
        xt = fill(x2, inBetween, c.err{3}, 'edgecolor', 'none'); 
        set(xt, 'facealpha', 0.3);

    dots = plot(interval1, mu_Pright1/(ntrials_calib/length(interval1)), 'o', 'MarkerSize', 5, 'MarkerFaceColor', c.corr{3}, 'MarkerEdgeColor', c.corr{3});
    errorbar(interval1, mu_Pright1/(ntrials_calib/length(interval1)), sem_Pright1/(ntrials_calib/length(interval1)), 'Color', c.corr{3},'LineStyle', 'none', 'LineWidth', 2)
    set(gca, 'YLim', [0 1],'YTick', [0:0.2:1], 'XTick', interval1, 'FontSize',13);

    title([(nat)], 'FontSize', 30);
    xlabel('Evidence level * dir.', 'FontSize', 20)
    ylabel('P(right)', 'FontSize', 20);
    
    subplot(2,3, n+3)
    hold all; grid on
        x2 = linspace(min(interval2),max(interval2), ntrials/length(interval2));
        lines = plot(x2, mu_pred2, 'LineWidth', 3, 'Color', 'r');  
        dots = plot(interval2, mu_Pright2/(ntrials/length(interval2)), 'o', 'MarkerSize', 5, 'MarkerFaceColor', c.corr{3}, 'MarkerEdgeColor', c.corr{3});
        errorbar(interval2, mu_Pright2/(ntrials/length(interval2)), sem_Pright2/(ntrials/length(interval2)), 'Color', c.corr{3},'LineStyle', 'none', 'LineWidth', 2)
        set(gca, 'XLim', [-3 3], 'XTick', interval2, 'YLim', [0 1],'YTick', 0:0.2:1, 'FontSize',15);       
        title(nat, 'FontSize', 30);
        xlabel('Evidence level * dir.', 'FontSize', 20)
        ylabel('P(right)', 'FontSize', 20);
    
end
%disp the stats
%%accuracy contrasts calibration
disp(['average accuracy calibration in PKU : ' num2str(mean(acc1(1,:)))]);
disp(['average accuracy calibration in UCL : ' num2str(mean(acc1(2,:)))]);
disp(['average accuracy calibration in NYU : ' num2str(nanmean(acc1(3,:)))]);
[H,P,CI,STATS] = ttest2(acc1(1,:), acc1(3,:)); 
disp(['P-value: ' num2str(P)]);

%%coherence contrasts main task
disp(['average coherence main task in PKU : ' num2str(mean(coherence{1}))]);
disp(['average cohernce main task in UCL : ' num2str(mean(coherence{2}))]);
disp(['average cohernce main task in NYU : ' num2str(nanmean(coherence{3}))]);
[H,P,CI,STATS] = ttest(coherence{1}, coherence{2}); 
disp(['P-value: ' num2str(P)]);

%%accuracy main task 
disp(['average accuracy main task in PKU : ' num2str(mean(accuracy{1}))]);
disp(['average accuracy main task in UCL : ' num2str(mean(accuracy{2}))]);
disp(['average accuracy main task in NYU : ' num2str(mean(accuracy{3}))]);
[H,P,CI,STATS] = ttest(accuracy{1}, accuracy{2}); 
disp(['P-value: ' num2str(P)]);

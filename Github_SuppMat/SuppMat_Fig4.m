%% supplementary material 2.3
fs = filesep;

culture = {'PKU'; 'UCL'};
sj_mat = {[403:418, 420:432, 434:435, 437:443, 445:459]; [25:76, 79]};
c.corr =  {[0.117, 0.447, 0.662],[0.560, 0.803, 0.898]};
c.err = {[0.909, 0.784, 0.772], [0.709, 0.090, 0.031]};
  
s05 = figure;
set(s05,'units','points','position',[10,10,1300,500])
    
for n = 1:length(culture)
    nat = culture{n};
    
    baseDir =  ['~' fs 'Dropbox' fs 'CulturalMetacognition_2020' fs 'DATA' fs 'EXP2' fs];
    dirData = [baseDir nat '_data' fs nat '_data' fs];

    filename = 'fMRI_pilotData_sub_';   
    suffix = '';
    subjects = sj_mat{n};
    cwd = pwd;
    
    a_corr_s{n} = []; %to store the conf for every ID
    a_corr_a{n} = [];
    a_err_s{n} = [];
    a_err_a{n} = [];
    da_corr_s{n} = [];
    da_corr_a{n} = [];
    da_err_s{n} = [];
    da_err_a{n} = [];
    
    for s = 1:length(subjects)
        datafile = [filename num2str(subjects(s)) suffix '_2.mat']; %actual task data
        cd(dirData);
        load(datafile, 'locDATA');
        
        precoh = locDATA.dots_coherence';
        postcoh = locDATA.post_coherence';
        dir = locDATA.dots_direction/360;
        dir(dir==0.5) = -1;
        action = locDATA.button_response - 1;
        conf = locDATA.mouse_response;
        conf_dir = locDATA.mouse_response_dir;
        transformed_action = action;
        transformed_action(action == 0) = -1;
        acc = dir == transformed_action;
        conf_adv = locDATA.conf_adv;
        a_adv = locDATA.a_adv;
        a_adv(a_adv == 99) = NaN;
        a_adv = a_adv-1; %back from 1/2 to 0/1
        agree = a_adv == action;
        
        % Add indicator variables for pre/post confidence
        coherence = unique(precoh);
        
        for i = 1:3
            precoh_index(locDATA.dots_coherence==coherence(i))=i;
            postcoh_index(locDATA.post_coherence==coherence(i))=i;
        end

        conf_adv(conf_adv == 99.000) = NaN;
       
        for i = 1:length(conf_dir)
            if dir(i) == -1
                conf_inv(i) = 1 - conf_dir(i);
                conf_adv_inv(i) = 1 - conf_adv(i);
            else
                conf_inv(i) = conf_dir(i);
                conf_adv_inv(i) = conf_adv(i);
            end
            
            if conf_adv(i) < 0.5
                conf_adv2(i) = 1 - conf_adv(i);
            else
                conf_adv2(i) = conf_adv(i);
                end
        end
  
        quant_conf_adv = quantile(conf_adv2, [.333, .666]);
        conf_adv_level = discretize(conf_adv2, [0.5 quant_conf_adv 1]); 
        
        a_corr_s{n} = [a_corr_s{n}, conf_inv(acc == 1 & agree == 1 & locDATA.condition == 1)];
        a_corr_a{n} = [a_corr_a{n}, conf_adv_inv(acc == 1 & agree == 1 & locDATA.condition == 1)];
        
        da_corr_s{n}= [da_corr_s{n}, conf_inv(acc == 1 & agree == 0 & locDATA.condition == 1)];
        da_corr_a{n} = [da_corr_a{n}, conf_adv_inv(acc == 1 & agree == 0 & locDATA.condition == 1)];
        
        a_err_s{n} = [a_err_s{n}, conf_inv(acc == 0 & agree == 1 & locDATA.condition == 1)];
        a_err_a{n} = [a_err_a{n}, conf_adv_inv(acc == 0 & agree == 1 & locDATA.condition == 1)];
        
        da_err_s{n} = [da_err_s{n}, conf_inv(acc == 0 & agree == 0 & locDATA.condition == 1)];
        da_err_a{n} = [da_err_a{n}, conf_adv_inv(acc == 0 & agree == 0 & locDATA.condition == 1)];
        
    end
    filename = 'regression_betas_';
    betaDir = [dirData nat '_betas' fs];
    
    for cor = 1:2
        correct = {'err_'; 'corr_'};
        for agr = 1:2
            agree = {'da'; 'a'};
            suffix = correct{cor}; 
            suffix2 = agree{agr};
 
            cd(betaDir);
            datafile = [filename suffix suffix2 '_' nat '.csv'];
            dat{n,cor, agr} = readtable(datafile);
            cd(cwd);
        end
    end
end

subplot(1,2,1)
    
%plot COM-line
hline1 = line([0,1],[0.5,0.5], 'linestyle', '--', 'color', [0 0 0], 'linewidth', 0.2);
hline2 = line([0.5,0.5],[0,1], 'linestyle', '--', 'color', [0 0 0], 'linewidth', 0.2);
hold all

width = [9,5]; 
for n = 1:2 
    y = erf(a_corr_a{n});
    p = polyfit(a_corr_a{n},a_corr_s{n},1);
    f = polyval(p, y);
    ae{n} = plot(y,f,'-', 'Color', c.corr{1}, 'LineWidth', width(n));
    
    y = erf(da_corr_a{n});
    p = polyfit(da_corr_a{n}, da_corr_s{n},1);
    f = polyval(p, y);
    ae{n} = plot(y, f, '-','Color', c.corr{2}, 'LineWidth', width(n));

    y = erf(a_err_a{n});
    p = polyfit(a_err_a{n}, a_err_s{n},1);
    f = polyval(p, y);
    ae{n} = plot(y, f, '-', 'Color', c.err{1}, 'LineWidth', width(n));
    
    y = erf(da_err_a{n});
    p = polyfit(da_err_a{n}, da_err_s{n},1);
    f = polyval(p, y);
    ae{n} = plot(y, f,'-', 'Color', c.err{2}, 'Linewidth',width(n));
    
    fake_line{n} = plot(y+20, f+20, '-', 'Color', 'k', 'Linewidth', width(n)); 
end
text(0.54, 0.97, 'correct, agree', 'Color', c.corr{1}, 'FontSize', 20)
text(0.18, 0.85, 'correct, disagree', 'Color', c.corr{2}, 'FontSize', 20)
text(0.25, 0.4, 'error, agree', 'Color', c.err{1}, 'FontSize', 20)
text(0.52, 0.55, 'error, disagree', 'Color', c.err{2}, 'FontSize', 20)

ylabel([{'Participant confidence in'}; {'objectively correct direction'}], 'Fontsize', 30);
xlabel([{'Adviser confidence in'}; {'objectively correct direction'}], 'FontSize', 30);
set(gca, 'YLim', 0:1, 'XLim', 0:1, 'XTick', 0:0.2:1, 'XTickLabel', [0:0.2:1],'YTick', [0:0.2:1], 'YTickLabel', [0:0.2:1], 'FontSize', 20);  

[lgd, handle] = legend([fake_line{1}, fake_line{2}], {'PKU', 'UCL'},'location', 'SouthEast');
linehandle = findobj(handle, 'type', 'line');
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);

subplot(1,2,2)
xpos = [1 2 4 5];
hold all
errorbar(xpos,  [dat{1,1,2}{3,2} dat{1,1,1}{3,2} dat{2,1,2}{3,2} dat{2,1,1}{3,2}],[dat{1,1,2}{3+5,2} dat{1,1,1}{3+5,2} dat{2,1,2}{3+5,2} dat{2,1,1}{3+5,2}], '.', 'Color', c.err{2}*0.3, 'LineWidth', 2);
errorbar(xpos,  [dat{1,2,2}{3,2} dat{1,2,1}{3,2} dat{2,2,2}{3,2}  dat{2,2,1}{3,2}] , [dat{1,2,2}{3+5,2} dat{1,2,1}{3+5,2} dat{2,2,2}{3+5,2}  dat{2,2,1}{3+5,2}], '.', 'Color', c.corr{1}*0.3, 'LineWidth', 2);
x = [1 4];
err_a = plot(x,[dat{1,1,2}{3,2} dat{2,1,2}{3,2}] , 'o', 'MarkerSize',10,  'MarkerFaceColor', c.err{1}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{1});%dat{PKU, ERR, AGR},dat{UCL, ERR, AGR}  
corr_a = plot(x, [dat{1,2,2}{3,2} dat{2,2,2}{3,2}] , 'o', 'MarkerSize', 10,'MarkerFaceColor', c.corr{1}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{1});%dat{PKU, CORR, AGR}, dat{UCL, CORR, AGR}
x = [2 5];
err_da = plot(x,[dat{1,1,1}{3,2} dat{2,1,1}{3,2}] , 'o', 'MarkerSize',10,  'MarkerFaceColor', c.err{2}, 'Linewidth', 3, 'MarkerEdgeColor', c.err{2});%dat{PKU,ERR, DISAGR}, dat{UCL, ERR, DISAGR},
corr_da = plot(x, [dat{1,2,1}{3,2} dat{2,2,1}{3,2} ] , 'o', 'MarkerSize', 10,'MarkerFaceColor', c.corr{2}, 'Linewidth', 3, 'MarkerEdgeColor', c.corr{2});%dat{PKU,CORR, DISAGR},  dat{UCL, CORR, DISAGR}, 
    
grid off; box off
hline1 = line([0 6], [0,0], 'linestyle', '-', 'color', [0 0 0], 'linewidth', 0.7);

[lgd, handle] = legend([corr_a, err_da],{'correct', 'error'},'location', 'SouthWest');
linehandle = findobj(handle, 'type', 'line');
set(linehandle, 'LineWidth',6)
legend boxoff
texthandle = findobj(handle, 'type', 'text');
set(texthandle,'FontSize',20);

set(gca, 'XLim', [0 6], 'XTick', xpos, 'XTickLabel',{'agree', 'disagree'}, 'YLim', [-0.2 0.2],'YTick', [-0.2:0.1:0.2], 'FontSize',20);
ylabel('Impact confidence adviser (a.u.)', 'FontSize', 20);

text(xpos(2)-0.75, -0.24, 'PKU', 'FontSize', 20)
text(xpos(4)-0.75, -0.24, ' UCL', 'FontSize', 20)


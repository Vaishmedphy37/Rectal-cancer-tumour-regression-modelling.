clc;
clear;
close all;
%% =========================================================
% BREAKPOINTS vs REPLANNING FRACTION
% ---------------------------------------------------------
% Evaluates the relationship between the optimised tumour
% regression breakpoints (BP1 and BP2) and the identified
% re-planning fraction across patients.
% =========================================================

%% =========================================================
% LOAD RESULTS
% =========================================================

T = readtable('Tumour_Kinetics_Results.xlsx');


%% =========================================================
% CORRELATION ANALYSIS
% TOOLBOX-FREE PEARSON + SPEARMAN
% =========================================================


%% ---------------------------------------------------------
% BP1 vs IdealReplanFx
% ---------------------------------------------------------

idx1 = ~isnan(T.BP1) & ~isnan(T.IdealReplanFx);

x1 = T.BP1(idx1);
y1 = T.IdealReplanFx(idx1);

n1 = length(x1);
df1 = n1 - 2;


% =========================================================
% PEARSON CORRELATION: BP1 vs IdealReplanFx
% =========================================================

x1_mean = mean(x1);
y1_mean = mean(y1);

R_pearson_BP1 = ...
    sum((x1-x1_mean).*(y1-y1_mean)) / ...
    sqrt(sum((x1-x1_mean).^2) * ...
         sum((y1-y1_mean).^2));


% Pearson t-statistic

t_pearson_BP1 = ...
    R_pearson_BP1 * sqrt(df1/(1-R_pearson_BP1^2));


% Two-tailed p-value

P_pearson_BP1 = ...
    betainc(df1/(df1+t_pearson_BP1^2), ...
            df1/2, ...
            0.5);


% =========================================================
% SPEARMAN CORRELATION: BP1 vs IdealReplanFx
% =========================================================

% Proper average ranks for tied observations

rank_x1 = averageRank(x1);
rank_y1 = averageRank(y1);


% Pearson correlation of the ranks

rank_x1_mean = mean(rank_x1);
rank_y1_mean = mean(rank_y1);

R_spearman_BP1 = ...
    sum((rank_x1-rank_x1_mean).*(rank_y1-rank_y1_mean)) / ...
    sqrt(sum((rank_x1-rank_x1_mean).^2) * ...
         sum((rank_y1-rank_y1_mean).^2));


% Spearman t-statistic

t_spearman_BP1 = ...
    R_spearman_BP1 * sqrt(df1/(1-R_spearman_BP1^2));


% Two-tailed p-value

P_spearman_BP1 = ...
    betainc(df1/(df1+t_spearman_BP1^2), ...
            df1/2, ...
            0.5);



%% ---------------------------------------------------------
% BP2 vs IdealReplanFx
% ---------------------------------------------------------

idx2 = ~isnan(T.BP2) & ~isnan(T.IdealReplanFx);

x2 = T.BP2(idx2);
y2 = T.IdealReplanFx(idx2);

n2 = length(x2);
df2 = n2 - 2;


% =========================================================
% PEARSON CORRELATION: BP2 vs IdealReplanFx
% =========================================================

x2_mean = mean(x2);
y2_mean = mean(y2);

R_pearson_BP2 = ...
    sum((x2-x2_mean).*(y2-y2_mean)) / ...
    sqrt(sum((x2-x2_mean).^2) * ...
         sum((y2-y2_mean).^2));


% Pearson t-statistic

t_pearson_BP2 = ...
    R_pearson_BP2 * sqrt(df2/(1-R_pearson_BP2^2));


% Two-tailed p-value

P_pearson_BP2 = ...
    betainc(df2/(df2+t_pearson_BP2^2), ...
            df2/2, ...
            0.5);


% =========================================================
% SPEARMAN CORRELATION: BP2 vs IdealReplanFx
% =========================================================

% Proper average ranks for tied observations

rank_x2 = averageRank(x2);
rank_y2 = averageRank(y2);


% Pearson correlation of the ranks

rank_x2_mean = mean(rank_x2);
rank_y2_mean = mean(rank_y2);

R_spearman_BP2 = ...
    sum((rank_x2-rank_x2_mean).*(rank_y2-rank_y2_mean)) / ...
    sqrt(sum((rank_x2-rank_x2_mean).^2) * ...
         sum((rank_y2-rank_y2_mean).^2));


% Spearman t-statistic

t_spearman_BP2 = ...
    R_spearman_BP2 * sqrt(df2/(1-R_spearman_BP2^2));


% Two-tailed p-value

P_spearman_BP2 = ...
    betainc(df2/(df2+t_spearman_BP2^2), ...
            df2/2, ...
            0.5);



%% =========================================================
% DISPLAY RESULTS IN COMMAND WINDOW
% =========================================================

fprintf('\n');
fprintf('=====================================================\n');
fprintf('       BP vs IDEAL REPLAN CORRELATION ANALYSIS\n');
fprintf('=====================================================\n');


fprintf('\nBP1 vs Ideal Replan Fraction\n');
fprintf('N = %d\n', n1);
fprintf('-----------------------------------------------------\n');
fprintf('Pearson  : r = %.3f, p = %.4f\n', ...
    R_pearson_BP1, P_pearson_BP1);
fprintf('Spearman : r = %.3f, p = %.4f\n', ...
    R_spearman_BP1, P_spearman_BP1);


fprintf('\nBP2 vs Ideal Replan Fraction\n');
fprintf('N = %d\n', n2);
fprintf('-----------------------------------------------------\n');
fprintf('Pearson  : r = %.3f, p = %.4f\n', ...
    R_pearson_BP2, P_pearson_BP2);
fprintf('Spearman : r = %.3f, p = %.4f\n', ...
    R_spearman_BP2, P_spearman_BP2);


fprintf('\n=====================================================\n');



%% =========================================================
% PAGE 1 — SCATTER PLOTS
% =========================================================

figure('Color','w',...
       'Position',[100 100 1400 600]);


phenotypes = unique(T.Phenotype);

colors = lines(length(phenotypes));


%% =========================================================
% PANEL 1: BP1 vs IdealReplanFx
% =========================================================

subplot(1,2,1);

hold on;


% ---------------------------------------------------------
% Plot each phenotype
% ---------------------------------------------------------

for k = 1:length(phenotypes)

    idx = strcmp(T.Phenotype,phenotypes{k});

    scatter(T.BP1(idx),...
            T.IdealReplanFx(idx),...
            120,...
            colors(k,:),...
            'filled',...
            'DisplayName',phenotypes{k});

end


% ---------------------------------------------------------
% Patient labels
% ---------------------------------------------------------

for i = 1:height(T)

    if ~isnan(T.BP1(i)) && ...
       ~isnan(T.IdealReplanFx(i))

        text(T.BP1(i)+0.15,...
             T.IdealReplanFx(i),...
             char(T.Patient(i)),...
             'FontSize',8);

    end

end


xlabel('First Breakpoint Fraction (BP1)');
ylabel('Ideal Replan Fraction');

title('BP1 vs Ideal Replan Fraction');

grid on;

legend('Location','eastoutside');


% ---------------------------------------------------------
% Reference line: BP1 = IdealReplanFx
% ---------------------------------------------------------

lims1 = [ ...
    min([T.BP1;T.IdealReplanFx],[],'omitnan')-1,...
    max([T.BP1;T.IdealReplanFx],[],'omitnan')+1 ...
    ];

plot(lims1,...
     lims1,...
     '--k',...
     'LineWidth',1,...
     'HandleVisibility','off');


% ---------------------------------------------------------
% Correlation box
% ---------------------------------------------------------

corrText1 = sprintf( ...
    ['Pearson:  r = %.3f, p = %.4f\n' ...
     'Spearman: r = %.3f, p = %.4f'], ...
    R_pearson_BP1,...
    P_pearson_BP1,...
    R_spearman_BP1,...
    P_spearman_BP1);


text(0.05,...
     0.95,...
     corrText1,...
     'Units','normalized',...
     'VerticalAlignment','top',...
     'FontSize',10,...
     'BackgroundColor','white',...
     'EdgeColor','black',...
     'Margin',6);


% ---------------------------------------------------------
% Phenotype / patient ID explanation
% ---------------------------------------------------------

annotation('textbox',...
    [0.08 0.80 0.10 0.10],...
    'String',{...
    'Colour = Phenotype';...
    'Label = Patient ID'},...
    'FitBoxToText','on',...
    'BackgroundColor','white');



%% =========================================================
% PANEL 2: BP2 vs IdealReplanFx
% =========================================================

subplot(1,2,2);

hold on;


% ---------------------------------------------------------
% Plot each phenotype
% ---------------------------------------------------------

for k = 1:length(phenotypes)

    idx = strcmp(T.Phenotype,phenotypes{k});

    scatter(T.BP2(idx),...
            T.IdealReplanFx(idx),...
            120,...
            colors(k,:),...
            'filled',...
            'DisplayName',phenotypes{k});

end


% ---------------------------------------------------------
% Patient labels
% ---------------------------------------------------------

for i = 1:height(T)

    if ~isnan(T.BP2(i)) && ...
       ~isnan(T.IdealReplanFx(i))

        text(T.BP2(i)+0.15,...
             T.IdealReplanFx(i),...
             char(T.Patient(i)),...
             'FontSize',8);

    end

end


xlabel('Second Breakpoint Fraction (BP2)');
ylabel('Ideal Replan Fraction');

title('BP2 vs Ideal Replan Fraction');

grid on;


% ---------------------------------------------------------
% Reference line: BP2 = IdealReplanFx
% ---------------------------------------------------------

lims2 = [ ...
    min([T.BP2;T.IdealReplanFx],[],'omitnan')-1,...
    max([T.BP2;T.IdealReplanFx],[],'omitnan')+1 ...
    ];

plot(lims2,...
     lims2,...
     '--k',...
     'LineWidth',1,...
     'HandleVisibility','off');


% ---------------------------------------------------------
% Correlation box
% ---------------------------------------------------------

corrText2 = sprintf( ...
    ['Pearson:  r = %.3f, p = %.4f\n' ...
     'Spearman: r = %.3f, p = %.4f'], ...
    R_pearson_BP2,...
    P_pearson_BP2,...
    R_spearman_BP2,...
    P_spearman_BP2);


text(0.05,...
     0.95,...
     corrText2,...
     'Units','normalized',...
     'VerticalAlignment','top',...
     'FontSize',10,...
     'BackgroundColor','white',...
     'EdgeColor','black',...
     'Margin',6);



%% =========================================================
% LOCAL FUNCTION — AVERAGE RANKS
% Handles tied observations correctly
% =========================================================

function ranks = averageRank(x)

    % Sort data
    [sortedX, order] = sort(x);

    ranks = zeros(size(x));

    i = 1;

    while i <= length(x)

        % Find all observations tied with current value
        j = i;

        while j < length(x) && sortedX(j+1) == sortedX(i)
            j = j + 1;
        end

        % Average rank for tied observations
        avgRank = (i+j)/2;

        % Assign average rank
        ranks(order(i:j)) = avgRank;

        i = j + 1;

    end

end
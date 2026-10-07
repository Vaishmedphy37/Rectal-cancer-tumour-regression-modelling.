clc;
clear;
close all;

%% =========================================================
% LOAD RESULTS
%% =========================================================
T = readtable('Tumour_Kinetics_Results.xlsx');

phenotypes = unique(T.Phenotype);
colors = lines(length(phenotypes));

%% =========================================================
% CORRELATIONS: k1 vs k2, k1 vs k3, k2 vs k3
% (toolbox-free: corrcoef + manual tied-rank Spearman + betainc p-value)
%% =========================================================
pairs = {'k1','k2'; 'k1','k3'; 'k2','k3'};
pairLabels = {'k_1 vs k_2','k_1 vs k_3','k_2 vs k_3'};

fprintf('\n=========================================\n');
fprintf('PHASE-RATE PAIRWISE CORRELATIONS\n');
fprintf('=========================================\n');
fprintf('%-10s %10s %10s %10s %10s %6s\n','Pair','Pearson r','p (Pearson)','Spearman rho','p (Spearman)','n');

corrResults = table('Size',[3 7], ...
    'VariableTypes',{'string','string','double','double','double','double','double'}, ...
    'VariableNames',{'Var1','Var2','PearsonR','PearsonP','SpearmanRho','SpearmanP','n'});

for j = 1:3
    x = T.(pairs{j,1});
    y = T.(pairs{j,2});
    n = numel(x);

    Rmat = corrcoef(x, y);
    rP = Rmat(1,2);
    pP = pValueFromR(rP, n);

    xr = tiedRank(x);
    yr = tiedRank(y);
    Rmat2 = corrcoef(xr, yr);
    rS = Rmat2(1,2);
    pS = pValueFromR(rS, n);

    corrResults.Var1(j) = pairs{j,1};
    corrResults.Var2(j) = pairs{j,2};
    corrResults.PearsonR(j) = rP;
    corrResults.PearsonP(j) = pP;
    corrResults.SpearmanRho(j) = rS;
    corrResults.SpearmanP(j) = pS;
    corrResults.n(j) = n;

    fprintf('%-10s %10.3f %10.4f %10.3f %10.4f %6d\n', ...
        [pairs{j,1} ' vs ' pairs{j,2}], rP, pP, rS, pS, n);
end

writetable(corrResults, 'k_pairwise_Correlation_Results.xlsx');

%% =========================================================
% 3-PANEL SIDE-BY-SIDE SCATTER PLOT
%% =========================================================
fig = figure('Color','w', 'Position',[100 100 1500 500]);

for j = 1:3
    subplot(1,3,j)
    hold on
    xvar = pairs{j,1};
    yvar = pairs{j,2};

    for k = 1:length(phenotypes)
        idx = strcmp(T.Phenotype, phenotypes{k});
        scatter(T.(xvar)(idx), T.(yvar)(idx), 100, colors(k,:), ...
            'filled', 'DisplayName', phenotypes{k});
    end

    xlabel(strrep(xvar,'k','k_'))
    ylabel(strrep(yvar,'k','k_'))
    title(pairLabels{j})
    grid on

    r = corrResults.PearsonR(j);
    p = corrResults.PearsonP(j);
    rho = corrResults.SpearmanRho(j);
    ps = corrResults.SpearmanP(j);
    annotation('textbox', [0.02+(j-1)*0.33, 0.75, 0.28, 0.12], ...
        'String', {sprintf('Pearson r = %.3f (p=%.3f)', r, p), ...
                    sprintf('Spearman \\rho = %.3f (p=%.3f)', rho, ps)}, ...
        'FitBoxToText','on', 'BackgroundColor','white', 'FontSize', 8);

    if j == 3
        legend('Location','eastoutside')
    end
end

sgtitle('Relationships Between Phase-Specific Regression Rates')

exportgraphics(fig, 'Fig_k1k2k3_Pairwise.pdf', 'ContentType','vector');
exportgraphics(fig, 'Fig_k1k2k3_Pairwise.png', 'Resolution',300);

disp('Saved: Fig_k1k2k3_Pairwise.pdf / .png');
disp('Saved: k_pairwise_Correlation_Results.xlsx');

%% =========================================================
% LOCAL FUNCTIONS (no Statistics Toolbox required)
%% =========================================================

function p = pValueFromR(r, n)
    if n <= 2
        p = NaN;
        return
    end
    df = n - 2;
    t  = r * sqrt(df / (1 - r^2));
    x  = df / (df + t^2);
    p  = betainc(x, df/2, 0.5);
end

function r = tiedRank(x)
    [~, sortIdx]    = sort(x);
    [~, invSortIdx] = sort(sortIdx);
    ranks = (1:numel(x))';
    r = ranks(invSortIdx);

    uniqueVals = unique(x);
    for i = 1:numel(uniqueVals)
        tieIdx = (x == uniqueVals(i));
        if sum(tieIdx) > 1
            r(tieIdx) = mean(r(tieIdx));
        end
    end
end
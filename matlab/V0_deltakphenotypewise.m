clc;
clear;
close all;
%% =========================================================
% BASELINE TUMOUR VOLUME vs REGRESSION DYNAMICITY
% ---------------------------------------------------------
% Evaluates the relationship between baseline tumour volume
% (V0) and regression dynamicity (\Delta k), including
% overall and phenotype-wise correlation analyses.
% =========================================================
%% =========================================================
% LOAD RESULTS
%% =========================================================
T = readtable('Tumour_Kinetics_Results.xlsx');

%% =========================================================
% SCATTER PLOT: V0 vs Delta k, coloured by Phenotype
%% =========================================================
fig = figure('Color','w', 'Position',[100 100 1200 700]);
hold on

phenotypes = unique(T.Phenotype);
colors = lines(length(phenotypes));

for k = 1:length(phenotypes)
    idx = strcmp(T.Phenotype,phenotypes{k});
    scatter(T.V0(idx), T.Delta_k(idx), 120, colors(k,:), ...
        'filled', 'DisplayName', phenotypes{k});
end
for i = 1:height(T)
    text(T.V0(i)+2, T.Delta_k(i), char(T.Patient(i)), 'FontSize', 8);
end

xlabel('Initial Tumor Volume (V0)')
ylabel('\Delta k')
title('V0 vs Dynamicity by Phenotype')
grid on
legend('Location','eastoutside')

%% =========================================================
% OVERALL (POOLED) CORRELATION — all patients, ignoring phenotype
%% =========================================================
x_all = T.V0;
y_all = T.Delta_k;
n_all = numel(x_all);

R_all = pearsonR(x_all, y_all);
p_all = pValueFromR(R_all, n_all);

fprintf('\n================================\n');
fprintf('OVERALL (POOLED) V0 vs DELTA K\n');
fprintf('================================\n');
fprintf('n = %d\n', n_all);
fprintf('R = %.3f\n', R_all);
fprintf('p = %.4f\n', p_all);

annotation('textbox', [0.14 0.80 0.15 0.08], ...
    'String', sprintf('Pooled: R = %.2f (p = %.3f, n = %d)', R_all, p_all, n_all), ...
    'FitBoxToText','on', 'BackgroundColor','white');

exportgraphics(fig, 'Fig_V0_vs_DeltaK_byPhenotype.pdf', 'ContentType','vector');
exportgraphics(fig, 'Fig_V0_vs_DeltaK_byPhenotype.png', 'Resolution',300);

%% =========================================================
% PHENOTYPE-WISE CORRELATION
% Minimum n = 3 required to compute a meaningful R and p-value;
% smaller subgroups are flagged rather than reported.
%% =========================================================
disp(' ')
disp('================================')
disp('PHENOTYPE-WISE V0 vs DELTA K')
disp('================================')

MIN_N = 3;

for k = 1:length(phenotypes)
    idx = strcmp(T.Phenotype,phenotypes{k});
    x = T.V0(idx);
    y = T.Delta_k(idx);
    n = sum(idx);

    fprintf('\n%s\n', phenotypes{k});
    fprintf('n = %d\n', n);

    if n < MIN_N
        fprintf('R = NOT COMPUTED (n < %d; insufficient for a meaningful correlation)\n', MIN_N);
        continue
    end

    R = pearsonR(x, y);
    p = pValueFromR(R, n);
    fprintf('R = %.3f\n', R);
    fprintf('p = %.4f\n', p);
end

disp(' ')
disp('Note: subgroups with n < 3 are reported as sample size only;')
disp('no correlation coefficient or p-value is computed for them.')

%% =========================================================
% LOCAL FUNCTIONS (no Statistics Toolbox required)
%% =========================================================

function R = pearsonR(x, y)
    R = sum((x-mean(x)).*(y-mean(y))) / ...
        sqrt(sum((x-mean(x)).^2) * sum((y-mean(y)).^2));
end

function p = pValueFromR(r, n)
    if n <= 2
        p = NaN;
        return
    end
    df = n - 2;
    t  = r * sqrt(df / (1 - r^2));
    xval = df / (df + t^2);
    p  = betainc(xval, df/2, 0.5);
end
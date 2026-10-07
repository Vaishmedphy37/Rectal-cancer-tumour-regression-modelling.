%% =========================================================
% LOAD RESULTS
%% =========================================================
T = readtable('Tumour_Kinetics_Results.xlsx');
 
%% =========================================================
% COMPUTE V12 FROM PIECEWISE-EXPONENTIAL FIT
% v(fx) = V(fx)/V0, using each patient's k1,k2,k3,BP1,BP2
% =========================================================
fx_target = 12;
V12 = zeros(height(T),1);
for i = 1:height(T)
    k1  = T.k1(i);
    k2  = T.k2(i);
    k3  = T.k3(i);
    BP1 = T.BP1(i);
    BP2 = T.BP2(i);
 
    if fx_target <= BP1
        v = exp(-k1*fx_target);
    elseif fx_target <= BP2
        v = exp(-k1*BP1) * exp(-k2*(fx_target-BP1));
    else
        v = exp(-k1*BP1) * exp(-k2*(BP2-BP1)) * exp(-k3*(fx_target-BP2));
    end
    V12(i) = v;
end
T.V12 = V12;
 
phenotypes = unique(T.Phenotype);
colors = lines(length(phenotypes));
 
%% =========================================================
% FIGURE 1: V0 vs IdealReplanFx
%% =========================================================
fig1 = figure('Color','w', 'Position',[100 100 800 650]);
hold on
for k = 1:length(phenotypes)
    idx = strcmp(T.Phenotype,phenotypes{k});
    scatter(T.V0(idx), T.IdealReplanFx(idx), 120, colors(k,:), ...
        'filled', 'DisplayName', phenotypes{k});
end
for i = 1:height(T)
    text(T.V0(i)+1, T.IdealReplanFx(i), char(T.Patient(i)), 'FontSize', 8);
end
xlabel('Initial Tumour Volume (V0)')
ylabel('Ideal Replan Fraction')
title('V0 vs Ideal Replan Fraction')
grid on
legend('Location','eastoutside')
yline(20, '--r', 'Fx20 Boundary', 'LineWidth', 1.5);
annotation('textbox', [0.08 0.80 0.10 0.10], ...
    'String', {'Colour = Phenotype';'Label = Patient ID'}, ...
    'FitBoxToText','on', 'BackgroundColor','white');
 
exportgraphics(fig1, 'Fig1_V0_vs_IdealReplanFx.pdf', 'ContentType','vector');
exportgraphics(fig1, 'Fig1_V0_vs_IdealReplanFx.png', 'Resolution',300);
 
%% =========================================================
% FIGURE 2: V12 vs IdealReplanFx
%% =========================================================
fig2 = figure('Color','w', 'Position',[100 100 800 650]);
hold on
for k = 1:length(phenotypes)
    idx = strcmp(T.Phenotype,phenotypes{k});
    scatter(T.V12(idx), T.IdealReplanFx(idx), 120, colors(k,:), ...
        'filled', 'DisplayName', phenotypes{k});
end
for i = 1:height(T)
    text(T.V12(i)+0.01, T.IdealReplanFx(i), char(T.Patient(i)), 'FontSize', 8);
end
xlabel('Mid-course Relative Volume at Fx12 (V12)')
ylabel('Ideal Replan Fraction')
title('V12 vs Ideal Replan Fraction')
grid on
legend('Location','eastoutside')
yline(20, '--r', 'Fx20 Boundary', 'LineWidth', 1.5);
 
exportgraphics(fig2, 'Fig2_V12_vs_IdealReplanFx.pdf', 'ContentType','vector');
exportgraphics(fig2, 'Fig2_V12_vs_IdealReplanFx.png', 'Resolution',300);
 
%% =========================================================
% FIGURE 3: Vmid vs IdealReplanFx
%% =========================================================
fig3 = figure('Color','w', 'Position',[100 100 800 650]);
hold on
for k = 1:length(phenotypes)
    idx = strcmp(T.Phenotype,phenotypes{k});
    scatter(T.Vmid(idx), T.IdealReplanFx(idx), 120, colors(k,:), ...
        'filled', 'DisplayName', phenotypes{k});
end
for i = 1:height(T)
    text(T.Vmid(i)+0.01, T.IdealReplanFx(i), char(T.Patient(i)), 'FontSize', 8);
end
xlabel('Mid-treatment Relative Volume (Vmid)')
ylabel('Ideal Replan Fraction')
title('Vmid vs Ideal Replan Fraction')
grid on
legend('Location','eastoutside')
yline(20, '--r', 'Fx20 Boundary', 'LineWidth', 1.5);
 
exportgraphics(fig3, 'Fig3_Vmid_vs_IdealReplanFx.pdf', 'ContentType','vector');
exportgraphics(fig3, 'Fig3_Vmid_vs_IdealReplanFx.png', 'Resolution',300);
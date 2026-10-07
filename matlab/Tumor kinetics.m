clc;
clear;
close all;

%% =========================================================
% TUMOUR REGRESSION KINETICS - PIECEWISE EXPONENTIAL MODEL
% ---------------------------------------------------------
% This script analyses anonymised tumour volume data from patients
% undergoing long-course chemoradiotherapy for rectal cancer.
%
% The workflow includes:
%   - Tumour volume normalisation
%   - PCHIP smoothing
%   - Breakpoint optimisation
%   - Phase-wise exponential regression
%   - Regression kinetic biomarkers (k1, k2, k3)
%   - Dynamicity and phenotype classification
%   - Model performance evaluation using RMSE and R^2
%   - Exploration of potential adaptive replanning windows
%
% Input:
%   P1.csv ... P25.csv
%   Column 1 = treatment fraction
%   Column 2 = tumour volume
%
% Output:
%   Tumour_Kinetics_Results.xlsx
%   PatientPage_*.pdf
%
% Patient data are represented using anonymised IDs (P1, P2, ...).
%% =========================================================

%% =========================================================
% FILES (anonymised)
%% =========================================================
nPatients = 25;

files = arrayfun(@(i) sprintf('P%d.csv', i), 1:nPatients, ...
                 'UniformOutput', false);
names = arrayfun(@(i) sprintf('P%d', i), 1:nPatients, ...
                 'UniformOutput', false);

Results = table();

figCount = 0;
%% =========================================================
% MAIN LOOP
%% =========================================================

for p = 1:nPatients

    filename = files{p};

    patientName = names{p};

    fprintf('\nProcessing %s\n',patientName);
    %% ==========================================
    % READ CSV
    %% ==========================================
    T = readtable(filename);

    t = T{:,1};
    V = T{:,2};
    [t,idx] = sort(t);
    V = V(idx);
    %% ==========================================
    % NORMALIZATION
    %% ==========================================

    V0 = V(1);

    Vn = V./V0;

    lnV = log(Vn);
    %% ==========================================
    % PCHIP SMOOTHING
    %% ==========================================

    tt = linspace(min(t),max(t),300);

    lnV_s = interp1(t,lnV,tt,'pchip');
    %% ==========================================
    % BREAKPOINT OPTIMIZATION
    %% ==========================================
    bp1_range = 4:0.1:10;
    bp2_range = 10:0.1:18;

    bestSSE = inf;

    best_bp1 = NaN;
    best_bp2 = NaN;

    for bp1 = bp1_range

        for bp2 = bp2_range

            if bp2 <= bp1
                continue;
            end
            %% Phase masks

            id1 = tt <= bp1;
            id2 = tt > bp1 & tt <= bp2;
            id3 = tt > bp2;

            if sum(id1)<5 || sum(id2)<5 || sum(id3)<5
                continue;
            end
            %% ----- Phase 1 -----

            c1 = polyfit(tt(id1),lnV_s(id1),1);

            k1 = -c1(1);

            lnV_fit = nan(size(tt));

            lnV_fit(id1) = polyval(c1,tt(id1));
            %% ----- Phase 2 -----

            lnV_bp1 = polyval(c1,bp1);

            x2 = tt(id2)-bp1;

            y2 = lnV_s(id2)-lnV_bp1;

            c2 = polyfit(x2,y2,1);

            k2 = -c2(1);

            lnV_fit(id2) = lnV_bp1 + c2(1)*x2;

            %% ----- Phase 3 -----

            lnV_bp2 = lnV_bp1 + c2(1)*(bp2-bp1);

            x3 = tt(id3)-bp2;

            y3 = lnV_s(id3)-lnV_bp2;

            c3 = polyfit(x3,y3,1);

            k3 = -c3(1);

            lnV_fit(id3) = lnV_bp2 + c3(1)*x3;

            %% ----- SSE -----

            SSE = sum((lnV_s - lnV_fit).^2,'omitnan');

            if SSE < bestSSE

                bestSSE = SSE;

                best_bp1 = bp1;
                best_bp2 = bp2;

                best_k1 = k1;
                best_k2 = k2;
                best_k3 = k3;

                best_fit = lnV_fit;

            end

        end

    end
    RMSE = sqrt(mean((lnV_s - best_fit).^2,'omitnan'));

    SS_res = sum((lnV_s - best_fit).^2,'omitnan');

    SS_tot = sum((lnV_s - mean(lnV_s)).^2,'omitnan');

    R2 = 1 - SS_res/SS_tot;
    %% ==========================================
    % KINETIC BIOMARKERS
    %% ==========================================

    k1 = best_k1;

    k2 = best_k2;

    k3 = best_k3;

    %% ----- Phase changes -----

    k12 = abs(k2-k1);

    k23 = abs(k3-k2);

    %% ----- Overall dynamicity -----

    Delta_k = max([k1 k2 k3]) - min([k1 k2 k3]);
    %% ==========================================
    % DYNAMICITY
    %% ==========================================

    if  Delta_k < 0.0354

        Dynamicity = "Stable";

    elseif (0.0354 < Delta_k) && (Delta_k < 0.0475)

        Dynamicity = "Moderate";

    elseif Delta_k > 0.0475

        Dynamicity = "Highly Dynamic";

    end
    %% ==========================================
    % PHENOTYPE
    %% ==========================================

    tol = 0.005;

    if abs(k1-k2)<tol && abs(k2-k3)<tol

        Phenotype = "Uniform";

    elseif (k3 < 0.01) && (k1 > k3) && (k2 > k3)

        Phenotype = "Plateau";

    elseif (k1 > k2) && (k2 > k3)

        Phenotype = "Decelerating";

    elseif (k1 < k2) && (k2 > k3)

        Phenotype = "Delayed";

    elseif (k1 < k2) && (k2 < k3)

        Phenotype = "Accelerating";

    elseif (k1 > k2) && (k3 > k2)

        Phenotype = "Reactivation";

    elseif abs(k1-k2)<tol && (k2 > k3)

        Phenotype = "Biphasic";

    else

        Phenotype = "Complex";

    end
    if abs(k1)>1e-6
        R1 = k2/k1;
    else
        R1 = NaN;
    end

    if abs(k2)>1e-6
        R2ratio = k3/k2;
    else
        R2ratio = NaN;
    end
    V_fit = exp(best_fit);

    Vmid = interp1(tt,V_fit,16,'linear','extrap');
    if Vmid < 0.5

        Regression = "Major";

    elseif Vmid < 0.7

        Regression = "Substantial";

    elseif Vmid < 0.9

        Regression = "Mild";

    else

        Regression = "Minimal";

    end
    % Default values
    ReplanFx = NaN;
    IdealReplanFx = NaN;
    IdealReplan = "NO";

    %% ---------------------------------------------------------
    % 1. Find first WHOLE treatment fraction where V_fit < 0.7
    % ----------------------------------------------------------

    % Whole treatment fractions
    wholeFx = 1:ceil(max(tt));

    % Evaluate fitted curve at whole treatment fractions
    % using shape-preserving PCHIP interpolation
    V_wholeFx = interp1(tt, V_fit, wholeFx, 'pchip');

    % First whole fraction where volume falls below 70%
    idx_replan = find(V_wholeFx < 0.7, 1, 'first');

    if ~isempty(idx_replan)
        ReplanFx = wholeFx(idx_replan);
    end


    %% ---------------------------------------------------------
    % 2. REPLAN CANDIDATE
    % Vmid < 0.7
    % ----------------------------------------------------------

    if Vmid < 0.7
        ReplanCandidate = "YES";
    else
        ReplanCandidate = "NO";
    end


    %% ---------------------------------------------------------
    % 3. IDEAL REPLAN
    % Vmid < 0.7 AND ReplanFx is between Fx 4 and Fx 20
    % ----------------------------------------------------------

    if Vmid < 0.7 && ...
            ~isnan(ReplanFx) && ...
            ReplanFx >= 4 && ReplanFx <= 20

        IdealReplan = "YES";
        IdealReplanFx = ReplanFx;

    end


    %% =========================================================
    % CREATE RESULTS TABLE ROW
    %% =========================================================

    PatientID = patientName;   % anonymous ID (P1, P2, ...) - no MR number

    NewRow = table( ...
        string(PatientID),...
        V0,...
        best_bp1,...
        best_bp2,...
        k1,...
        k2,...
        k3,...
        R1,...
        R2ratio,...
        k12,...
        k23,...
        Delta_k,...
        Dynamicity,...
        Phenotype,...
        Regression,...
        Vmid,...
        ReplanCandidate,...
        ReplanFx,...
        IdealReplanFx,...
        IdealReplan,...
        RMSE,...
        R2,...
        'VariableNames',{ ...
        'Patient',...
        'V0',...
        'BP1',...
        'BP2',...
        'k1',...
        'k2',...
        'k3',...
        'R1',...
        'R2',...
        'k12',...
        'k23',...
        'Delta_k',...
        'Dynamicity',...
        'Phenotype',...
        'Regression',...
        'Vmid',...
        'ReplanCandidate',...
        'ReplanFx',...
        'IdealReplanFx',...
        'IdealReplan',...
        'RMSE',...
        'R2_fit'});

    Results = [Results; NewRow];
    %% ==========================================
    % 4 PATIENTS PER PAGE
    %% ==========================================

    if mod(p-1,4)==0

        figure('Color','w',...
               'Position',[100 100 1200 800]);

        figCount = figCount + 1;

    end

    subplot(2,2,mod(p-1,4)+1)

    plot(tt,lnV_s,'b','LineWidth',2);
    hold on

    plot(tt,best_fit,'r--','LineWidth',2);

    xline(best_bp1,...
          '--k',...
          sprintf('BP1=%.1f',best_bp1));

    xline(best_bp2,...
          '--k',...
          sprintf('BP2=%.1f',best_bp2));

    grid on

    title(PatientID,'Interpreter','none')

    xlabel('Fraction')
    ylabel('ln(V/V0)')

    txt = sprintf([...
    'k1=%.3f\n' ...
    'k2=%.3f\n' ...
    'k3=%.3f\n' ...
   '\\Delta k=%.3f\n' ...
    'RMSE=%.4f\n' ...
    'R^2=%.4f\n\n' ...
    '%s'],...
    best_k1,...
    best_k2,...
    best_k3,...
    Delta_k,...
    RMSE,...
    R2,...
    Phenotype);
    text(0.58,...
         0.95,...
         txt,...
         'Units','normalized',...
         'VerticalAlignment','top',...
         'BackgroundColor','white',...
         'EdgeColor','black',...
         'FontSize',8);


    %% Save every 4 patients

    if mod(p,4)==0 || p==nPatients

        sgtitle(sprintf('Patients %d - %d',...
                        max(1,p-3),p));

        exportgraphics(gcf,...
            sprintf('PatientPage_%d.pdf',figCount),...
            'ContentType','vector');

    end
end

disp(' ');
disp('======================================');
disp('FINAL RESULTS');
disp('======================================');
disp(Results);
writetable(Results,'Tumour_Kinetics_Results.xlsx');

fprintf('\nResults saved as:\n');
fprintf('Tumour_Kinetics_Results.xlsx\n');
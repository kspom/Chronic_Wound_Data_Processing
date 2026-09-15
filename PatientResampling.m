%% PatientResampling: patient-level bootstrap, permutation test and leave-one-patient-out
% PLOS Medicine revision, September 2026 (Academic Editor, point 2)
%
% Input:  Patient_transitions.xlsx  (written by markov_chain_model.m; S15 Table):
%         per patient, outcome group, number of samples in each wound state and
%         number of transitions between consecutive weekly samples for each pair
%         of states (Imp = Impairment, Prol = Proliferation, Infl = Inflammation).
% Output: Resampling_results.xlsx  (S16 Table), three sheets:
%         Transitions      - transition probabilities per group with patient-level
%                            bootstrap 95% CI and leave-one-patient-out range
%         GroupDifferences - selected quantities, non-healers minus healers, with
%                            bootstrap 95% CI, leave-one-patient-out range and
%                            permutation p-value
%         Settings         - number of resamples, seed, additional summaries
%
% Cluster membership and wound-state assignment are held fixed at the values of
% the main analysis; only the set of patients (bootstrap, leave-one-out) or the
% outcome labels of the patients (permutation) are varied. The patient, not the
% sample, is the unit of resampling.

NBOOT = 1000;     % bootstrap resamples (stratified by outcome group)
NPERM = 10000;    % label permutations
rng(1);           % for reproducibility

%% Read per-patient counts
P = readtable('Patient_transitions.xlsx');
isH = strcmp(P.group,'Healer');
N = P{:, {'n_Imp','n_Prol','n_Infl'}};                       % patients x 3
stateNames = {'Imp','Prol','Infl'};
tn = cell(1,9); c = 0;
for a = 1:3, for b = 1:3, c = c+1; tn{c} = [stateNames{a} '_to_' stateNames{b}]; end, end
T = P{:, tn};                                                 % patients x 9 (from-major order)
np = size(P,1); nH = sum(isH); nN = sum(~isH);
fprintf('%d patients: %d healers, %d non-healers\n', np, nH, nN);

% quantities of interest: 9 transition probabilities + 3 occupancies per group (12 x 2),
% and 4 group differences (non-healers minus healers)
iInflImp = 7;  iProlProl = 5;  iImpInfl = 3;                  % column indices in T
qNames = {'P(Infl -> Imp)', 'P(Prol -> Prol)', 'P(Imp -> Infl)', 'Impairment occupancy'};

%% Observed values
[qH, prH, ocH] = summarise(T(isH,:),  N(isH,:),  iInflImp, iProlProl, iImpInfl);
[qN, prN, ocN] = summarise(T(~isH,:), N(~isH,:), iInflImp, iProlProl, iImpInfl);
dObs = qN - qH;

%% Bootstrap (resample patients with replacement within each group)
prBH = nan(NBOOT,9); prBN = nan(NBOOT,9); ocBH = nan(NBOOT,3); ocBN = nan(NBOOT,3); dB = nan(NBOOT,4);
hIdx = find(isH); nIdx = find(~isH);
for b = 1:NBOOT
    sH = hIdx(randi(nH, nH, 1));
    sN = nIdx(randi(nN, nN, 1));
    [q1, prBH(b,:), ocBH(b,:)] = summarise(T(sH,:), N(sH,:), iInflImp, iProlProl, iImpInfl);
    [q2, prBN(b,:), ocBN(b,:)] = summarise(T(sN,:), N(sN,:), iInflImp, iProlProl, iImpInfl);
    dB(b,:) = q2 - q1;
end

%% Permutation test (shuffle outcome labels among patients, group sizes fixed)
dP = nan(NPERM,4);
for k = 1:NPERM
    lab = false(np,1); lab(randperm(np, nH)) = true;
    q1 = summarise(T(lab,:),  N(lab,:),  iInflImp, iProlProl, iImpInfl);
    q2 = summarise(T(~lab,:), N(~lab,:), iInflImp, iProlProl, iImpInfl);
    dP(k,:) = q2 - q1;
end
pPerm = (1 + sum(abs(dP) >= abs(dObs), 1)) ./ (NPERM + 1);   % two-sided

%% Leave-one-patient-out
prLH = nan(np,9); prLN = nan(np,9); ocLH = nan(np,3); ocLN = nan(np,3); dL = nan(np,4);
for i = 1:np
    keep = true(np,1); keep(i) = false;
    [q1, prLH(i,:), ocLH(i,:)] = summarise(T(keep & isH,:),  N(keep & isH,:),  iInflImp, iProlProl, iImpInfl);
    [q2, prLN(i,:), ocLN(i,:)] = summarise(T(keep & ~isH,:), N(keep & ~isH,:), iInflImp, iProlProl, iImpInfl);
    dL(i,:) = q2 - q1;
end

%% Sheet 1: transition probabilities with bootstrap CI and LOPO range
rows = cell(18, 10); r = 0;
for g = 1:2
    if g == 1, gname = 'Healer';     Tg = T(isH,:);  pr = prH; prB = prBH; prL = prLH;
    else,      gname = 'Non-healer'; Tg = T(~isH,:); pr = prN; prB = prBN; prL = prLN; end
    Tsum = sum(Tg, 1);
    for a = 1:3
        tot = sum(Tsum((a-1)*3 + (1:3)));
        for b2 = 1:3
            c = (a-1)*3 + b2; r = r + 1;
            rows(r,:) = {gname, stateNames{a}, stateNames{b2}, Tsum(c), tot, pr(c), ...
                pct(prB(:,c), 2.5), pct(prB(:,c), 97.5), min(prL(:,c)), max(prL(:,c))};
        end
    end
end
Ttrans = cell2table(rows, 'VariableNames', {'Group','From_state','To_state','Transition_count', ...
    'Total_from_state','Probability','Bootstrap_CI95_lower','Bootstrap_CI95_upper','LOPO_min','LOPO_max'});

%% Sheet 2: group differences
Tdiff = table(qNames', qH', qN', dObs', pct(dB, 2.5)', pct(dB, 97.5)', min(dL,[],1)', max(dL,[],1)', pPerm', ...
    'VariableNames', {'Quantity','Healers','Non_healers','Difference_NH_minus_H', ...
    'Bootstrap_CI95_lower','Bootstrap_CI95_upper','LOPO_min','LOPO_max','Permutation_p'});

%% Sheet 3: settings and additional summaries
zeroFrac = mean(prBH(:, iInflImp) == 0);
Tset = table({'Bootstrap resamples'; 'Permutations'; 'Random seed'; ...
              'Bootstrap samples with healer P(Infl -> Imp) = 0 (fraction)'; ...
              'Leave-one-out runs'}, ...
             [NBOOT; NPERM; 1; zeroFrac; np], 'VariableNames', {'Setting','Value'});

%% Write
fname = 'Resampling_results.xlsx';
if isfile(fname), delete(fname); end
writetable(Ttrans, fname, 'Sheet', 'Transitions');
writetable(Tdiff,  fname, 'Sheet', 'GroupDifferences');
writetable(Tset,   fname, 'Sheet', 'Settings');
disp(Tdiff)
fprintf('Healer P(Infl -> Imp) = 0 in %.1f%% of bootstrap samples\n', 100*zeroFrac);
fprintf('Results written to %s\n', fname);

%% Local functions
function [q, pr, oc] = summarise(Tg, Ng, iInflImp, iProlProl, iImpInfl)
% transition probabilities (9, from-major order) and state occupancy (3) of a group,
% plus the four quantities of interest
    Tsum = sum(Tg, 1);
    pr = nan(1,9);
    for a = 1:3
        tot = sum(Tsum((a-1)*3 + (1:3)));
        if tot > 0, pr((a-1)*3 + (1:3)) = Tsum((a-1)*3 + (1:3)) / tot; end
    end
    Nsum = sum(Ng, 1);
    oc = Nsum / sum(Nsum);
    q = [pr(iInflImp), pr(iProlProl), pr(iImpInfl), oc(1)];
end

function v = pct(x, p)
% percentile of each column of x, ignoring NaN (linear interpolation, as in prctile)
    v = nan(1, size(x,2));
    for j = 1:size(x,2)
        s = sort(x(~isnan(x(:,j)), j));
        n = numel(s);
        if n == 0, continue; end
        pos = p/100 * n + 0.5;
        lo = floor(pos); hi = ceil(pos);
        lo = min(max(lo,1),n); hi = min(max(hi,1),n);
        v(j) = s(lo) + (s(hi) - s(lo)) * (pos - floor(pos));
    end
end

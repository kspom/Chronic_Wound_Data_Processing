%% RemoveLowRIN: sensitivity analysis - exclusion of low-RIN samples
% PLOS Medicine revision, September 2026 (Academic Editor, point 2)
%
% Requires in the current folder:
%   ChronicWoundData.mat           - output of Read_Data.m
%   S1_Table_Sample_Quality.xlsx   - S1 Table (per-sample RIN)
%
% Creates ChronicWoundData_lowRIN.mat: the same variables as in
% ChronicWoundData.mat, restricted to samples with RIN >= RINcut. Samples whose
% RIN could not be measured (empty cells in S1 Table; 17 Batch 2 samples) are
% removed as well. The gene filter of Read_Data.m is not repeated, i.e. the
% gene set is the same as in the main analysis.
%
% Usage:
%   RINcut = 4; RemoveLowRIN
%   then switch the load line in BatchEffectCorrection.m to
%   ChronicWoundData_lowRIN.mat and run BatchEffectCorrection,
%   markov_chain_model and SaveRINresults as usual. Z-scoring, cluster
%   construction and state assignment are all recomputed on the retained
%   samples.

if ~exist('RINcut','var')
    RINcut = 0;
end

load ChronicWoundData.mat

%% RIN per sample from S1 Table
Q = readtable('S1_Table_Sample_Quality.xlsx','PreserveVariableNames',true);
[tf,loc] = ismember(Metadata_filename, Q.sample);
if ~all(tf)
    error('%d samples of the metadata are not found in S1 Table', sum(~tf));
end
RIN = Q.RIN;
if iscell(RIN)
    RIN = str2double(RIN);      % in case empty cells were read as text
end
RIN = RIN(loc);

%% Select samples
keep = ~isnan(RIN) & RIN >= RINcut;

% all sample-level variables saved by Read_Data.m
sampleVars = {'Metadata_filename','patient','week_num','healed_label', ...
    'treatment_label','batch_label','healed_by_15_label','healed_visit_label'};
for i = 1:numel(sampleVars)
    if exist(sampleVars{i},'var')
        tmp = eval(sampleVars{i});
        eval([sampleVars{i} ' = tmp(keep);']);
    end
end
data = data(:,keep);
RIN  = RIN(keep);

%% Report
fprintf('RINcut = %g: %d of %d samples retained\n', RINcut, sum(keep), numel(keep));
fprintf('   Batch 1: %d,  Batch 2: %d\n', sum(batch_label==1), sum(batch_label==2));
fprintf('   healer samples: %d,  non-healer samples: %d\n', ...
    sum(healed_label==1), sum(healed_label==0));
up = unique(patient);
fprintf('   patients: %d', numel(up));
nps = cellfun(@(x) sum(strcmp(patient,x)), up);
fprintf('  (with <2 samples: %d)\n', sum(nps<2));

clear Q tf loc keep i tmp sampleVars up nps
save('ChronicWoundData_lowRIN.mat')

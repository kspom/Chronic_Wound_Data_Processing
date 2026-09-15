%% NormalizeCounts: sensitivity analysis - library-size normalization
% PLOS Medicine revision, September 2026 (Academic Editor, point 2)
%
% Requires in the current folder:
%   ChronicWoundData.mat   - output of Read_Data.m (raw Salmon read counts)
%
% Creates ChronicWoundData_normalized.mat: the same variables as in
% ChronicWoundData.mat, with counts divided by a per-sample size factor
% estimated by the median-of-ratios method (Anders & Huber 2010; the same
% estimator used inside rnaseqde and DESeq2). Each gene of a sample is divided
% by the same factor, so that differences in sequencing depth between samples
% are removed before log-transformation and batch correction.
%
% Usage:
%   NormalizeCounts
%   then switch the load line in BatchEffectCorrection.m to
%   ChronicWoundData_normalized.mat and run BatchEffectCorrection and
%   markov_chain_model as usual.

load ChronicWoundData.mat

%% Median-of-ratios size factors
ok = all(data > 0, 2);                    % genes with non-zero counts in every sample
lg = log(data(ok,:));                     % log counts of these genes
sf = exp(median(lg - mean(lg,2), 1));     % per-sample size factor (1 x samples)
data = data ./ sf;

%% Report
fprintf('Size factors estimated from %d genes without zero counts\n', sum(ok));
fprintf('   Batch 1: median %.2f, range %.2f - %.2f\n', median(sf(batch_label==1)), min(sf(batch_label==1)), max(sf(batch_label==1)));
fprintf('   Batch 2: median %.2f, range %.2f - %.2f\n', median(sf(batch_label==2)), min(sf(batch_label==2)), max(sf(batch_label==2)));

clear ok lg
save('ChronicWoundData_normalized.mat')

% goes through all cell-odor pairs, and calculates whether there is a
% response using:
% 1 t-test comparison for individual breaths; these reports tonic firing changes
% 2 KS comparison of PSTHs; this reports changes in phasic behaviour
%  differ from v2 in that we now group together  breaths 2-4
% skip v4 to align with MP_CompareResponseAfterLight5
% differ from v3 in structure of stats
% differ from v5 by breaking PCA into its own subfunction, and having flexible trial #s
% v7 gets a lookup table for trial #s
% v8 adds postodor ANOVA tests
% v9 adds better handling of empty breaths
% v10 adds ability to delete bad trials
function exp_out = MP_CalcPercentResponsiveCells10(exp_in, p_test )
if nargin < 2
    p_test = 0.01;
end

short_odor_threshold = 0.15;
exp_in = IdentifyShortOdorTrials(exp_in, short_odor_threshold); % get rid of trials with breaths before odor stabilizes
exp_in = CalcPostOdorSpkTimes(exp_in);
exp_in = MP_KSWholeExp8(exp_in, p_test);% this calculates all of the Kolmogorov Smirnov changes
exp_in = test_tonic_response(exp_in); % also does PCA
exp_in = PCA_response(exp_in);

exp_out = MP_ReportStatsResponse2(exp_in);


% calculate average # of spikes in a breath before odor, and during 4
% breaths of odor
function exp_out = test_tonic_response(exp_in)

exp_out = exp_in;
num_odors = length(exp_in.odor);
num_cells = length(exp_in.cells);
breath_length = exp_in.odor(1).breath_window / 1000;
breath_bins = -5*breath_length:breath_length:4*breath_length;
num_breath = length(breath_bins);

for o=1:num_odors
    trials = GetTrialsForExp(exp_in, o);
    num_trials = length(trials);
    spknum_breath = [];
    po_spknum_breath = [];
    for c=1:num_cells
        for t=1:num_trials
            cur_trial = trials(t);
            % get the individual trial average firing rates
           cur_spktimes = exp_in.odor(o).odor_SpkTimes_aligned{cur_trial,c};
           if isempty(cur_spktimes) % if cur_spktimes is all zeros, histc will return empty vector; reassign to zeros
               spknum_breath(o,c,t,:) = zeros(1,1,1,num_breath);
           else
               spknum_breath(o,c,t,:) = histc(cur_spktimes, breath_bins); % get number of spikes in each breath bin
           end

           % repeat for postodor
           po_spktimes = exp_in.odor(o).postodor_SpkTimes_aligned{cur_trial,c};
           if isempty(po_spktimes) % if cur_spktimes is all zeros, histc will return empty vector; reassign to zeros
               po_spknum_breath(o,c,t,:) = zeros(1,1,1,num_breath);
           else
               po_spknum_breath(o,c,t,:) = histc(po_spktimes, breath_bins); % get number of spikes in each breath bin
           end
        end
 
        %ANOVA!!!!!!!!!!!
        % breath 5 before odor is "ctl pre breath", then 4 baseline breaths, then 4 odor responsive breaths
        % for post-odor, use breaths 3-5
%         anova_labels = [repmat({'pre_ctl'},num_trials,1), repmat({'baseline'}, num_trials,4),...
%                         repmat({'br1'}, num_trials,1), repmat({'br2'}, num_trials,1), repmat({'br3-4'}, num_trials,2)];
        anova_labels = [repmat({'pre_ctl'},num_trials,1), repmat({'baseline'}, num_trials,4),...
                        repmat({'br1'}, num_trials,1), repmat({'br2'}, num_trials,1), repmat({'br3-4'}, num_trials,2),...
                        repmat({'brPO3-5'},num_trials,3)];
        anova_values = [squeeze(spknum_breath(o,c,:,1:9)) squeeze(po_spknum_breath(o,c,:,7:9))]; % use 1:9 for breaths because histc creates an extra bin
        [~, ~, stats] = anova1( anova_values(:), anova_labels(:), 'off' );
        shittyMultCompareOut = multcompare(stats);
        % multcompare outputs an odd pairwise comparison between groups;
        % the below line strips out the 4 breathgn cycles of odor + pre_ctl negative control
        % see multcompare help, and view shittyMultCompareOut to see it
        pos_check = shittyMultCompareOut([1 6:7],3)>0; % check for increases of firing rate after odor
        neg_check = shittyMultCompareOut([1 6:7],5)<0; % check for decreases in firing rate after odor
%         decreases of firing rate after odor
        ANOVA_hvalue(c,1:3,o) = (pos_check | neg_check)';
        br12_tonic_decrease(c,1:2, o) = neg_check(2:3);
        br12_tonic_increase(c,1:2, o) = pos_check(2:3);
        ANOVA_br12_hvalue(c,o) = (shittyMultCompareOut(10, 3) > 0) | (shittyMultCompareOut(10,5)<0);
        po_decrease(c,o) = shittyMultCompareOut(9,3)>0;
        po_increase(c,o) = shittyMultCompareOut(9,5)<0;
    end
end

exp_out.stats.ANOVA_preL_odor_response_hvalue = ANOVA_hvalue;
exp_out.stats.breaths.ANOVA_br1_dif_br2 = ANOVA_br12_hvalue;
exp_out.stats.postodor.po_ANOVA_increase = po_increase;
exp_out.stats.postodor.po_ANOVA_decrease = po_decrease;
exp_out.stats.br12_tonic_increase = br12_tonic_increase;
exp_out.stats.br12_tonic_decrease = br12_tonic_decrease;

%% PCA stuff
% this is untested since I pulled it out of v5
function exp_out = PCA_response(exp_in)

exp_out = exp_in;
num_odors = length(exp_in.odor);
num_cells = length(exp_in.cells);
breath_length = exp_in.odor(1).breath_window / 1000;
pre_breath = 5; post_breath = 5;
tot_breath = pre_breath + post_breath;
       PCA_bins = -pre_breath*breath_length:breath_length/8:post_breath*breath_length;
for o=1:num_odors
    trials = GetTrialsForExp(exp_in, o);
    for c=1:num_cells
       % subdivide into bins for classification algorithm
       temp_odor = SumOdorTrials(exp_in.odor(o), trials);
       cur_spktimes_sum = temp_odor.odor_SpkTimes_aligned_sum{c};
       hist_out = histc(cur_spktimes_sum, PCA_bins);
       hist_out2 = reshape(hist_out(1:tot_breath * 8), [8 tot_breath]);
       breath_odor_indices = (o*tot_breath-tot_breath+1):(o*tot_breath);
       cell_bin_indices = (c*8-7):(c*8);
       spksum_cellbin_breathodor(cell_bin_indices,breath_odor_indices) = hist_out2;
       
       % also create matrix for PCA time cycles - not updated to 5 post-odor breaths
       breathbin_odor = (o*72-71):(o*72);
       spksum_cell_breathbinodor(c, breathbin_odor) = hist_out(1:72)'; % matrix cell x (breath*bin*odor)
    end
end

% test distance between first breath and other breaths
spkdistance_matrix = dist(spksum_cellbin_breathodor);
[~, PCA_cellbin_breathodor, latents] = princomp(spksum_cellbin_breathodor'); % also try distance on PCA (when'm I gonna do a MANOVA?!?!)
PCA_distance_matrix = dist(PCA_cellbin_breathodor(:,1:5)');

% perform PCA with each cell as variable, and breath*phase*odor as observations
ctl_breath_indices = [1:5 10:14 19:23 28:32 37:41];
[~,PCA_scores, b] = princomp(spksum_cell_breathbinodor');

exp_out.stats.spksum_cellbin_breathodor = spksum_cellbin_breathodor; % for spkdistance
exp_out.stats.spkdistance_matrix = spkdistance_matrix;
exp_out.stats.PCA_cellbin_breathodor = PCA_cellbin_breathodor; % for PCA_distance
exp_out.stats.PCA_distance_matrix = PCA_distance_matrix;
exp_out.stats.spksum_cell_breathbinodor = spksum_cell_breathbinodor; % for PCA timecourse
exp_out.stats.wholeexp_PCA_scores = PCA_scores;
exp_out.stats.PCA_latents = latents;

% called by Script_CompareLowConc
% calculates distance between breaths for a single odor
% v2 takes the odor_index as input; v3 takes the odor stub as input
function out_dist_matrix = BreathSpikeDistanceBoot3(PCA_flag, odor_stub, varargin)
% list of odors to look at, e.g. 'Aa_', 'Eb_'
% varargin is a cell with the experiments 

num_bins = 8; % number of breath bins
max_PCA_coef = 16; % number of PCs to use
post_breath = 5; % number of breaths before and after odor onset
pre_breath = 5;
% reset these variables
cell_counter = 0;

% build average vector and individual trial vectors
for x = 1:length(varargin) % for each experiment
    % get current experiment
    cur_exp = varargin{x};
    trials = GetTrials(cur_exp.exp_info.exp_name); % get trials for the exp
    odor_ind = MP_GetOdorIndexFromStub(cur_exp, odor_stub);
    cur_odor = SumOdorTrials(cur_exp.odor(odor_ind), trials);
    cells = GetCellsForExp(cur_exp); % get cells with good data during the trials
    
    num_cells = length(cells);
    
    % set up PCA bins for the breath length of this exp
    breath_length = cur_odor.breath_window / 1000; %;CalcBreath(cur_exp, trials);
    PCA_bins = 0:breath_length/num_bins:post_breath*breath_length;
    ctl_bins = -pre_breath*breath_length:breath_length/num_bins:0;
    
    % get all individual trials for each odor
    for c_ind = 1:num_cells
        c = cells(c_ind); % cell index for accessing good cells in exp structure
        cell_counter = cell_counter + 1;
        
        % for 8 bins, should be 1:8, 9:16, etc
        cellbin_index = (cell_counter-1)*num_bins+1:cell_counter*num_bins;
            
        % build population avg vector for PCA
        sum_spktimes = cur_odor.odor_SpkTimes_aligned_sum{c};
        hist_out = histc(sum_spktimes, PCA_bins);
        hist_out2 = reshape(hist_out(1:post_breath * num_bins), [num_bins post_breath]);
        spksum_matrix(cellbin_index, 1:post_breath) = hist_out2 / length(trials);
        
        % get ctl average
        hist_out = histc(sum_spktimes, ctl_bins);
        hist_out2 = reshape(hist_out(1:pre_breath * num_bins), [num_bins pre_breath]);
        ctlsum_matrix(cellbin_index, 1:pre_breath) = hist_out2 / length(trials);
                    
    end % end cells
end % end experiments

if PCA_flag
    % now transform to PCA space, then 
    [PCA_coefs, PCA_spksum_matrix, latents] = princomp([ctlsum_matrix, spksum_matrix]');
    PCA_coefs = PCA_coefs(:,1:max_PCA_coef);

    weightings = repmat((latents(1:max_PCA_coef) / sum(latents))', size(PCA_spksum_matrix,1), 1);
    weighted_PCA_spksum_matrix = PCA_spksum_matrix(:,1:max_PCA_coef) .* weightings;
    
    % calculate spk distance
    spkdist_matrix = dist(weighted_PCA_spksum_matrix');
else
    spkdist_matrix = dist([ctlsum_matrix, spksum_matrix]);
end

spkdist_matrix = triu(spkdist_matrix); % delete bottom half (probably doesn't matter, but why not?)

% calculate mean distance between control breaths
ctl_breaths = spkdist_matrix(1:pre_breath, 2:pre_breath);
ctl_breaths = ctl_breaths(ctl_breaths>0);
mean_ctl_breath = mean(ctl_breaths);

% calculate mean distance between ctl breaths and odor breaths
for i=pre_breath+1:pre_breath+post_breath
   ctl_v_odor(i-pre_breath) = mean(spkdist_matrix(1:pre_breath,i));
end
tot_breaths = pre_breath + post_breath;
% calculate mean distance between breaths
   odor_v_odor = spkdist_matrix((pre_breath+1):(tot_breaths-1),(pre_breath+2):(tot_breaths))';
   odor_v_odor = odor_v_odor(odor_v_odor>0);
   
   out_dist_matrix = [mean_ctl_breath, ctl_v_odor, odor_v_odor'];
   
   fprintf([' Num_cells = ' num2str(cell_counter) '\n']);

%% in subset of experiments, animal is awake at beginning, then put to
% sleep; for those exp's, use trials 13-24; elsewise use 1-12
function trials_out2 = GetTrials(exp_name)
    trials12list = '111220'; % for first odor presentation of exp't, use trials 13-24
    trials15list = '130515A, 130515B, 130513A, 130513B, 130513C';
    if strfind(trials12list, exp_name)
        trials_out = [1 12];
    elseif strfind(trials15list, exp_name )
        trials_out = [3, 15];
    else
        trials_out = [1 10]; % for 2nd+ odor presentations, use trials 1-12
    end
    num_trials = 12;
    trials_out2 = randi(trials_out, 1, num_trials );
%     trials_out2 = 2:10; 

% in subset of experiments, animal is awake at beginning, then put to
% sleep; lose some cells between wake/sleep
function cells_out2 = GetCellsForExp(cur_exp)
exp_name = cur_exp.exp_info.exp_name;

    cells_out = 1:length(cur_exp.cells);
    num_cells = length(cells_out);
    cells_out2 = cells_out( randi([1 num_cells], 1, randi( [round(num_cells*2/3), num_cells],1) )); % random number of cells, averaging 2/3rds of them
%     cells_out2 = cells_out( randi([1 num_cells], 1, ceil(num_cells *rand(1)) ) ); % random number of cells, averaging half of them
%     cells_out2 = cells_out( randi([1 num_cells], 1, ceil(num_cells * rand(1)*3/4) ) ); % random number of cells, averaging 1/3rd of them

% troubleshooting code
% if exp_name == '130515A'
%     cells_out2 = 8;
% end
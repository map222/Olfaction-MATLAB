% called by MP_PredictOdorAnesthAllExp
% adapted from MP_PredictOdorTTX
% v2 adds post-odor breaths
% v4 adds pre-odor prediction and weighted PCA
% v5 adds option to use PCA or not
% v6 uses precalculated post-odor spktimes
% v7 gets passed the trials
% v8 is has more flexible odor handling
function [all_odor_pred cell_counter] = MP_PredictOdorBreath8(PCA_flag, num_bins, trials, varargin)
% list of odors to look at, e.g. 'Aa_', 'Eb_'
% varargin is a cell with the experiments 

max_PCA_coef = 25; % number of PCs to use
pre_breath = 5; post_breath = 4; % number of breaths before and after odor onset
post_odor_breath = 10;
tot_breath = pre_breath + post_breath;
% reset these variables
cell_counter = 0;
spksum_cellbin_breathodor = []; spknum_cellbin_breathodor = [];

odor_list = 1:3;
num_odor = length(odor_list);
% build average vector and individual trial vectors
for x = 1:length(varargin) % for each experiment
    % get current experiment
    cur_exp = varargin{x};
    cells = GetCellsForExp(cur_exp); % get cells with good data during the trials
    
    num_cells = length(cells);
    
    % set up PCA bins for the breath length of this exp
    breath_length = cur_exp.odor(1).breath_window / 1000;
    PCA_bins = -pre_breath*breath_length:breath_length/num_bins:post_breath*breath_length;
    po_PCA_bins = 0:breath_length/num_bins:post_odor_breath*breath_length;

    % get all individual trials for each odor
    for c_ind = 1:num_cells
        c = cells(c_ind); % cell index for accessing good cells in exp structure
        cell_counter = cell_counter + 1;
        
        % for 8 bins, should be 1:8, 9:16, etc
        cellbin_index = ((cell_counter-1)*num_bins+1):cell_counter*num_bins;
        for o = 1:num_odor 
            % get current odor
            odor_ind = odor_list(o);
            % dear god the next line was important
            cur_odor = SumOdorTrials(cur_exp.odor(odor_ind), trials);
            cur_odor = SumPostOdorTrials(cur_odor, trials);
            
            % for 10 breaths, should be 1:10, 11:20, etc
            breathodor_indices = (o*tot_breath-tot_breath+1):(o*tot_breath);
            
            % build population avg vector for PCA
            sum_spktimes = cur_odor.odor_SpkTimes_aligned_sum{c};
            hist_out = histc(sum_spktimes, PCA_bins);
            hist_out2 = reshape(hist_out(1:tot_breath * num_bins), [num_bins tot_breath]);
            spksum_cellbin_breathodor(cellbin_index, breathodor_indices) = hist_out2 / length(trials);
            
            for t=1:length(trials) % all trials
                % get spikes in this trial
                cur_spktimes = cur_odor.odor_SpkTimes_aligned{trials(t),c};
                hist_out = histc(cur_spktimes, PCA_bins);
                if isempty(hist_out)
                    hist_out = zeros(tot_breath* num_bins,1);
                end
                hist_out2 = reshape(hist_out(1:tot_breath * num_bins), [num_bins tot_breath]);
                spknum_cellbin_breathodor(cellbin_index, breathodor_indices, t) = hist_out2;
                
                % get average vector minus the test trial
                temp_trials = trials(trials~=trials(t));
%                 temp_trials = trials_postTTX;
                temp_odor = SumOdorTrials(cur_odor, temp_trials);
                extrial_spktimes = temp_odor.odor_SpkTimes_aligned_sum{c};
                hist_out = histc(extrial_spktimes, PCA_bins);
                hist_out2 = reshape(hist_out(1:tot_breath * num_bins), [num_bins tot_breath]);
                spksum_extrial(cellbin_index, breathodor_indices,t) = hist_out2 / length(temp_trials);
            end % end trials
            
            %% repeat above for post-odor response
            
            po_breathodor_indices = (o*post_odor_breath-post_odor_breath+1):(o*post_odor_breath);
            for t=1:length(trials) % all trials
                cur_spktimes = cur_odor.postodor_SpkTimes_aligned{trials(t),c};
                hist_out = histc(cur_spktimes, po_PCA_bins);
                if isempty(hist_out)
                    hist_out = zeros(post_odor_breath* num_bins,1);
                end
                hist_out2 = reshape(hist_out(1:post_odor_breath * num_bins), [num_bins post_odor_breath]);
                po_spknum_cellbin_breathodor(cellbin_index, po_breathodor_indices, t) = hist_out2;
                
                % get average vector minues the test trial
                temp_trials = trials(trials~=trials(t));
%                 temp_trials = trials_postTTX;
                temp_odor = SumPostOdorTrials(cur_odor, temp_trials);
                extrial_spktimes = temp_odor.postodor_SpkTimes_aligned_sum{c};
                hist_out = histc(extrial_spktimes, po_PCA_bins);
                hist_out2 = reshape(hist_out(1:post_odor_breath * num_bins), [num_bins post_odor_breath]);
                po_spksum_extrial(cellbin_index, po_breathodor_indices,t) = hist_out2 / length(temp_trials);
            end
            
            % build population avg vector for PCA
            temp_odor = SumPostOdorTrials(cur_odor, trials);
            sum_spktimes = temp_odor.postodor_SpkTimes_aligned_sum{c};
            hist_out = histc(sum_spktimes, po_PCA_bins);
            hist_out2 = reshape(hist_out(1:post_odor_breath * num_bins), [num_bins post_odor_breath]);
            po_spksum_cellbin_breathodor(cellbin_index, po_breathodor_indices) = hist_out2 / length(trials);
            
        end % end odors
    end % end cells
    
end % end experiments

if PCA_flag
    % now transform to PCA space, then 
    [PCA_coefs, PCA_cellbin_breathodor, latents] = princomp([spksum_cellbin_breathodor, po_spksum_cellbin_breathodor]');
    PCA_coefs = PCA_coefs(:,1:max_PCA_coef);

    weightings = repmat((latents(1:max_PCA_coef) / sum(latents))', size(PCA_cellbin_breathodor,1), 1);
    weighted_PCA_cellbin_breathodor = PCA_cellbin_breathodor(:,1:max_PCA_coef) .* weightings;
end

%% prediction time
% for each trial and odor-breath, make a prediction
for bo=1:tot_breath*num_odor % num_breathodors
   % get breathodor position
   allbreathodors = (0:(num_odor-1)) * tot_breath + mod(bo-1, tot_breath) +1;
    for t= 1:length(trials)
        cur_trial = spknum_cellbin_breathodor(:, bo, t);
        if PCA_flag
           % transform spikes into PCA space using coefs from average
    %        cur_trial = mean(squeeze(spknum_cellbin_breathodor(:, bo, :)), 2);
           PCA_cur_trial = cur_trial' * PCA_coefs;
           weighted_PCA_cur_trial = PCA_cur_trial .*weightings(1,:);

           % some error checking stuff
    %        a = dist([PCA_cur_trial', PCA_cellbin_breathodor(:,1:10)']); % this is for checking that breath is close to itself
    %        a = dist([cur_trial, spksum_cellbin_breathodor(:,1:10)]);
           %  mean(squeeze(po_spknum_cellbin_breathodor(:, 1, 1:10)), 2) calculate distance for all three odors

           % check the prediction
%             temp_dist = dist([PCA_cur_trial', PCA_cellbin_breathodor(allbreathodors,1:max_PCA_coef)']); % unweighted
            temp_dist = dist([weighted_PCA_cur_trial', weighted_PCA_cellbin_breathodor(allbreathodors,:)']);
            
        else % not PCA
            temp_dist = dist([cur_trial, spksum_extrial(:, allbreathodors,t)]);
        end
        
        [~, pred_odor(bo, t)] = min(temp_dist(1, 2:(num_odor+1)));
    end
end

%% post-odor prediction time

for bo=1:post_odor_breath*num_odor % num_breathodors
    % get breathodor position
   po_allbreathodors = (0:(num_odor-1)) * post_odor_breath + mod(bo-1, post_odor_breath) +num_odor * tot_breath + 1;
    for t= 1:length(trials)
       % transform spikes into PCA space using coefs from average
       po_cur_trial = po_spknum_cellbin_breathodor(:, bo, t);
       if PCA_flag
           
    %        cur_trial = mean(squeeze(spknum_cellbin_breathodor(:, bo, :)), 2);
           po_PCA_cur_trial = po_cur_trial' * PCA_coefs;
           weighted_po_PCA_cur_trial = po_PCA_cur_trial .*weightings(1,:);

           % some error checking stuff
    %        a = dist([PCA_cur_trial', PCA_cellbin_breathodor(:,1:10)']); % this is for checking that breath is close to itself
    %        a = dist([cur_trial, spksum_cellbin_breathodor(:,1:10)]);

            % calculate distance for all three odors
    %         temp_dist = dist([po_PCA_cur_trial', PCA_cellbin_breathodor(po_allbreathodors,1:max_PCA_coef)']);
            temp_dist = dist([weighted_po_PCA_cur_trial', weighted_PCA_cellbin_breathodor(po_allbreathodors,:)']);
        else % not PCA
            temp_dist = dist([po_cur_trial, po_spksum_extrial(:, po_allbreathodors-num_odor * tot_breath,t)]);
%             temp_dist = dist([po_cur_trial, po_spksum_cellbin_breathodor(:, po_allbreathodors-30)]);
        end
        [~, po_pred_odor(bo, t)] = min(temp_dist(1, 2:(num_odor+1)));
    end
end


%% arrange for output
for o = 1:num_odor
    all_odor_pred(o,1:tot_breath) = sum(pred_odor(tot_breath*(o-1) + 1: (tot_breath*o), :) == o,2) / length(trials);
    all_odor_pred(o,tot_breath+1:(tot_breath+1+post_odor_breath-1)) = sum(po_pred_odor(((o-1)*post_odor_breath+1):(o*post_odor_breath), :) == o,2) / length(trials);
end

% in subset of experiments, animal is awake at beginning, then put to
% sleep; for those exp's, use trials 13-24; elsewise use 1-12
function trials_out = GetTrials(exp_name)
    trials13list = '120118A, Exp120308A_02, Exp120308B_02, Exp120308C_02, Exp120305A_01, Exp120305B_02, Exp120305C_02, Exp120413A_01';
    if strfind(trials13list, exp_name)
        trials_out = 13:24;
    elseif strcmp(exp_name, '120113A')
        trials_out = 25:36;
    else
        trials_out = 1:12;
    end

% in subset of experiments, animal is awake at beginning, then put to
% sleep; lose some cells between wake/sleep
function cells_out = GetCellsForExp(cur_exp)
exp_name = cur_exp.exp_info.exp_name;

    if strcmp(exp_name, '120118A')
        cells_out = [1,2,4,5,7,9,10,12,13:21, 23:29];
%     elseif strcmp(exp_name, '120520C')
%         cells_out = [5];
%     elseif strcmp(exp_name, '120520D')
%         cells_out = [1,3,5,11,14,18,19, 24];
    else
        cells_out = 1:length(cur_exp.cells);
    end

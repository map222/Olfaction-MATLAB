function exp_out = ReadBehavCSV(filename)
if nargin <=1
   [filename, pathname]= uigetfile('*.csv'); 
   filename = [pathname filename];
end

% get starts of protocols
exp_out.behaviour_session = get_protocol_indices(filename);
num_protocols = length(exp_out.behaviour_session);
for p = 1:num_protocols
    cur_info = exp_out.behaviour_session(p);
    csv_in = dlmread(filename, ';', ['B', num2str(cur_info.start_index), '..R', num2str(cur_info.end_index)]);
    csv_in2 = csv_in(1:2:length(csv_in), :); % skip every other because it has an odd time component

    % separate the input into trials
    num_trials = max(csv_in2(:, 2));
    for t = 1:num_trials % go through all trials

        data_trials{t} = csv_in2(csv_in2(:,2)==t, :);
        data_trials2(t) = process_trial(data_trials{t});

    end

    exp_out.behaviour_session(p).data = data_trials2;
end

function struct_out = process_trial( cur_trial)
    % use find below due to diff offset
    lick_starts = find(diff(cur_trial(:,6))==1)+1; % index for start of each lick
    struct_out.lick_starts = cur_trial(lick_starts, 3); % turn index into time
    
    lick_ends = find(diff(cur_trial(:,6))==-1)+1; % time for end of each lick
    struct_out.lick_ends = cur_trial(lick_ends, 3);
    
    H20_times = find(diff(cur_trial(:,9))==1)+1; % time for start of water reward
    struct_out.H20_times = cur_trial(H20_times, 3);
    
    % get row/column for protocol starts
    % actually read in protocol as well
function protocol_out = get_protocol_indices(filename)
    % read in the csv file, with all data as strings
    fstuff = fopen(filename);
    csv_in = textscan(fstuff, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'Delimiter', ';', 'EndOfLine', '\n');
    first_column = csv_in{1};
    
    % there is no find function for cells, so step through all values to
    % find the protocol starts
    start_indices = [];
    for i=1:length(first_column)
        if strcmp(first_column{i}, 'Time') % lines in file that have protocol info
            start_indices = [start_indices, i];
        end
    end
    start_indices = start_indices +1; % add one for offset to actual start
    num_protocols = length(start_indices);
        
    % get end indices
    end_indices = zeros(1, num_protocols);
    for i=1: (num_protocols-1)
        end_indices(i) = start_indices(i+1) - 6;
    end
    end_indices(num_protocols) = length(first_column);
    
    % now identify the protocol
    fourth_column = csv_in{4}; % this has "seqence" label
    for p=1:num_protocols
        protocol_labels{p} = fourth_column{start_indices(p) -3};
    end
    
    % refactor data structure
    for i=1:num_protocols
       protocol_out(i).start_index = start_indices(i); 
       protocol_out(i).end_index = end_indices(i);
       protocol_out(i).protocol_label = protocol_labels{i};
    end    
    
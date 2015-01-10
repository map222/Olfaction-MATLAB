
% this function compares two given breaths against each other
function [h_value p_value] = MP_KSTwoOdorBreaths( odor_SpkTimes_in, breath_length, breaths_in, plot_flag, p_test)
% odor_SpkTimes_in is aligned_sum SpkTimes from ExpYYMMDD.mat file, for a specific odor
% breath length is the breath length in seconds
% breaths_in is which breaths you want to test (e.g. [1 2])
% plot flag determines whether to plot the distributions
% p_test is the p value you want to be under (sent as alpha to kstest2)

if nargin < 4 % default to no plotting
    plot_flag = 0;
end
if nargin < 5 % default to no plotting
    p_test = 0.05;
end

breath1 = breaths_in(1);
breath2 = breaths_in(2);
% get baseline spike distribution from TWO breaths before odor to avoid any
% odor leakage
nth_breath_spk_times = odor_SpkTimes_in( odor_SpkTimes_in > (breath_length * (breath1 -1)) & odor_SpkTimes_in < (breath_length * breath1) );
nth_breath_spk_times = nth_breath_spk_times - breath_length * (breath1 -1); % move them to same time frame
breath2_odor_spk_times = odor_SpkTimes_in( odor_SpkTimes_in > (breath_length * (breath2 -1)) & odor_SpkTimes_in < (breath_length * breath2) );
breath2_odor_spk_times = breath2_odor_spk_times - breath_length * (breath2 -1); % move them to same time frame

% if you have a very low firing rate, some breaths may have no spikes
if (isempty(breath2_odor_spk_times) || isempty(nth_breath_spk_times) )
    h_value = 0; p_value = 1;
    fprintf('One cell had no spikes during a breath');
    return;
end
[h_value, p_value] = kstest2(breath2_odor_spk_times, nth_breath_spk_times, p_test);

if plot_flag

    figure; hold all;
    cdfplot(breath2_odor_spk_times);
    cdfplot(nth_breath_spk_times);
end
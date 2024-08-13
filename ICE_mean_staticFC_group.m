% OCT 22, 2023
% compute the mean of all subjects'static FC in each ICE group
% -------------------------------------------------------------------------
% Folder with subj directories
% NORM 22 HC
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/2023_NORM_HOCo_output';

% ICE- Bitemporal
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/HOCo_outputs_entire_brain/Bitemporal';
% ICE- Left-TLE
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/HOCo_outputs_entire_brain/Left_TLE';
% ICE- Right-TLE
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/HOCo_outputs_entire_brain/Right_TLE';
% ICE- Entire-TLE
par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/HOCo_outputs_entire_brain/Entire_TLE';


subjs = dir(par_dir);
dir_flag = [subjs.isdir];
subjs(~dir_flag) = [];
subjs(1:2) = [];
NofROIs = 56;
% 

% =========================================================================
% Run across subjects
% -------------------------------------------------------------------------
group_static_FC = zeros(NofROIs,NofROIs,length(subjs));
for i=1:length(subjs) % across subjects
    subjName = subjs(i).name;
    subj_FC = dir(fullfile(par_dir,subjName,'CorrMtx', 'FC_*'));
    subj_FC_orig  = load(fullfile(par_dir,subjName,'CorrMtx', subj_FC.name));
%     fullfile(subj_FC_orig.FC_orig.StaticFCArray)
    subj_static_FC = subj_FC_orig.FC_orig.StaticFCArray;
    % replace nan with zero in the cell array
%     subj_static_FC = cellfun(@(M) subsasgn(M, substruct('()', {isnan(M)}), 0), subj_static_FC, 'uniform',0);
    true_empties = cellfun('isempty', subj_static_FC);
%     subj_static_FC(true_empties) = {0};     % replace [] with zero in the cell array
    subj_static_FC(true_empties) = {nan};     % replace [] with zero in the cell array
    subj_static_FC = cell2mat(subj_static_FC(2:end, 2:end));
    % concatenate all the static FC from all the subject in a 3D array
    group_static_FC(:,:,i) = subj_static_FC;
end

mean_group_static_FC = mean(group_static_FC,3);
% to sort the ROIs in networks
ROIs_network_order = [2,3,6,7,11,17,20,30,35,10,38,48,22,24,27,29,32,37,51,14,15,19,28,31,36,44,50,53,56,...
5,8,9,23,25,33,34,41,42,45,46,47,49,54,55,1,4,12,13,16,18,21,26,39,40,43,52];

sorted_mean_group_static_FC = zeros(NofROIs);
sorted_mean_group_static_FC= mean_group_static_FC(ROIs_network_order,:);
sorted_mean_group_static_FC= sorted_mean_group_static_FC(:,ROIs_network_order);
% figure;imagesc(z_trans_sorted_networks); colorbar
save('sorted_mean_group_static_FC','sorted_mean_group_static_FC')

figure; colormap('Parula'); colorbar
h = imagesc(sorted_mean_group_static_FC); title('sorted-networks-Entire-TLE-mean-static-FC'); colorbar
% make the diagonal white/transparent
set(h, 'AlphaData', 1-isnan(sorted_mean_group_static_FC))

% draw horizontal lines to separate the five networks
network_split=[9.5,12.5,19.5,29.5,44.5,56.5]; % #of ROIs for each of the five networks: 9,12,19,29,44,56
yline(network_split, '-',LineWidth=2)
yticks([])
% draw vertical lines to separate the five networks
xline(network_split, '-',LineWidth=2)
xticks([])
axis square



            
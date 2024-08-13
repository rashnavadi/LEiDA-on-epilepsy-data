% Oct 23, 2023, Tahereh
% ICE groups timeseries for LEiDA
% -------------------------------------------------------------------------
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence'
addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/for_LEiDA/phase_coherence'

% 200 TR
% Folder with subj'TS directories
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/
% 200TR--Bitemporal
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/subjects_TSs/Bitemporal_TS';
% 200TR--Left_TLE
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/subjects_TSs/Left_TLE_TS';
% 200TR--Right_TLE
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/subjects_TSs/Right_TLE_TS';
% 200TR--Entire TLE (42 subjects, 12 Bitemporal + 17 Left-TLE + 13 Right-TLE)
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/subjects_TSs/Entire_TLE_TS';

% 150 TR
% Folder with subj'TS directories
addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR'

% 150TR--Bitemporal
par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/subjects_TSs/Bitemporal_TS';
% 150TR--Left_TLE
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/subjects_TSs/Left_TLE_TS';
% 150TR--Right_TLE
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/subjects_TSs/Right_TLE_TS';
% 150TR--Entire TLE (42 subjects, 12 Bitemporal + 17 Left-TLE + 13 Right-TLE)
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/subjects_TSs/Entire_TLE_TS';



cd(fullfile(par_dir))
subjs = dir([par_dir, '/*.mat']);
NofSubjs = numel(subjs);
allSubjsNames = {subjs.name}';
Group = ones(NofSubjs,1);

% TS_length = 200;
TS_length = 150;
NofROIs   = 56;

% =========================================================================
% Run across subjects
% -------------------------------------------------------------------------
%  ORDER THE ROIS BASED ON THE NETWORKS MIST64 (64-8=56 ROIS)
ROIs_network_order = [2,3,6,7,11,17,20,30,35,10,38,48,22,24,27,29,32,37,51,14,15,19,28,31,36,44,50,53,56,...
5,8,9,23,25,33,34,41,42,45,46,47,49,54,55,1,4,12,13,16,18,21,26,39,40,43,52]; % ROIs for each network
sorted_networks = zeros(TS_length,56);

% MIST-20 (20-2 cerebellum=18)
% ROIs_network_order = [3,5,9,2,10,8,14,4,12,13,16,6,15,18,1,7,11,17];
% MIST-12 (12-1cerebellum=11)
% ROIs_network_order = [3,7,5,4,6,8,9,10,1,2,11];
% MIST-7 (7-1cerebellum=6)
% ROIs_network_order = [1,2,3,4,5,6];

Subjects = cell(size(allSubjsNames,1),1);
for i=1:NofSubjs % across subjects
%     Subjects = subjs.name{i,1}(1:end-4);
    load(cell2mat(fullfile(allSubjsNames(i)))) % after loading
    aa = TS(2,:);
    subj_TS = cell2mat(aa);
    sorted_subj_TS = subj_TS(:,ROIs_network_order);
%     sorted_subj_TS = TS(:,ROIs_network_order); % sort columns based on the ROIs'networks%     
    subj_name_i = allSubjsNames{i,1}(1:end-4);
    eval(sprintf('%s=sorted_subj_TS', subj_name_i))
    Subjects (i,1)= {subj_name_i};
%     save(sprintf('%s_sorted_ROIs',subj_name_i), 'sorted_subj_TS'); % it saves each of the file as a .mat file
end

% save -append sorted_networks_Bitemporal_for_LEiDA.mat Group Subjects  
  
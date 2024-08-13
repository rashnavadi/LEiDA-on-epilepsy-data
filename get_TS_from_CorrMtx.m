% NORM timeseries for LEiDA
% -------------------------------------------------------------------------
% Folder with subj directories
addpath '/Users/trashnavadi/Documents/Data_Analysis/dFC_HOCo_2020/from_Brad/Normdata'
% par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/MIST_20_HOCo_output';
par_dir = '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/MIST_64_HOCo_output';

subjs = dir(par_dir);
% dir_flag = [subjs.isdir];
% subjs(~dir_flag) = [];
subjs(1:3) = [];

stringToRemove = '_std'; % replace with the string you want to remove from folder names
allSubjsNames = {subjs.name}';
subjs = subjs(~contains(allSubjsNames, stringToRemove));
Subjects = cell(size(subjs,1), 1);
for i=1:size(subjs,1)
    Subjects(i) = {subjs(i).name};
end
Group = ones(size(subjs,1),1);

% =========================================================================
% Run across subjects
% -------------------------------------------------------------------------
% IF YOU WANT TO ORDER THE ROIS BASED ON THE NETWORKS MIST64 (64-8=56 ROIS)
ROIs_network_order = [2,3,6,7,11,17,20,30,35,10,38,48,22,24,27,29,32,37,51,14,15,19,28,31,36,44,50,53,56,...
5,8,9,23,25,33,34,41,42,45,46,47,49,54,55,1,4,12,13,16,18,21,26,39,40,43,52]; % ROIs for each network
sorted_networks = zeros(150,56);

% MIST-20 (20-2 cerebellum=18)
% ROIs_network_order = [3,5,9,2,10,8,14,4,12,13,16,6,15,18,1,7,11,17];
% MIST-12 (12-1cerebellum=11)
% ROIs_network_order = [3,7,5,4,6,8,9,10,1,2,11];
% MIST-7 (7-1cerebellum=6)
% ROIs_network_order = [1,2,3,4,5,6];


for i=1:size(subjs,1) % across subjects
    func_file = dir(fullfile(par_dir,subjs(i).name));
    load(fullfile(par_dir,subjs(i).name,'CorrMtx', 'Time_series.mat'))
    aa = TSs(2,:);
    subj_TS = cell2mat(aa);
    sorted_subj_TS = subj_TS(:,ROIs_network_order); % sort columns based on the ROIs' networks
    eval(sprintf('%s=sorted_subj_TS', subjs(i).name))

%     save NORM_Subjs_for_LEiDA.mat subj_TS
end

% save -append NORM_Subjs_for_LEiDA.mat Group Subjects  
  
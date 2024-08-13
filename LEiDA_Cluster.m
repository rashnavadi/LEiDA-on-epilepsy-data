function LEiDA_Cluster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function computes a k-means clustering algorith to the eigenvectors
%
% - Loads all the eigenvectors saved in LEiDA_data
% - Computes k-means clustering
% - Calculates the Dunn's score
% - Selects the optimal solution
%
% - Saves the optimal solution into LEiDA_Clusters.mat
%   Clusters =
%      IDX: Cluster index for each observation
%        C: Cluster centroids
%        D: Distance from each observation to every centroid
%     SUMD: Within-cluster sums of point-to-centroid distances
%

% modified by Tahereh Rashnavadi, 2024
% Written by 
% Joana Cabral joana.cabral@psych.ox.ac.uk
% Paulo Marques paulo.c.g.marques@gmail.com
% Last edited May 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/for_LEiDA/'))
% CHANGE THE PATH FOR THE CD
% NORM, HC
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/for_LEiDA/MIST_20
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/for_LEiDA/MIST_64
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/for_LEiDA

% Bitemporal
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Bitemporal/2024/
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/for_LEiDA/phase_coherence/Bitemporal/

% Left-TLE
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Left_TLE/2024
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Left_TLE

% Right-TLE
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Right_TLE
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Right_TLE/2024

% Entire-TLE
% cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Entire_TLE
cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Entire_TLE/2024



load LEiDA_data Leading_Eig

N_sub=size(Leading_Eig,1);

% Generate vector X concatenating all eigenvectors from all subjects and
% time points, where rows are observations and collumns are variables.
% load(fullfile('/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/finding_k/Left_TLE/X.mat'))
X=[];
for s=1:N_sub
    X=cat(1,X,squeeze(Leading_Eig(s,:,:))); 
end
% clear Leading_Eig
% N_sub= 17; % Left-TLE
% N_sub= 12; % 
% N_sub=13;

%Kmeans clustering
maxk=15;
% opt= statset('UseParallel',1); %,'UseSubstreams',1);
% The options may vary according to the Matlab version
Kmeans_results=cell(1,15);

parfor k=1:maxk  
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','MaxIter',50,'Replicates',20,'Display','final'); %,'Options',opt);   
    Kmeans_results{k}.IDX=IDX;
    Kmeans_results{k}.C=C; 
    Kmeans_results{k}.SUMD=SUMD; 
    Kmeans_results{k}.D=D; 
    % added by Tahereh, OCT 2023
    Elbow_opt_nC(k) = mean(Kmeans_results{k}.SUMD);
    Sil_opt_nC(k)   = mean(silhouette(X, Kmeans_results{k}.IDX));
    Kmeans_results{k}.Elbow(k) = Elbow_opt_nC(k);
    Kmeans_results{k}.Sils(k)  = Sil_opt_nC(k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate Clustering performance using Dunn's score
distM_fcd=squareform(pdist(X,'cityblock'));
dunn_score=zeros(maxk,1);
for j=2:maxk
    dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(j) 'clusters'])
end
[~,ind_max]=max(dunn_score);
disp(['Best clustering solution: ' num2str(ind_max) ' clusters']);

Clusters= Kmeans_results{ind_max};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best k based on results from silhouette
[~,ind_max_sils]=max(Sil_opt_nC);
Clusters_sils = Kmeans_results{ind_max_sils};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dunn's score and Dave-Bouldin
Dav_Bld    = evalclusters(X,'kmeans','DaviesBouldin','klist',1:maxk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('LEiDA_Clusters','Clusters', 'Clusters_sils', 'Elbow_opt_nC', 'Sil_opt_nC', 'Dav_Bld','dunn_score','X')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting Dunn's score and Dave-Bouldin
figure, clf
subplot(221)
plot(dunn_score, 'k^-', 'markerfacecolor', 'k', 'markersize', 15)
set(gca, 'xtick', 1:maxk, 'xlim', [0.5 maxk+1])
title('Dunn score method')
xlabel('Number of brain states')

subplot(222)
plot(Dav_Bld.CriterionValues, 'k^-', 'markerfacecolor', 'k', 'markersize', 15)
set(gca, 'xtick', 1:maxk, 'xlim', [0.5 maxk+1])
title('Davies-Bouldin method')
xlabel('Number of brain states')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting elbow and silhouette
% figure, clf
subplot(223)
plot(Elbow_opt_nC, 'k^-', 'markerfacecolor', 'k', 'markersize', 15)
set(gca, 'xtick', 1:maxk, 'xlim', [0.5 maxk+1])
title('Elbow method')
xlabel('Number of brain states')
% 
subplot(224)
plot(Sil_opt_nC, 'k^-', 'markerfacecolor', 'k', 'markersize', 15)
set(gca, 'xtick', 1:maxk, 'xlim', [0.5 maxk+1])
title('Silhouette method')
xlabel('Number of brain states')
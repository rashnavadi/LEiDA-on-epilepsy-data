function LEiDA_analysis
% modified by Tahereh Rashnavadi, October 17, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function plots the cluster solutions on the cortex, and is matrix 
% format, in order to visualize the FC patterns.
% It alsso plots the probability of each cluster and the weighted sum of
% FC patterns in order to compare with the mean BOLD FC matrix
%
% Written by 
% Joana Cabral joana.cabral@psych.ox.ac.uk
% Last edited May 2016 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORM
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/for_LEiDA
% Bitmeporal
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Bitemporal/2024'
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Bitemporal'
% Left-TLE
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Left_TLE/2024'
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Left_TLE'
% Right-TLE
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Right_TLE/2024'
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Right_TLE'
% Entire-TLE
addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Entire_TLE/2024'
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Right_TLE'

% November 8th, 2023, finding the best k value when using LEiDA on HOCo FC matrices (not the phase coherencies), Tahereh
% Bitmeporal
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/for_LEiDA/phase_coherence/Bitemporal'
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Bitemporal'
% Left-TLE
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/for_LEiDA/phase_coherence/Left_TLE'
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/finding_k/Left_TLE
% Right-TLE
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/for_LEiDA/phase_coherence/Right_TLE'
% addpath '/Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Right_TLE'



load LEiDA_Clusters Clusters
[N_Cl, N_ba]=size(Clusters.C);
h=hist(Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Clusters.C(ind,:);
VVT_mean=zeros(N_ba);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this reordering is not needed as i have already reorder the ROIs for
% MIST56
% % To reorder matrix plots
% Order=[1:2:N_ba N_ba:-2:2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% colormap(jet)    
colormap(parula)    

% Pannel A
% Plot the cluster centroids over the cortex 
% Pannel C 
% Plot the centroids outerproduct
% normalize the centroids

for c=1:N_Cl
%     subplot(2,N_Cl+1,c)
%     plot_nodes_in_cortex(V(c,:))
%     title(['#' num2str(c)])
%     axis off
%     subplot(2,N_Cl+1,c+N_Cl+1)
    subplot(2,N_Cl,c)
    VVT=V(c,:)'*V(c,:);   
%     imagesc(VVT(Order,Order)); 
    imagesc(VVT); 
    clim([prctile(VVT(:),2) prctile(VVT(:),98)]); colorbar
    axis square
    title('Outer product') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    VVT_mean=VVT_mean+squeeze(VVT)*y(c);
    network_split=[9.5,12.5,19.5,29.5,44.5,56.5]; % #of ROIs for each of the five networks: 9,12,19,29,44,56
    yline(network_split, '-',LineWidth=2)
%     yticks([])
    % draw vertical lines to separate the five networks
    xline(network_split, '-',LineWidth=2)
%     xticks([])
end


VVT_mean=VVT_mean/sum(y);

% Pannel B
% Plot the probabilities of each cluster
% The Colormap is adjusted for 5 clusters 
% figure;
% mymap=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1];
mymap=[1 0 0; 1 1 0; 0 1 0];

x=y/sum(y);
% subplot(2,N_Cl+1,N_Cl+1)
subplot(2,N_Cl,N_Cl+1)
hold on
for c=1:N_Cl
    bar(c,x(c),'Facecolor',mymap(c,:))
end
set(gca,'xtick',1:N_Cl)
box off
xlabel('State Number')
ylabel('Probability')
xlim([0 6])

% Pannel D plot the weighted sum of outer products
subplot(2,N_Cl,N_Cl+2)
% imagesc(VVT_mean(Order,Order))
imagesc(VVT_mean)
clim([prctile(VVT_mean(:),2) prctile(VVT_mean(:),98)]); colorbar
title({'Weighted sum of VV^T'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square
y=[9.5,12.5,19.5,29.5,44.5,56.5]; % #of ROIs for each of the five networks: 9,12,19,29,44,56
yline(y, '-',LineWidth=2)
% yticks([])
% draw vertical lines to separate the five networks
xline(y, '-',LineWidth=2)
% xticks([])


load sorted_mean_group_static_FC.mat
I_sup_diag=find(triu(ones(N_ba),1));
[cc p]=corrcoef(sorted_mean_group_static_FC(I_sup_diag),VVT_mean(I_sup_diag))

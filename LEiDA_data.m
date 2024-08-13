function LEiDA_data
% modified bt Tahereh Rashnavadi, October 17, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function processes the data for LEiDA
%
% - Loads the BOLD data from all subjects in a matrix of Time x Brain areas
% - Computes the BOLD phases via Hilbert transform
% - Calculates the iFC (or instantaneous BOLD coherence matrix)
% - Calculates the instantaneous Leading Eigenvector
% - Calculates the instantaneous % of variance
% - Calculates de time x time FCD matrices
%
% Saves into LEiDA_data.mat
% Leading_Eig - Leading Eigenvector at each timepoint & each subject 
% Var_Eig     - % of variance of the leading eigenvector 
% FCD_eig     - Cosine similarity of eigenvectors over time 
%
% Written by 
% Joana Cabral joana.cabral@psych.ox.ac.uk
% Last edited May 2016 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORM HC
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/Healthy_Controls/from_Brad_in2020/for_LEiDA

% Bitemporal
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Bitemporal

% Left-TLE
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Left_TLE/2024/
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Left_TLE/

% Right-TLE
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/HOCo_outputs_entire_brain/Right_TLE/2024/
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Right_TLE/


% Entire-TLE
addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/150TR/MIST_64/for_LEiDA/phase_coherence/Entire_TLE/2024/
% addpath /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs/HOCo_inputs/200TR/for_LEiDA/phase_coherence/Entire_TLE/


% Get the BOLD data Subjects names
load 'sorted_networks_Entire_TLE_for_LEiDA.mat' Subjects 
N_sub=length(Subjects);

% Set parameters and preallocate variables
N=56; % Number of MIST brain areas considered
% Tmax=150; % TR=2s
Tmax=150; % TR=1.5s
Phases=zeros(N,Tmax);
Leading_Eig=zeros(N_sub,Tmax,N);
Var_Eig=zeros(N_sub,Tmax);
FCD_eig =zeros(N_sub,Tmax-2,Tmax-2);

for s=1:N_sub
    % Load BOLD data and save it as "signaldata'
    disp(['  Subject ' num2str(s) ' from ' num2str(N_sub)])
    Subject=char(Subjects(s));
    load('sorted_networks_Entire_TLE_for_LEiDA.mat',Subject) %%%%%%%%******* CHANGE THIS LINE FOR EACH DATA ******%%%%%%%%
    eval(['signaldata = ' Subject '(:,1:N);']);
    eval(['clear ' Subject ';']);
    
    % Get the Phase of BOLD data
    for seed=1:N
        timeseriedata=demean(detrend(signaldata(:,seed)));
        Phases(seed,:) = angle(hilbert(timeseriedata));
    end
    
    % Get the iFC leading eigenvector at each time point 
    iFC=zeros(N);
    
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(n,p)=cos(adif(Phases(n,t),Phases(p,t)));
            end
        end
        [eVec,eigVal]=eig(iFC);
        eVal=diag(eigVal);
        [val1, i_vec_1] = max(eVal);
        Leading_Eig(s,t,:)=eVec(:,i_vec_1);
        % Calculate the variance explained by the leading eigenvector
        Var_Eig(s,t)=val1/sum(eVal);
    end
    
    % Compute the FCD for each subject
    for t=1:Tmax
        eig1=squeeze(Leading_Eig(s,t,:));
        for t2=1:Tmax
            eig2=squeeze(Leading_Eig(s,t2,:));
            FCD_eig(s,t,t2)=dot(eig1,eig2)/norm(eig1)/norm(eig2);
        end
    end
end

save('LEiDA_data','Leading_Eig','FCD_eig','Var_Eig')

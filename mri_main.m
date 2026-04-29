%% Reset all
clear all
close all
restoredefaultpath
addpath("C:\Github\WM_Distractors_OPM\Matlab")

%% Set up fieldtrip

A_Initialize_environment_OPM
global STATIC;
STATIC = StaticVarClass_WM_shift_OPM();
ft_default.showcallinfo = 'no';

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/21099_opm/';
    base_save_path = '/home/chrpfe/Documents/21099_opm/';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/21099_opm/phalanges';
    on_server = true;
else
    on_server = false;
    base_data_path = 'G:\SV10_OPM_Distractor_Processing\SV10_RawData\MEG\raw';
    base_save_path = 'G:\SV10_OPM_Distractor_Processing\ProcessedData_Natmeg';
    base_matlab_path = 'C:\Github\WM_Distractors_OPM\Matlab';
    project_scripts_path = 'C:\toolbox\opm_general';
end

params.src_density = '32'; % Sourcemodel density ('4', '8' or '32') = approximate number of sourmes per hemisphere


%% Subjects + dates
% subses = {'0005' '240208';
%     '0905' '240229';
%     '0916' '240320';
%     '0953' '241104';
%     '1096' '241022';
%     '1153' '240321';
%     '1167' '240425';
%     '1186' '240925';
%     '1190' '241023';
%     '1191' '241024';
%     '1193' '241029';
%     '1194' '241029';
%     '1195' '241030'};
% mri_files = {'00000001.dcm' 
%     '/mri/sub-15931_T1w.nii.gz'  
%     '/nifti/anat/sub-15985_T1w.nii.gz'};
% 
% if on_server
%     subs_to_run = 1:size(subses,1);
% else
%     subs_to_run = 2; %1:size(subses,1)
% end
% excl_subs = [1];

subs_to_run = 1:height(STATIC.subjectData);

%% Loop over subjects
for i_sub = subs_to_run
    params.sub = ['sub_' num2str(i_sub,'%02d')];

    % % Paths
    % raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    % mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    save_path_mri = fullfile(base_save_path,'MRI',params.sub);
    
    % Create folders if they do not exist yet
    if ~exist(fullfile(base_save_path,'MRI'), 'dir')
        mkdir(fullfile(base_save_path,'MRI'))
    end
    if ~exist(save_path_mri, 'dir')
        mkdir(save_path_mri)
    end
    meg_file = STATIC.subjectData.OPM_FilePath(i_sub);
    wsl_path_subj = char(fullfile(STATIC.wsl_path,"WM_Distractors_OPM",string(STATIC.subjectID(i_sub))));
    prepare_mri(wsl_path_subj,meg_file,save_path_mri,params);
    close all
end
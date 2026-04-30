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

base_data_path = 'G:\SV10_OPM_Distractor_Processing\SV10_RawData\MEG\raw';
base_save_path = 'G:\SV10_OPM_Distractor_Processing\ProcessedData_Natmeg';
base_matlab_path = 'C:\Github\WM_Distractors_OPM\Matlab';
project_scripts_path = 'C:\toolbox\opm_general';
params.src_density = '32'; % Sourcemodel density ('4', '8' or '32') = approximate number of sourmes per hemisphere

%% Loop over subjects
for i_sub = 1:height(STATIC.subjectData)
    params.sub = ['sub-' num2str(i_sub,'%02d')];

    save_path_mri = fullfile(base_save_path,params.sub ,'MRI');

    % Create folders if they do not exist yet
    if ~exist(save_path_mri, 'dir')
        mkdir(save_path_mri)
    end

    % meg_file is the file that contains the headshape points
    meg_file = STATIC.subjectData.MEGFilePath(i_sub);
    wsl_path_subj = char(fullfile(STATIC.wsl_path,"WM_Distractors_OPM",string(STATIC.subjectID(i_sub))));
    prepare_mri(wsl_path_subj,meg_file,save_path_mri,params);
    close all
end
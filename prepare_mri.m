function prepare_mri(mri_file,opm_file,meg_file,save_path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global STATIC;
    if ~exist(mri_file,'file')
        error(['Did not find MRI file: ' mri_file])
    end

    if ~exist(meg_file,'file')
        error(['Did not find MEG file with headshape: ' meg_file])
    end

    %% Read data
    headshape = ft_read_headshape(meg_file);
    grad    = ft_read_sens(opm_file,'senstype','meg'); % Load MEG sensors
    % elec    = ft_read_sens(meg_file,'senstype','eeg'); % Load EEG electrodes
    mri = ft_read_mri(mri_file);
    
    %% Align fiducials
    ft_sourceplot([], mri);
    mri_coordsys = ft_determine_coordsys(mri);
    cfg = [];
    cfg.method   = 'interactive';
    cfg.coordsys = 'neuromag';
    mri_realigned_1 = ft_volumerealign(cfg, mri_coordsys);

    %% ICP align
    cfg = [];
    cfg.method              = 'headshape';
    cfg.headshape.headshape = headshape;
    cfg.headshape.icp       = 'yes';
    cfg.headshape.interactive    = 'no';
    mri_realigned_2 = ft_volumerealign(cfg, mri_realigned_1);

    % Check co-registration
    cfg.headshape.icp       = 'no';        % Do not fit points again
    cfg.headshape.interactive    = 'yes';
    mri_realigned_2 = ft_volumerealign(cfg, mri_realigned_2);
    
    %% Reslice MRI
    %cfg = [];
    %cfg.resolution = 1;
    %mri_resliced = ft_volumereslice(cfg, mri_realigned_2);
    %mri_resliced = ft_convert_units(mri_resliced, 'cm');
    mri_resliced = ft_convert_units(mri_realigned_2, 'cm');
    
    save(fullfile(save_path, 'mri_resliced.mat'), 'mri_resliced'); disp('done')

    %% Segment MRI
    cfg = [];
    cfg.output = {'brain' 'skull' 'scalp'};
    mri_segmented = ft_volumesegment(cfg, mri_resliced);
    save(fullfile(save_path, 'mri_segmented.mat'), 'mri_segmented'); disp('done')
    
    %% Correct compartments
    binary_brain = mri_segmented.brain;
    binary_skull = mri_segmented.skull | binary_brain;
    binary_scalp = mri_segmented.scalp | binary_brain | binary_skull;
    
    % use boolean logic together with IMERODE
    binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
    binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull

    % Copy MRI
    mri_segmented_2 = mri_segmented;                    % Copy stucture
    mri_segmented_2.anatomy = mri_resliced.anatomy;  % Copy anatomical data
      
    % insert the updated binary volumes, taking out the center part for skull and scalp
    mri_segmented_2.brain    = binary_brain;
    mri_segmented_2.skull    = binary_skull & ~binary_brain;
    mri_segmented_2.scalp    = binary_scalp & ~binary_brain & ~binary_skull;

    %% Make meshes
    cfg = [];
    cfg.method = 'projectmesh';
    cfg.tissue = 'brain';
    cfg.numvertices = 3000;
    mesh_brain = ft_prepare_mesh(cfg, mri_segmented_2);
    
    cfg.tissue = 'skull';
    cfg.numvertices = 2000;
    mesh_skull = ft_prepare_mesh(cfg, mri_segmented_2);
    
    cfg.tissue = 'scalp';
    cfg.numvertices = 2000;
    mesh_scalp = ft_prepare_mesh(cfg, mri_segmented_2);
    
    % Collect meshes into a single structure
    meshes = [mesh_brain mesh_skull mesh_scalp];
    
    figure
    ft_plot_mesh(mesh_brain,'EdgeAlpha',0,'FaceAlpha',1,'FaceColor','r')
    ft_plot_mesh(mesh_skull,'EdgeAlpha',0,'FaceAlpha',0.5)
    ft_plot_mesh(mesh_scalp,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[229 194 152]/256)

    %% Headmodels
    cfg = [];
    cfg.method = 'singleshell';
    headmodel_meg = ft_prepare_headmodel(cfg, mesh_brain);
    
    try
        cfg = [];
        cfg.method = 'bemcp';
        cfg.conductivity = [1 1/20 1] .* (1/3);  % Standard values     
        headmodel_eeg = ft_prepare_headmodel(cfg, meshes);
    catch 
        headmodel_eeg = [];
    end
    %% Plot
    figure
    ft_plot_sens(grad)
    ft_plot_headshape(headshape)
    ft_plot_mesh(mesh_scalp,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[229 194 152]/256)
    ft_plot_headmodel(headmodel_meg)
    ft_plot_axes([], 'unit', 'cm','coordsys','neuromag');
   
    
    %%
    figure
    % ft_plot_sens(elec, 'style', 'ok','elecsize',10);
    ft_plot_headshape(headshape);
    ft_plot_headmodel(headmodel_eeg,'facealpha', 0.5,'FaceColor',[229 194 152]/256)
    ft_plot_axes([], 'unit', 'cm');

    %%
    %headshape = ft_read_headshape(aux_file);
    %hpi_polhemus = headshape.pos(find(contains(headshape.label,'hpi')),:);
    %headshape_tf = headshape;
    %headshape_tf.pos = opm_trans.transformPointsForward(headshape.pos);
    %hpi_polhemus_tf = headshape_tf.pos(find(contains(headshape.label,'hpi')),:);
    
    %meshes_opm = meshes;
    %meshes_opm(1).pos = opm_trans.transformPointsForward(meshes_opm(1).pos);
    %meshes_opm(2).pos = opm_trans.transformPointsForward(meshes_opm(2).pos);
    %meshes_opm(3).pos = opm_trans.transformPointsForward(meshes_opm(3).pos);

    %headmodel_opm = headmodel_meg;
    %headmodel_opm.bnd.pos = opm_trans.transformPointsForward(headmodel_opm.bnd.pos);

    %headmodel_opmeeg = headmodel_eeg;
    %headmodel_opmeeg.bnd(1).pos = opm_trans.transformPointsForward(headmodel_opmeeg.bnd(1).pos);
    %headmodel_opmeeg.bnd(2).pos = opm_trans.transformPointsForward(headmodel_opmeeg.bnd(2).pos);
    %headmodel_opmeeg.bnd(3).pos = opm_trans.transformPointsForward(headmodel_opmeeg.bnd(3).pos);

    %%
    headmodels = [];
    headmodels.headmodel_meg = headmodel_meg;
    %headmodels.headmodel_opm = headmodel_opm;
    headmodels.headmodel_eeg = headmodel_eeg;
    %headmodels.headmodel_opmeeg = headmodel_opmeeg;

    %meshes = [];
    %meshes.meg = meshes_meg;
    %meshes.opm = meshes_opm;

    %% Save
    save(fullfile(save_path, 'headmodels.mat'), 'headmodels');
    save(fullfile(save_path, 'meshes.mat'), 'meshes');

    %% Extract file to Freesurfer
    % mri_segmented_2_mm = ft_convert_units(mri_segmented_2,'mm');
    % 
    % cfg            = [];
    % cfg.resolution = 1;
    % cfg.dim        = [256 256 256];
    % mri_fs            = ft_volumereslice(cfg, mri_segmented_2_mm);
    % 
    % 
    % mri_fs = ft_convert_coordsys(mri_fs,'acpc');
    % 
    % 
    [parentDir, ~, ~] = fileparts(save_path);
    if ~exist(fullfile(parentDir, 'freesurfer'),'dir')
        mkdir(fullfile(parentDir, 'freesurfer'));
    end
    cfg              = [];
    cfg.filename     = fullfile(parentDir, 'freesurfer', 'sub02');
    cfg.filetype     = 'mgh'; 
    cfg.parameter    = 'anatomy';
    ft_volumewrite(cfg, mri_realigned_2);
end
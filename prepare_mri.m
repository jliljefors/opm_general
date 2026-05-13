function prepare_mri(mri_path,meg_file,save_path, params)
%Prepare_mri Prepare headmodel, sourcemodels and mri for MEG analysis
%   Requires the following inputs:
%   mri_path: path to freesurfer folder containing MRIs and workbench
%             cortical surface meshes
%   meg_file: path to a meg file containing polhemus points for
%             co-registration
%   save_path: path where to save the results

    if ~isfield(params,'src_density')
        params.src_density = '8'; % default source density to use. Results in ~16k sources
    end

    mri_file = fullfile(mri_path, 'mri', 'T1.mgz');
    if ~exist(mri_file,'file')
        error(['Did not find MRI file: ' mri_file])
    end

    if ~exist(meg_file,'file')
        error(['Did not find MEG file with headshape: ' meg_file])
    end

    %% Read data
    headshape = ft_read_headshape(meg_file);
    mri = ft_read_mri(mri_file);
mri.coordsys = 'ras';

mri = ft_volumereslice([], mri);     % the real reslice — fixes parity once and for all

    %% Prealign with fiducials
    cfg = [];
    cfg.method   = 'interactive';
    cfg.coordsys = 'neuromag';
    mri_realigned_1 = ft_volumerealign(cfg, mri);

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
    
    % cfg.tissue = 'skull';
    % cfg.numvertices = 2000;
    % mesh_skull = ft_prepare_mesh(cfg, mri_segmented_2);
    % 
    % cfg.tissue = 'scalp';
    % cfg.numvertices = 2000;
    % mesh_scalp = ft_prepare_mesh(cfg, mri_segmented_2);
    % 
    % % Collect meshes into a single structure
    % meshes = [mesh_brain mesh_skull mesh_scalp];
    % 
    % figure
    % ft_plot_mesh(mesh_brain,'EdgeAlpha',0,'FaceAlpha',1,'FaceColor','r')
    % ft_plot_mesh(mesh_skull,'EdgeAlpha',0,'FaceAlpha',0.5)
    % ft_plot_mesh(mesh_scalp,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[229 194 152]/256)

    %% Headmodels
    cfg = [];
    cfg.method = 'singleshell';
    headmodel_meg = ft_prepare_headmodel(cfg, mesh_brain);
    
    % try
    %     cfg = [];
    %     cfg.method = 'bemcp';
    %     cfg.conductivity = [1 1/20 1] .* (1/3);  % Standard values     
    %     headmodel_eeg = ft_prepare_headmodel(cfg, meshes);
    % catch 
    %     headmodel_eeg = [];
    % end

    %%
    headmodels = [];
    headmodels.headmodel_meg = headmodel_meg;
    % headmodels.headmodel_eeg = headmodel_eeg;

    %% Save
    save(fullfile(save_path, 'headmodels.mat'), 'headmodels');
    save(fullfile(save_path, 'mesh_brain.mat'), 'mesh_brain');

    %% Sourcemodel
    % Read and transform cortical restrained source model
    clear sourcemodel sourcemodel_inflated
    files = dir(fullfile(mri_path,'workbench'));
    for i = 1:length(files)
        if endsWith(files(i).name,['.L.midthickness.' params.src_density 'k_fs_LR.surf.gii'])
            filename = fullfile(mri_path,'workbench',files(i).name);
        elseif endsWith(files(i).name,['.L.aparc.' params.src_density 'k_fs_LR.label.gii'])
            filename2 = fullfile(mri_path,'workbench',files(i).name);
        end
    end
    sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
    % sourcemodel.coordsys = 'scanras';
    aparc_L = ft_read_atlas({filename2,filename});
    aparc_R = ft_read_atlas({strrep(filename2,'.L.','.R.'),strrep(filename,'.L.','.R.')});
    tmp = ft_read_atlas(strrep(filename2, '.L.', '.R.'),'format','caret_label');
    n_labels = length(aparc_L.parcellationlabel);
    atlas = [];
    atlas.parcellationlabel = [aparc_L.parcellationlabel; aparc_R.parcellationlabel];
    atlas.parcellation = [aparc_L.parcellation; aparc_R.parcellation + n_labels];
    atlas.rgba = [aparc_L.rgba; aparc_R.rgba; [0 0 0 1]];
    n_labels = length(atlas.parcellationlabel);
    atlas.parcellation(isnan(atlas.parcellation))=n_labels+1;
    sourcemodel.brainstructure = atlas.parcellation;
    sourcemodel.brainstructurelabel = atlas.parcellationlabel;
    sourcemodel.brainstructurecolor = atlas.rgba;
    clear atlas aparc_L aparc_R

    T = mri_resliced.transform/mri.transform;
    sourcemodel = ft_transform_geometry(T, sourcemodel);
    sourcemodel.inside = true(size(sourcemodel.pos,1),1);

     figure; hold on
    ft_plot_mesh(mesh_brain, 'EdgeAlpha', 0, 'FaceAlpha', 0.3, 'FaceColor', 'r')
    ft_plot_mesh(sourcemodel, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.6)
    title('Source model after ICP realignment to brain mesh')
    view(3); axis equal

    for i = 1:length(files)
        if endsWith(files(i).name,['.L.inflated.' params.src_density 'k_fs_LR.surf.gii'])
            filename = fullfile(mri_path,'workbench',files(i).name);
        end
    end
    sourcemodel_inflated = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
    sourcemodel_inflated = ft_transform_geometry(T, sourcemodel_inflated);    
    sourcemodel_inflated.inside = true(size(sourcemodel_inflated.pos,1),1);
    sourcemodel_inflated.brainstructure = sourcemodel.brainstructure;
    sourcemodel_inflated.brainstructurelabel = sourcemodel.brainstructurelabel;
    sourcemodel_inflated.brainstructurecolor = sourcemodel.brainstructurecolor;

    sourcemodel.unit = 'cm';
    sourcemodel_inflated.unit = 'cm';
    save(fullfile(save_path, 'sourcemodel'), 'sourcemodel', '-v7.3');
    save(fullfile(save_path, 'sourcemodel_inflated'), 'sourcemodel_inflated', '-v7.3');
end
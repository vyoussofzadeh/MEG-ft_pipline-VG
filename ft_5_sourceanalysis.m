% clc; clear; close all;
% 
% ft_defaults
% 
% sub = input('subject number (e.g = 1)?');
% p = ['sub',num2str(sub)];
% load(['.\data\',p]);

% Calling params
headmodel_time = output.headmodel.headmodel_time;
vol            = output.mri.vol;
timeAll        = output.timelockanalysis.Verbs_Data;
timeBaseline   = output.timelockanalysis.Verbs_post;
timePost       = output.timelockanalysis.Verbs_Baseline;
mri_aligned    = output.mri.mri_aligned;

% freqAll = output.freqanalysis.Verbs_Data;
% freqBaseline = output.freqanalysis.Verbs_post;
% freqPost = output.freqanalysis.Verbs_Baseline;


%%
% cfg                 = [];
% cfg.method          = 'lcmv';
% cfg.grid            = headmodel_time; % leadfield, which has the grid information
% cfg.vol             = vol; % volume conduction model (headmodel)
% % cfg.keepfilter      = 'yes';
%  cfg.keeptrials     = 'yes';
% cfg.lcmv.keepfilter  = 'yes';
% cfg.lcmv.fixedori   = 'yes'; % project on axis of most variance using SVD
% %  cfg.frequency = [12,23];
% % cfg.lcmv.projectnoise = 'yes';
% cfg.lcmv.lambda = '0.1%';
% sourceAll = ft_sourceanalysis(cfg, timeAll);

%% Bandpass filter before source localization
% use wide-band 4-30 Hz reconstruction
% Fbp = [1 12];
% N   = 1;
% 
% filteredData = ft_preproc_bandpassfilter(timeAll.avg, 1200, Fbp, N);
% timeAll.avg = filteredData;
% 
% filteredData = ft_preproc_bandpassfilter(timePost.avg, 1200, Fbp, N);
% timePost.avg = filteredData;
% 
% filteredData = ft_preproc_bandpassfilter(timeBaseline.avg, 1200, Fbp, N);
% timeBaseline.avg = filteredData;
%%

% create spatial filter using the lcmv beamformer
cfg                 = [];
cfg.method          = 'lcmv';
cfg.grid            = headmodel_time; % leadfield, which has the grid information
cfg.vol             = vol; % volume conduction model (headmodel)
cfg.keepfilter      = 'yes';
cfg.lcmv.keepfilter  = 'yes';
 cfg.keeptrials     = 'yes';
cfg.lcmv.fixedori   = 'yes'; % project on axis of most variance using SVD
%  cfg.frequency = [12,23];
% cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.lambda = '0.1%';
sourceAll = ft_sourceanalysis(cfg, timeAll);

% remove center of head bias by using the neural activity index (NAI)
% sourceAll.avg.pow = sourceAll.avg.pow ./sourceAll.avg.noise;

cfg.grid.filter = sourceAll.avg.filter;
cfg.grad = timeBaseline.grad;
sourceBaseline = ft_sourceanalysis(cfg, timeBaseline);
% sourceBaseline.avg.pow = sourceBaseline.avg.pow ./sourceBaseline.avg.noise;

cfg.grad = timePost.grad;
sourcePost = ft_sourceanalysis(cfg, timePost);
% sourcePost.avg.pow = sourcePost.avg.pow ./sourcePost.avg.noise;

%% subtraction ...  
cfg = [];
cfg.parameter = 'pow';
cfg.operation = '(x1-x2)/(x1+x2)';
S1bl = ft_math(cfg,sourcePost,sourceBaseline);

cfg = [];
cfg.parameter = 'pow';
source_diff_int = ft_sourceinterpolate(cfg, S1bl, mri_aligned);

cfg = [];
% cfg.template = 'E:\My Matlab\My codes\My GitHub\fieldtrip\external\spm8\templates\T1.nii';
% cfg.template = 'E:\My Matlab\My codes\My GitHub\My tool\MEG_ft_pipline\child template\nihpd_sym_04.5-08.5_pdw.nii';
% cfg.spmversion = 'spm8';
% cfg.template = 'single_subj_T1.nii';
normalised = ft_volumenormalise(cfg, source_diff_int);

%%
% cfg = [];
% cfg.parameter = 'pow';
% sourceBaseline_int = ft_sourceinterpolate(cfg, sourceBaseline, mri_aligned);
% sourcePost_int = ft_sourceinterpolate(cfg, sourcePost, mri_aligned);
% 
% source_diff_int = sourcePost_int;
% source_diff_int.pow = (sourcePost_int.pow - sourceBaseline_int.pow) ./ (sourcePost_int.pow + sourceBaseline_int.pow);

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs';
cfg.opacitymap = 'rampup';
cfg.funcolorlim = [-0.15 0.15];
cfg.opacitylim = [-0.15 0.15];
ft_sourceplot(cfg, normalised);
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
ft_sourceplot(cfg, normalised);
view ([-70 20 50])
% view ([70 20 50])
light ('Position',[-70 20 50])
colormap jet

cfg = [];
% cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs';
cfg.opacitymap = 'rampup';
cfg.funcolorlim = [-0.15 0.15];
cfg.opacitylim = [-0.15 0.15];
% ft_sourceplot(cfg, normalised);
cfg.method = 'surface';
cfg.camlight       = 'no';
cfg.projthresh     = 0.6;
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
ft_sourceplot(cfg, normalised);
view ([70 20 50])
light ('Position',[70 20 50])
colormap jet

cfg = [];
% cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs';
cfg.opacitymap = 'rampup';
cfg.funcolorlim = [-0.15 0.15];
cfg.opacitylim = [-0.15 0.15];
% ft_sourceplot(cfg, normalised);
cfg.method = 'surface';
cfg.camlight       = 'no';
cfg.projthresh     = 0.6;
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
ft_sourceplot(cfg, normalised);
view ([-70 20 50])
light ('Position',[-70 20 50])
colormap jet

% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = 'pow';
% cfg.filename  = ['.\data\sa_',p];
% ft_volumewrite(cfg, normalised)


%% saving data
output.source.sourceAll       = sourceAll;
output.source.sourcePost      = sourcePost;
output.source.sourceBaseline  = sourceBaseline;
output.source.source_diff_int = source_diff_int;

% save(['.\data\',p], 'output');

% %% Laterality index
% Verticies  = S1bl.pos(S1bl.inside,:);
% 
% vol1 = ft_convert_units(vol, 'mm');
% 
% figure
% hold on
% ft_plot_vol(vol1, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% % ft_plot_mesh(Verticies);
% %
% % left hem
% idx_l = find(Verticies(:,2)<0);
% idx_r = find(Verticies(:,2)>0);
% ft_plot_mesh(Verticies(idx_l,:),'facecolor','b');
% ft_plot_mesh(Verticies(idx_r,:),'facecolor','b');
% 
% %%
% pow = S1bl.pow(S1bl.inside,:);
% pow_l = pow(idx_l);
% pow_r = pow(idx_r);
% 
% Lat = (mean(pow_l) - mean(pow_r)) / (mean(pow_l) + mean(pow_r))

% figure,
% plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'b*');
% box off
% set(gca,'color','none');
% axis off




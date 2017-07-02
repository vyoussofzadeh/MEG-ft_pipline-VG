% clc; clear; close all;
% 
% ft_defaults
% 
% sub = input('subject number (e.g = 1)?');
% p = ['sub',num2str(sub)];
% load(['.\data\',p]);

timeAll = output.timelockanalysis.Verbs_Data;
vol = output.mri.vol;


cfg = [];
cfg.grad = timeAll.grad;
cfg.headmodel = vol;
cfg.reducerank = 2;
cfg.normalize = 'yes';
cfg.normalizeparam = 1;
cfg.channel = {'MEG', '-MLC12'};
cfg.grid.resolution = 10;
cfg.grid.unit = 'mm';
headmodel_time = ft_prepare_leadfield(cfg);

% saving data
output.headmodel.headmodel_time     = headmodel_time;

%%
% construct the dipole grid in the template brain coordinates
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.resolution = 1;
cfg.grid.tight  = 'yes';
cfg.headmodel   = vol;
cfg.reducerank = 2;
cfg.channel = {'MEG', '-MLC12'};
cfg.normalize = 'yes';
cfg.grad = timeAll.grad;
sourcemodel   = ft_prepare_sourcemodel(cfg);
sourcemodel = ft_convert_units(sourcemodel, 'mm');


figure
hold on
ft_plot_vol(vol, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));

cfg             = [];
cfg.grid        = sourcemodel;
cfg.headmodel   = vol;
cfg.channel     = {'MEG'};
cfg.grad        = timeAll.grad;
cfg.normalize = 'yes';
cfg.normalizeparam = 1;
sourcemodel_lf  = ft_prepare_leadfield(cfg, timeAll);


% save(['.\data\',p], 'output');

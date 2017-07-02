% clc; clear; close all;
% 
% ft_defaults
% 
% sub = input('subject number (e.g = 1)?');
% p = ['sub',num2str(sub)];
% load(['.\data\',p]); 

% cd ..
% cd ..
% mri = ft_read_mri(['.\T1 scans\CTL',num2str(sub),'\T1.nii']);
mri = ft_read_mri([mripath,'\T1.nii']);
dataset = output.dataset;


cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
mri_aligned = ft_volumerealign(cfg, mri);

% mri_resliced = ft_volumereslice([], mri_aligned);
% cfg = [];
% mri_aligned = ft_volumenormalise(cfg, mri_aligned);

cfg = [];
cfg.output = {'brain', 'skull', 'scalp'};
segmentedmri = ft_volumesegment(cfg, mri_aligned);
disp(segmentedmri);

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

vol = ft_convert_units(vol, 'mm');
sens = ft_read_sens(dataset);
figure
ft_plot_sens(sens, 'style', '*b');
hold on
ft_plot_vol(vol);

% mri_resliced = ft_volumereslice([], segmentedmri);

%% saving data
% output.mri.mri_resliced     = mri_resliced;
output.mri.vol     = vol;
output.mri.mri_aligned = mri_aligned;

% save(['.\Time domain\Noun condition\data\',p], 'output');
% save(['.\data\',p], 'output');


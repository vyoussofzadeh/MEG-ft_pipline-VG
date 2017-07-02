clc; clear; close all; 

%%
restoredefaultpath
addpath(genpath('E:\My Matlab\My codes\My GitHub\fieldtrip'));
ft_defaults

%%
load .\data\sub1
%%
timeAll = output.timelockanalysis.Verbs_Data;
timeBaseline = output.timelockanalysis.Verbs_post;
timePost = output.timelockanalysis.Verbs_Baseline;

figure
cfg = [];
cfg.layout = 'CTF151.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, timeAll)

%%
% source_diff_int  = output.source.source_diff_int;
% vol              = output.mri.vol;
% ft_021_sa_vis
ft_5_sourceanalysis
disp('Source analysis was visualised');

%% laterality check
% ft_12_laterality

%% Source + Conn/netwrok analysis (stat-based)
ft_6_stat
% disp('Stat was computed & visualised');

%%
% pause,
% source_conn_diff =  output.conn.source_conn_diff;
% network          =  output.conn.network;
% network_int      =  output.conn.network_int;
% ft_022_ca_vis
ft_7_connanalysis
disp('Conn analysis was visualised');

%% parcellation
% pause,
output_par_conn  =  output.par.output_par_conn;
% output_par_conn.eigenvector_cent = output_par_conn.eigenvector_cent - min(output_par_conn.eigenvector_cent);
atlas = ft_read_atlas('ROI_MNI_V4.nii');
% parcellation = output.par;
ft_023_pa_vis
disp('Parcellation analysis was visualised');

%%
% ft_11_freq



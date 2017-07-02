% clc; clear;
%
% ft_defaults
%
% sub = input('subject number (e.g = 1)?');
% p = ['sub',num2str(sub)];

% load('.\data\ft_data_preprocess.mat');
% load(['.\data\',p]);

Verbs_Data = output.preprocess.Verbs_Data;
Verbs_post = output.preprocess.Verbs_post;
Verbs_Baseline = output.preprocess.Verbs_Baseline;

%%
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
% cfg.vartrllength     = 0;
cfg.keeptrials       = 'yes';
timeAll        = ft_timelockanalysis(cfg, Verbs_Data);

timeBaseline = ft_timelockanalysis(cfg, Verbs_Baseline);
timePost = ft_timelockanalysis(cfg, Verbs_post);

%%
% cfg           = [];
% cfg.output    = 'fourier';
% cfg.keeptrials= 'yes';
% cfg.method    = 'mtmfft';
% cfg.taper     = 'dpss';
% cfg.foilim    = [1 50];
% cfg.tapsmofrq = 3;
% freqAll        = ft_freqanalysis(cfg, Verbs_Data);
% freqPost       = ft_freqanalysis(cfg, Verbs_post);
% freqBaseline   = ft_freqanalysis(cfg, Verbs_Baseline);
% 
% fg = [];
% cfg.output     = 'pow';
% cfg.channel    = 'MEG';
% cfg.method     = 'mtmconvol';
% cfg.foi        = 1:1:30;
% cfg.t_ftimwin  = 5./cfg.foi;
% cfg.tapsmofrq  = 0.4 *cfg.foi;
% cfg.toi        = -0.3:0.05:1;
% TFRmult = ft_freqanalysis(cfg, Verbs_Data);
% 
% 
% cfg = [];
% cfg.baseline     = [-0.3 -0.1];  
% cfg.baselinetype = 'absolute';	        
% cfg.zlim         = [-3e-27 3e-27];	
% cfg.showlabels   = 'yes';	        
% cfg.layout       = 'CTF151.lay';
% figure
% ft_multiplotTFR(cfg, TFRmult)


%% saving data
output.timelockanalysis.Verbs_Data     = timeAll;
output.timelockanalysis.Verbs_post     = timeBaseline;
output.timelockanalysis.Verbs_Baseline = timePost;
output.sub = sub;

%%
% output.freqanalysis.Verbs_Data     = freqAll;
% output.freqanalysis.Verbs_post     = freqBaseline;
% output.freqanalysis.Verbs_Baseline = freqPost;
%%

% save('.\data\ft_timelock.mat', 'timeAll','timeBaseline','timePost');
% save(['.\data\',p], 'output');


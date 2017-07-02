% clc; clear;
% 
% ft_defaults
% 
% % dataset  = 'Y:\controls\CTL01\MEG\CTL01_verbs.ds';
% dataset = 'Y:\controls\CTL02\MEG\CTL02_R-kadis-verbsER_20150403_01.ds';
% % dataset  = 'Y:\controls\CTL05\MEG\CTL0520150522_R-kadis-verbsER_20150522_01.ds';
% % dataset  = 'Y:\controls\CTL07\MEG\CTL0720150507_R-kadis-verbsER_20150507_01.ds';
% % dataset  = 'Y:\controls\CTL08\MEG\CTL08_R-kadis-verbsER_20150518_01.ds';
% dataset  = 'Y:\controls\CTL09\MEG\CTL09_R-kadis-verbsER_20150514_01.ds';
% dataset = 'Y:\controls\CTL15\MEG\CTL15_20160219\CTL15_R-kadis-verbsER_20160219_01.ds';
% dataset = 'Y:\controls\CTL19\MEG\CTL19_R-kadis-verbsER_20150910_02.ds';
% dataset = 'Y:\controls\CTL14\MEG\CTL1420150618_R-kadis-verbsER_20150618_01.ds';
% dataset = 'Y:\controls\CTL23\MEG\CTL23_R-kadis-verbsER_20150904_01.ds';
% disp(dataset(13:17))
% 
% sub = input('subject number (e.g = 1)?');

cfg = [];
cfg.dataset = dataset;
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'noun';
cfg.trialdef.prestim = 0.45;
cfg.trialdef.poststim = 1.25;
cfg = ft_definetrial(cfg);

cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfiltord = 3;
cfg.dftfreq = [60 120 180];
% cfg.hpfreq = 0.1;
% cfg.lpfreq = 100;
cfg.hpfreq = 10;
cfg.lpfreq = 30;
cfg.channel = {'MEG', '-MLC12'};
cfg.demean = 'yes';
cfg.baselinewindow = [-0.45 0.0];
Verbs_Data = ft_preprocessing(cfg);


% Verbs_Data = output.preprocess.Verbs_Data;

% correct for jump artificats
[cfg, artifact] = ft_artifact_jump(cfg);    % id jumps
Verbs_Data = ft_rejectartifact(cfg, Verbs_Data);        % remove scanner (MEG) jumps

% cfg = [];
% Verbs_Data2 = ft_rejectvisual(cfg, Verbs_Data);
% 
% cfg = [];
% cfg.artfctdef.zvalue.fltpadding
% cfg.artfctdef.zvalue.trlpadding
% [cfg, artifact_EOG] = ft_artifact_zvalue(cfg);
% 
% 
% cfg = [];
% cfg = ft_databrowser(cfg, Verbs_Data);

% cfg = [];
% cfg.method = 'pca';
% cfg.updatesens = 'no';
% cfg.channel = {'MEG', '-MLC12'};
% comp = ft_componentanalysis(cfg, Verbs_Data);
% 
% cfg = [];
% cfg.updatesens = 'no';
% cfg.component = comp.label;
% Verbs_Data = ft_rejectcomponent(cfg, comp);

cfg = [];
cfg.toilim = [-0.4 0.0];
Verbs_Baseline = ft_redefinetrial(cfg, Verbs_Data);
cfg.toilim = [0.6 1.0];
Verbs_post = ft_redefinetrial(cfg, Verbs_Data);

%%
% fg = [];
% cfg.output     = 'pow';
% cfg.channel    = 'MEG';
% % cfg.channel      = 'MLT33';
% cfg.method     = 'mtmconvol';
% cfg.foi        = 1:1:50;
% cfg.t_ftimwin  = 5./cfg.foi;
% cfg.tapsmofrq  = 0.4 *cfg.foi;
% cfg.toi        = -0.5:0.05:1.2;
% TFRmult = ft_freqanalysis(cfg, Verbs_Data);
% 
% cfg = [];
% cfg.baseline     = [-0.5 -0];  
% cfg.baselinetype = 'absolute';	        
% cfg.zlim         = [-3e-27 3e-27];	
% cfg.showlabels   = 'yes';	        
% cfg.layout       = 'CTF151.lay';
% % cfg.channel      = 'MLT33';
% figure
% ft_multiplotTFR(cfg, TFRmult)

%%
fg = [];
cfg.output     = 'pow';
cfg.channel    = 'MEG';
cfg.layout       = 'CTF151.lay';
cfg.channel      = 'MLT33';
cfg.method     = 'mtmconvol';
cfg.foi        = 1:1:50;
cfg.t_ftimwin  = 3./cfg.foi;
cfg.tapsmofrq  = 0.8 *cfg.foi;
cfg.toi        = -0.4:0.05:1;
TFRmult = ft_freqanalysis(cfg, Verbs_Data);

cfg = [];
cfg.baseline     = [-0.9 0];  
cfg.baselinetype = 'absolute';	        
% cfg.zlim         = [-3e-27 3e-27];	
cfg.showlabels   = 'yes';	        
cfg.layout       = 'CTF151.lay';
% cfg.channel      = 'MLT33';
figure
ft_multiplotTFR(cfg, TFRmult)
xlabel('Time (ms)');
ylabel('Freq (Hz)');

%%
% fg = [];
% cfg.output     = 'pow';
% cfg.channel    = 'MEG';
% % cfg.channel      = 'MLT33';
% cfg.method     = 'mtmconvol';
% cfg.foi        = 1:1:50;
% cfg.t_ftimwin  = 2./cfg.foi;
% cfg.tapsmofrq  = 0.6 *cfg.foi;
% cfg.toi        = 0.6:0.05:1;
% TFRmult = ft_freqanalysis(cfg, Verbs_post);
% 
% cfg = [];
% % cfg.baseline     = [-0.9 0];  
% cfg.baselinetype = 'absolute';	        
% % cfg.zlim         = [-3e-27 3e-27];	
% cfg.showlabels   = 'yes';	        
% cfg.layout       = 'CTF151.lay';
% % cfg.channel      = 'MLT33';
% figure
% ft_multiplotTFR(cfg, TFRmult)
% xlabel('Time (ms)');
% ylabel('Freq (Hz)');

%%
% cfg              = [];
% cfg.output       = 'pow';
% % cfg.channel      = 'MRC15';
% cfg.method       = 'mtmconvol';
% cfg.taper        = 'hanning';
% cfg.foi          = 1:1:30;
% cfg.t_ftimwin    = 7./cfg.foi;  % 7 cycles per time window
% cfg.toi          = -0.4:0.05:1.2;
% TFRhann = ft_freqanalysis(cfg, Verbs_post);

% cfg              = [];
% cfg.baseline     = [-0.5 -0.1]; 
% cfg.baselinetype = 'absolute'; 
% cfg.maskstyle    = 'saturation';	
% cfg.zlim         = [-3e-27 3e-27];	
% % cfg.channel      = 'MRC15';
% % cfg.interactive  = 'no';
% figure
% ft_singleplotTFR(cfg, TFRhann);

%% saving data
output.preprocess.Verbs_Data     = Verbs_Data;
output.preprocess.Verbs_post     = Verbs_post;
output.preprocess.Verbs_Baseline = Verbs_Baseline;
output.sub = sub;
output.dataset = dataset;

% p = ['sub',num2str(sub)];
% save(['.\data\',p], 'output'); 

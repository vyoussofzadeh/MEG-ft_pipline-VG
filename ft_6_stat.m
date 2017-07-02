%% tutorial link: http://www.fieldtriptoolbox.org/tutorial/salzburg

timeBaseline   = output.timelockanalysis.Verbs_post;
timePost       = output.timelockanalysis.Verbs_Baseline;
headmodel_time = output.headmodel.headmodel_time;
vol            = output.mri.vol;
timeAll        = output.timelockanalysis.Verbs_Data;
Verbs_Data     = output.preprocess.Verbs_Data;
Verbs_post     = output.preprocess.Verbs_post;
Verbs_Baseline = output.preprocess.Verbs_Baseline;

cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.vartrllength     = 0;
cfg.covariance       = 'yes';
cfg.keeptrials       = 'yes';
timeAll              = ft_timelockanalysis(cfg, Verbs_Data);
timeBaseline         = ft_timelockanalysis(cfg, Verbs_Baseline);
timePost             = ft_timelockanalysis(cfg, Verbs_post);

cfg = [];
cfg.method = 'lcmv';
cfg.lcmv.lambda = '0.1%';
cfg.grid = headmodel_time;
cfg.vol = vol;
cfg.lcmv.keepfilter = 'yes';
cfg.channel = timeAll.label;
sourceavg = ft_sourceanalysis(cfg, timeAll);

cfg = [];
cfg.method = 'lcmv';
cfg.lcmv.lambda = '0.1%';
cfg.grid = headmodel_time;
cfg.grid.filter = sourceavg.avg.filter;
cfg.rawtrial = 'yes';
cfg.vol = vol;
sourcepreS1 = ft_sourceanalysis(cfg, timeBaseline);
sourcepstS1 = ft_sourceanalysis(cfg, timePost);

%%
cfg = [];
cfg.parameter        = 'pow';
cfg.dim              = headmodel_time.dim;
cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.statistic        = 'depsamplesT';
% cfg.statistic         = 'indepsamplesT'; 
% cfg.correctm         = 'cluster';
cfg.correctm         = 'fdr';
cfg.clusteralpha     = 0.001;
% cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;

% cfg = [];
% cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability
% cfg.statistic         = 'indepsamplesT';        % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
% cfg.correctm          = 'cluster';
% cfg.clusteralpha      = 0.05;                   % alpha level of the sample-specific test statistic that will be used for thresholding
% cfg.clustertail       = 0;
% cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
% cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
% cfg.correcttail       = 'prob';                 % the two-sided test implies that we do non-parametric two tests
% cfg.alpha             = 0.05;                   % alpha level of the permutation test
% cfg.numrandomization  = 1000;                   % number of draws from the permutation distribution
% cfg.design            = TFRwavelet.trialinfo(:,1)'; % design matrix, note the transpose
% cfg.ivar              = 1;                      % the index of the independent variable in the design matrix
% cfg.channel           = {'MEG1243'};
% cfg.neighbours        = [];                     % there are no spatial neighbours, only in time and frequency
 
ntrials                       = numel(sourcepreS1.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = 1:ntrials;
design(2,ntrials+1:2*ntrials) = 1:ntrials;
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, timeAll);
 
cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat         = ft_sourcestatistics(cfg,sourcepstS1,sourcepreS1);
stat.pos     = headmodel_time.pos;% keep positions for plotting later
mri_aligned  = output.mri.mri_aligned;

%%
stat.inside = headmodel_time.inside;

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
% cfg.interpmethod = 'nearest';
cfg.interpmethod = 'sphere_avg';
statint  = ft_sourceinterpolate(cfg, stat, mri_aligned);
cfg.parameter    = 'mask';
maskint  = ft_sourceinterpolate(cfg, stat, mri_aligned);

cfg = [];
statint = ft_volumenormalise(cfg, statint);
cfg = [];
maskint = ft_volumenormalise(cfg, maskint);

statint.mask = maskint.mask;

atlas = ft_read_atlas('ROI_MNI_V4.nii');

statint.coordsys  = 'mni';
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
cfg.atlas         = atlas;
cfg.location      = 'max';
cfg.funcolorlim   = [-5 5];
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg,statint);

%%
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
view ([-70 20 50])
light ('Position',[-70 20 50])
% title('Network degrees')
colormap jet


cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
view ([70 20 50])
light ('Position',[70 20 50])
% title('Network degrees')
colormap jet
%% saving data
% output.stat.statint     = statint;
% output.stat.design      = design;
% save(['.\data\',p], 'output');
% save('output');

%% Conn analysis
cfg = [];
cfg.projectmom   = 'yes';
sourcepreS11 = ft_sourcedescriptives(cfg,sourcepreS1);
sourcepstS11 = ft_sourcedescriptives(cfg,sourcepstS1);
%%
sourcepreS12 = ft_source2sparse(sourcepreS11);
sourcepstS12 = ft_source2sparse(sourcepstS11);
%%
cfg = [];
% cfg.method       = 'coh';
% cfg.method         = 'wpli';
% cfg.method       = 'wpli_debiased';
% cfg.complex      = 'absimag';
cfg.method      = 'plv';
source_conn_pre = ft_connectivityanalysis(cfg, sourcepreS12);
source_conn_post = ft_connectivityanalysis(cfg, sourcepstS12);
%% subtraction ...
source_conn_diff = source_conn_post;
% source_conn_diff.cohspctrm = (source_conn_post.cohspctrm - source_conn_pre.cohspctrm) ...
%     ./ (source_conn_post.cohspctrm + source_conn_pre.cohspctrm);

% source_conn_diff.wpli_debiasedspctrm = (source_conn_post.wpli_debiasedspctrm - source_conn_pre.wpli_debiasedspctrm) ...
%     ./ (source_conn_post.wpli_debiasedspctrm + source_conn_pre.wpli_debiasedspctrm);

% source_conn_diff.wplispctrm = (source_conn_post.wplispctrm - source_conn_pre.wplispctrm) ...
%     ./ (source_conn_post.wplispctrm + source_conn_pre.wplispctrm);

source_conn_diff.plvspctrm = (source_conn_post.plvspctrm - source_conn_pre.plvspctrm) ...
    ./ (source_conn_post.plvspctrm + source_conn_pre.plvspctrm);
%%
source_conn_full = ft_source2full(source_conn_diff);
source_conn_full.dimord    = 'pos_pos';

% figure;imagesc(source_conn_full.cohspctrm); 
% figure;imagesc(source_conn_full.wpli_debiasedspctrm); 
% figure;imagesc(source_conn_full.wplispctrm);
figure;imagesc(source_conn_full.plvspctrm); 
%%
in2 = 2;
% computing graph metric
if in2 == 1;
    gtm    = 'degrees';
elseif in2 == 2
    gtm    = 'eigenvector_cent';
elseif in2 == 3
    gtm    = 'betweenness';
elseif in2 == 4
    gtm    = 'efficiency_bin';
end

cfg = [];
cfg.method    = gtm;
% cfg.parameter = 'cohspctrm';
cfg.parameter = 'plvspctrm';
% cfg.parameter = 'wpli_debiasedspctrm';
% cfg.parameter = 'wplispctrm';
cfg.threshold = .5;
network = ft_networkanalysis(cfg,source_conn_full);
% if in2 ==2
%     network.eigenvector_cent(network.eigenvector_cent<0.0001) = 0;
% end
network.pos     = source_conn_full.pos;
network.dim     = source_conn_full.dim;
network.inside  = source_conn_full.inside;

if in2 ==1
        figure;bar(network.degrees);
elseif in2 == 2
    figure;bar(network.eigenvector_cent(network.inside));
elseif in2 == 3
    figure;bar(network.betweenness);
end
%% interpolate on anatomical mri and plot the result
cfg              = [];
cfg.parameter    = gtm;
% cfg.interpmethod = 'nearest';
cfg.interpmethod = 'sphere_avg';
network_int  = ft_sourceinterpolate(cfg, network, mri_aligned);

cfg = [];
network_int = ft_volumenormalise(cfg, network_int);
%%
cfg               = [];
cfg.method        = 'ortho';
cfg.funcolormap   = 'jet';
cfg.location      = 'max';
cfg.funparameter  = gtm;
cfg.crosshair = 'no';
ft_sourceplot(cfg, network_int);
set(gca,'color','none');
%%
cfg               = [];
cfg.funparameter  = gtm;
cfg.method = 'surface';
cfg.funcolormap   = 'jet';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest';
cfg.camlight       = 'no';
cfg.crosshair = 'no';
ft_sourceplot(cfg, network_int);
view ([-70 20 50]);
% view ([180,0])
light ('Position',[-70 20 50]);
title(['Network', gtm]);
set(gca,'color','none');

cfg               = [];
cfg.funparameter  = gtm;
cfg.method = 'surface';
cfg.funcolormap   = 'jet';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest';
cfg.camlight       = 'no';
cfg.projthresh     = 0.3;
cfg.crosshair = 'no';
ft_sourceplot(cfg, network_int);
% view ([180,0])
view ([-70 20 50])
% view ([70 20 50])
light ('Position',[-70 20 50])
title(['Network', gtm])
set(gca,'color','none')

cfg               = [];
cfg.funparameter  = gtm;
cfg.method = 'surface';
cfg.funcolormap   = 'jet';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest';
cfg.camlight       = 'no';
cfg.projthresh     = 0.3;
cfg.crosshair = 'no';
ft_sourceplot(cfg, network_int);
% view ([180,0])
view ([70 20 50])
light ('Position',[70 20 50]);
title(['Network', gtm]);
set(gca,'color','none');






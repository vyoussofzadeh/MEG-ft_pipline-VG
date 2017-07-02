atlas = ft_read_atlas('ROI_MNI_V4.nii');
figure,
imagesc(atlas.tissue(:,:,45))

cfg = [];
cfg.parameter    = gtm;
cfg.interpmethod = 'sphere_avg';
network_int_par  = ft_sourceinterpolate(cfg, network_int, atlas);
% network_int_par2  = ft_sourceinterpolate(cfg, network, atlas);

cfg = [];
cfg.method      = 'mean';
output_par_conn = ft_sourceparcellate(cfg, network_int_par, atlas);

% cfg = [];
% cfg.funparameter = gtm;
% % cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% light ('Position',[-70 20 50])
% ft_sourceplot(cfg, output_par_conn)

% output_par_conn.eigenvector_cent(isnan(output_par_conn.eigenvector_cent))=0;
output_par_conn.eigenvector_cent = zscore(output_par_conn.eigenvector_cent);
output_par_conn.eigenvector_centdimord   = 'chan';

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = gtm;
ft_sourceplot(cfg, output_par_conn);


cfg = [];
cfg.funparameter = gtm;
cfg.funcolormap   = 'jet';
cfg.method = 'surface';
% cfg.projthresh     = 0.6;
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
ft_sourceplot(cfg, output_par_conn)

[L, idx] = sort(output_par_conn.eigenvector_cent,'descend')

output_par_conn.label(idx(1:10))


output.par.output_par_conn    = output_par_conn;
save(['.\data\',p], 'output_p');

%%
% clc; clear; close all; warning('off');
% 
% ft_defaults
% 
% sub = input('subject number (e.g = 1)?');
% p = ['sub',num2str(sub)];
% 
% subcat = input('Kids (=1)or Teens (= 2)?');
% if subcat == 1
%     load(['.\data\Kids\',p]);
% else
%     load(['.\data\Teens\',p]);
% end
% 
% 
% ft_7_connanalysis
% 
% % network_int = output.conn.network_int;
% % network = output.conn.network;
% 
% 
% atlas = ft_read_atlas('ROI_MNI_V4.nii');
% figure,
% imagesc(atlas.tissue(:,:,45))
% 
% % load('atlas_MMP1.0_4k.mat');
% % atlas.pos = source_diff_int_par.pos; % otherwise the parcellation won't work
% 
% % 
% % atlas.anatomy = network_int.anatomy;
% % cfg = [];
% % atlas = ft_volumenormalise(cfg, atlas);
% 
% % afni = ft_read_atlas('E:\My Matlab\My codes\My GitHub\fieldtrip\template\atlas\afni\TTatlas+tlrc');
% % atlas = ft_read_atlas('E:\My Matlab\My codes\My GitHub\fieldtrip\template\atlas\spm_anatomy\AllAreas_v17_MPM');
% 
% 
% 
% %% Source interpolation with atlas
% 
% sourcePost = output.source.sourcePost;
% sourceBaseline = output.source.sourceBaseline;
% 
% cfg = [];
% cfg.parameter = 'pow';
% cfg.operation = '(x1-x2)/(x1+x2)';
% S1bl = ft_math(cfg,sourcePost,sourceBaseline);
% 
% cfg = [];
% cfg.parameter = 'pow';
% cfg.interpmethod = 'nearest'; 
% % cfg.parameter = 'tissue'; 
% source_diff_int_par = ft_sourceinterpolate(cfg, S1bl, atlas);
% 
% % cfg = [];
% % % cfg.parameter = 'pow';
% % cfg.interpmethod = 'nearest'; 
% % cfg.parameter = 'tissue'; 
% % s2 = ft_sourceinterpolate(cfg, atlas, S1bl);
% 
% cfg = [];
% cfg.method = 'mean';
% output_par_source = ft_sourceparcellate(cfg, source_diff_int_par, atlas);
% 
% 
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% ft_sourceplot(cfg, output_par_source)
% % view ([-103 3])
% light ('Position',[-70 20 50])
% title('Source activations, power');
% colormap jet
% view ([-70 20 50])
% light ('Position',[-70 20 50])
% set(gca,'color','none')
% 
% %% Connectivity interpolation with atlas
% cfg = [];
% cfg.parameter = gtm;
% cfg.interpmethod = 'sphere_avg';
% network_int_par = ft_sourceinterpolate(cfg, network, atlas);
% 
% % network_int_par.anatomy = atlas.tissue;
% % cfg = [];
% % network_int_par = ft_volumenormalise(cfg, network_int_par);
% 
% cfg = [];
% cfg.method = 'mean';
% output_par_conn = ft_sourceparcellate(cfg, network_int_par, atlas);
% % output_par_conn.degreesdimord = 'chan';
% output_par_conn.eigenvector_centdimord   = 'chan';
% 
% cfg = [];
% cfg.funparameter = gtm;
% cfg.method = 'ortho';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% light ('Position',[-70 20 50])
% ft_sourceplot(cfg, output_par_conn)
% 
% cfg = [];
% cfg.funparameter = gtm;
% cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% light ('Position',[-70 20 50])
% ft_sourceplot(cfg, output_par_conn)
% title('Conn activations, degrees');
% colormap jet
% view ([-90 0])
% % view ([-70 20 50])
% light ('Position',[-70 20 50])
% set(gca,'color','none')
% 
% 
% % x = find(ismember(output_par_conn.label,'Heschl_L'));
% % indxHGL = find(stat_atlas.tissue==x);
% % 
% % x=find(ismember(atlas.tissuelabel,'Heschl_R'));
% % indxHGR = find(stat_atlas.tissue==x); 
% %  
% % x=find(ismember(atlas.tissuelabel,'Cingulum_Mid_L'));
% % indxCML = find(stat_atlas.tissue==x); 
% % 
% 
% %%
% con_gtm = output_par_conn.eigenvector_cent;
% con_gtm(isnan(con_gtm)) = 0 ;
% 
% % L = length(output_par_conn.degrees);
% L = length(con_gtm);
% id = table([1:L]','VariableNames',{'id'});
% ROI = [id,cell2table(output_par_conn.label')];
% ROI.Properties.VariableNames{'Var1'} = 'ROI';
% ROI
% figure,
% % barh(output_par_conn.degrees);
% barh(con_gtm);
% set(gca,'Ytick', 1:L,'YtickLabel',1:L);
% box off
% set(gca,'color','none');
% ylim([0,L+1])
% ylabel('ROI');
% xlabel('GTM');
% set(gcf, 'Position', [500   100   500   1200]);
% % legend({'Age','EVT','PPVT'})
% 
% 
% %% Frontal lobe
% idx = [3:16,23:26,79:90]';
% L = length(idx);
% figure,
% % barh(output_par_conn.degrees(idx));
% barh(con_gtm(idx));
% set(gca,'Ytick', 1:L,'YtickLabel',num2str(idx));
% box off
% set(gca,'color','none');
% ylim([0,L])
% ylabel('ROI');
% xlabel('GTM');
% set(gcf, 'Position', [500   100   500   1200]);
% ROI(idx,:)
% title(['BothHs, sum for sub ', num2str(sub), ' is ', num2str(sum(con_gtm(idx)))]);
% 
% 
% %% Left hemisphre lobe
% idx = [3:2:15,23:2:25,79:2:89]';
% L = length(idx);
% figure,
% % barh(output_par_conn.degrees(idx));
% barh(con_gtm(idx));
% set(gca,'Ytick', 1:L,'YtickLabel',num2str(idx));
% box off
% set(gca,'color','none');
% ylim([0,L])
% ylabel('ROI');
% xlabel('GTM');
% set(gcf, 'Position', [500   100   500   1200]);
% ROI(idx,:)
% title(['LH, sum for sub ', num2str(sub), ' is ', num2str(sum(con_gtm(idx)))]);
% 
% 
% %% Right hemisphre lobe
% idx = [4:2:16,24:2:26,80:2:90]';
% L = length(idx);
% figure,
% % barh(output_par_conn.degrees(idx));
% barh(con_gtm(idx));
% set(gca,'Ytick', 1:L,'YtickLabel',num2str(idx));
% box off
% set(gca,'color','none');
% ylim([0,L])
% ylabel('ROI');
% xlabel('GTM');
% set(gcf, 'Position', [500   100   500   1200]);
% ROI(idx,:)
% title(['RH, sum for sub ', num2str(sub), ' is ', num2str(sum(con_gtm(idx)))]);
% 
% %% Right - Left hemisphre lobe
% % idx_R = [4:2:16,24:2:26,80:2:90]';
% % L = length(idx_R);
% % figure,
% % % barh(output_par_conn.degrees(idx));
% % barh(con_gtm(idx));
% % set(gca,'Ytick', 1:L,'YtickLabel',[num2str(idx_R),'-', num2str(idx_L)]);
% % box off
% % set(gca,'color','none');
% % ylim([0,L])
% % ylabel('ROI');
% % xlabel('GTM');
% % set(gcf, 'Position', [500   100   500   1200]);
% % ROI(idx,:)
% % title(['RH, sum for sub ', num2str(sub), ' is ', num2str(sum(con_gtm(idx)))]);
% 
% 
% % Writing degrees of parcellated conn mat as an image 
% % cfg = [];
% % cfg.filetype  = 'nifti';
% % cfg.parameter = 'degrees';
% % cfg.filename  = ['.\data\par_',p];
% % % ['.\data\ft_output_par_conn'];
% % ft_volumewrite(cfg, output_par_conn)
% % 
% % % Writing degrees of conn mat as an image 
% % cfg = [];
% % cfg.filetype  = 'nifti';
% % cfg.parameter = 'degrees';
% % cfg.filename  = ['.\data\conn_',p];
% % % '.\data\ft_output_conn';
% % ft_volumewrite(cfg, network_int)
% % 
% % %% saving data
% % output.par    = atlas;
% % output.par.output_par_conn    = output_par_conn;
% % save(['.\data\',p], 'output');



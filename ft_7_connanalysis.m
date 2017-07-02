% clc; clear; close all;
%
% ft_defaults
%
% sub = input('subject number (e.g = 1)?');
% p = ['sub',num2str(sub)];
% load(['.\data\',p]);

sourcePost     = output.source.sourcePost;
sourceBaseline = output.source.sourceBaseline;
mri_aligned    = output.mri.mri_aligned;
vol            = output.mri.vol;
headmodel_time = output.headmodel.headmodel_time;

%% reduce the source reconstructed data to the dominant orientation
cfg = [];
cfg.projectmom   = 'yes';
source_proj_post = ft_sourcedescriptives(cfg,sourcePost);
source_proj_bsl  = ft_sourcedescriptives(cfg,sourceBaseline);

% %%
% spatialfilter = [];
% spatialfilter = cat(1,sourcePost.avg.filter{:});
% virtsens=[];
% for i=1:length(Verbs_post.trial)
%     virtsens.trial{i} = spatialfilter*Verbs_post.trial{i};
% end
% virtsens.time = Verbs_post.time;
% virtsens.fsample = Verbs_post.fsample;
% virtsens.label = cell(1627,1);
% A = sourcePost.pos(sourcePost.inside,:);
% for i=1:1627
%     virtsens.label{i} = num2str(A(i,:));
% end
% virtsens.trialinfo = ones(1,size(Verbs_post.trial,2));
% virtsens_post = virtsens;
% % 
% % % cfg=[];
% % % cfg.preproc.hpfilter='yes';
% % % cfg.preproc.hpfreq=1;
% % % cfg.preproc.lpfilter='yes';
% % % cfg.preproc.lpfreq=40;
% % % tlkvc=ft_timelockanalysis(cfg, virtsens);
% 
% spatialfilter = [];
% spatialfilter = cat(1,sourceBaseline.avg.filter{:});
% virtsens=[];
% for i=1:length(Verbs_post.trial)
%     virtsens.trial{i} = spatialfilter*Verbs_Baseline.trial{i};
% end
% virtsens.time = Verbs_post.time;
% virtsens.fsample = Verbs_post.fsample;
% virtsens.label = cell(1627,1);
% A = sourceBaseline.pos(sourceBaseline.inside,:);
% for i=1:1627
%     virtsens.label{i} = num2str(A(i,:));
% end
% virtsens.trialinfo = ones(1,size(Verbs_post.trial,2));
% virtsens_bsl = virtsens;
% %
% %%
% cfg           = [];
% cfg.output    = 'fourier';
% cfg.keeptrials= 'yes';
% cfg.method    = 'mtmfft';
% cfg.taper     = 'dpss';
% cfg.foilim    = [18 18];
% cfg.tapsmofrq = 3;
% verbs_post       = ft_freqanalysis(cfg, virtsens_post);
% verbs_bsl       = ft_freqanalysis(cfg, virtsens_bsl);
% 
% cfg           = [];
% % cfg.method    = 'wpli_debiased';
% cfg.complex = 'absimag';
% cfg.method  ='coh';
% verbs_wpli_post       = ft_connectivityanalysis(cfg, verbs_post);
% verbs_wpli_bsl       = ft_connectivityanalysis(cfg, verbs_bsl);
% 
% verbs_wpli = verbs_wpli_post;
% % verbs_wpli.wpli_debiasedspctrm = verbs_wpli_post.wpli_debiasedspctrm - verbs_wpli_bsl.wpli_debiasedspctrm ./ (verbs_wpli_post.wpli_debiasedspctrm + verbs_wpli_bsl.wpli_debiasedspctrm);
% % verbs_wpli.wpli_debiasedspctrm = verbs_wpli.wpli_debiasedspctrm./max(verbs_wpli.wpli_debiasedspctrm(:));
% 
% 
% verbs_wpli.cohspctrm = verbs_wpli_post.cohspctrm - verbs_wpli_bsl.cohspctrm ./ (verbs_wpli_post.cohspctrm + verbs_wpli_bsl.cohspctrm);
% verbs_wpli.cohspctrm = verbs_wpli.cohspctrm./max(verbs_wpli.cohspctrm(:));
% 
% figure;imagesc(verbs_wpli.cohspctrm); 
% colorbar

%% reduce memory demands and compute connectivity

% compute the sparse representation
source_sparse_post = ft_source2sparse(source_proj_post);
source_sparse_bsl  = ft_source2sparse(source_proj_bsl);

%%
data = source_sparse_post;
sizmom = size(data.avg.mom{data.inside(1)});
mom = zeros(size(data.pos,1), sizmom(2));
mom(data.inside, :) = cat(1, data.avg.mom{data.inside});
[nvox, nrpt]   = size(mom);
crsspctrm = (mom*mom')./nrpt;
figure;imagesc(crsspctrm);
colorbar

mmom = mean(mom,1);
figure,
plot(mmom)



% then compute connectivity
cfg = [];
% cfg.method  ='coh';
% cfg.complex = 'abs';
% cfg.complex = 'imag';
% cfg.complex = '-logabs';
% cfg.complex = 'absimag';
cfg.method       = 'plv';
% cfg.method = 'wpli_debiased';
% cfg.method = 'wpli';
% cfg.method = 'psi';
% cfg.bandwidth = 10:20;
source_conn_post = ft_connectivityanalysis(cfg, source_sparse_post);
source_conn_bsl  = ft_connectivityanalysis(cfg, source_sparse_bsl);


source_conn_diff = source_conn_post;
% source_conn_diff.cohspctrm = source_conn_post.cohspctrm - source_conn_bsl.cohspctrm ./ source_conn_post.cohspctrm + source_conn_bsl.cohspctrm;
source_conn_diff.plvspctrm = source_conn_post.plvspctrm - source_conn_bsl.plvspctrm ./ (source_conn_post.plvspctrm + source_conn_bsl.plvspctrm);
% source_conn_diff.wpli_debiasedspctrm = source_conn_post.wpli_debiasedspctrm - source_conn_bsl.wpli_debiasedspctrm ./ (source_conn_post.wpli_debiasedspctrm + source_conn_bsl.wpli_debiasedspctrm);
% source_conn_diff.wplispctrm = source_conn_post.wplispctrm - source_conn_bsl.wplispctrm ./ (source_conn_post.wplispctrm + source_conn_bsl.wplispctrm);


%%
% figure;imagesc(source_conn_diff.plvspctrm);
%%
aedge =  source_conn_diff.plvspctrm;
% aedge =  source_conn_diff.cohspctrm;
% aedge =  source_conn_diff.wpli_debiasedspctrm;
% aedge =  source_conn_diff.wplispctrm;
% aedge = abs(aedge);
% aedge = -aedge;
% aedge =  source_conn_diff.cohspctrm;
tedge = (aedge.* double(aedge > 0.9.*max(max(aedge))));
for k = 1:length(tedge), lab{k} = num2str(k); end

% figure,
% plot_conn(tedge,lab, 'thresholded plv');
% set(gcf, 'Position', [900   400   800  800]);
%
vol1 = ft_convert_units(vol, 'cm');
Verticies = vol.bnd.pos;
ROI = source_conn_diff.pos;

% figure(1),
% plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'b*');
% hold on
% plot3(ROI(:,1),ROI(:,2),ROI(:,3),'r*');
% box off
% set(gca,'color','none');
% axis off
% % set(gcf, 'Position', [900   400   800  800]);
% axis image
% rotate3d on
% hold on
% view([-90,90])
%
% for i = 1:length(tedge)
%     for j = 1:length(tedge)
%         if tedge(i,j)> 0.5
%             p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
%             p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
%             pts = [p1; p2];
%             line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
%         end
%     end
% end

%% network analysis
source_conn_full = ft_source2full(source_conn_diff);
source_conn_full.dimord    = 'pos_pos';

figure;imagesc(source_conn_diff.plvspctrm);

% figure;imagesc(source_conn_diff.cohspctrm);
% figure;imagesc(source_conn_diff.wpli_debiasedspctrm); 
% figure;imagesc(source_conn_diff.wplispctrm); 
colorbar
xlabel('voxels')
ylabel('voxels')
title('PLV')


% disp('=================');
% disp('degrees          = 1')
% disp('eigenvector_cent = 2')
% disp('betweenness      = 3');
% in2  = input('GT measure? ');
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
        figure;bar(network.degrees(network.inside));
elseif in2 == 2
    figure;bar(network.eigenvector_cent(network.inside));
    mm = network.eigenvector_cent(network.inside);
    mm = mm./max(mm);
    figure,
    counts = hist(mm,length(mm));
    binWidth = 0.1;
    binCtrs = 0:binWidth:1;
    paramEsts = wblfit(mm);
    n = length(mm);
    prob = counts / (n * binWidth);
    %     bar(prob,'hist');
    histogram(mm,100);
    h_gca = gca;
    h = h_gca.Children;
    h.FaceColor = [.9 .9 .9];
    xlabel('Eigenvector centrality');
    ylabel('Frequency');
    %     ylim([0 0.1]);
    xlim([0.02 1]);
    ylim([0 50]);
    xgrid = linspace(0,1,100);
    pdfEst = wblpdf(xgrid,paramEsts(1),paramEsts(2));
    %     line(xgrid,pdfEst,'color','k')
    
    paramEsts = wblfit(mm);
    xgrid = linspace(0,1800,200);
    pdfEst = wblpdf(xgrid,paramEsts(1),paramEsts(2));
    hold on
    line(xgrid,pdfEst,'color','k')
    
    
    normplot(mm)
    histfit(mm,200,'Normal')
    hold on
    qqplot(mm)
elseif in2 == 3
    figure;bar(network.betweenness(network.inside));
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
% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = gtm;
% cfg.filename  = ['.\data\Net_',p];
% ft_volumewrite(cfg, network_int);

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

%%
network_int1 = network_int;
thre = 0.5;
network_int1.eigenvector_cent(network_int1.eigenvector_cent < thre*max(network_int1.eigenvector_cent(:)))=0;
network_int1.eigenvector_cent(network_int1.eigenvector_cent > 0) = 1;
cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = gtm;
% cfg.filename  = '.\output';
cfg.filename  = 'output1.nii';
ft_volumewrite(cfg, network_int1)
%%
network_int2 = network_int;
thre = 0.7;
network_int2.eigenvector_cent(network_int2.eigenvector_cent < thre*max(network_int2.eigenvector_cent(:)))=0;
network_int2.eigenvector_cent(network_int2.eigenvector_cent > 0) = 1;
cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = gtm;
% cfg.filename  = '.\output';
cfg.filename  = 'output.nii';
ft_volumewrite(cfg, network_int2)

%%
% network_int3 = network_int;
% thre = 0.9;
% network_int3.eigenvector_cent(network_int3.eigenvector_cent < thre*max(network_int3.eigenvector_cent(:)))=0;
% network_int3.eigenvector_cent(network_int3.eigenvector_cent > 0) = 1;
% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = gtm;
% % cfg.filename  = '.\output';
% cfg.filename  = 'output3.nii';
% ft_volumewrite(cfg, network_int3)

% cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_left_caret.mat';
% cfg.projmethod     = 'nearest';
% cfg.camlight       = 'no';
% cfg.projthresh     = 0.5;
% cfg.crosshair = 'no';
% ft_sourceplot(cfg, network_int);
% view ([-70 20 50])
% light ('Position',[-70 20 50])
% title(['Network', gtm])
% set(gca,'color','none')
%%
figure,
plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'color',[0.7,0.7,0.7]);
hold on
h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.g');
% set(h(1),'MarkerEdgeColor',[1 0.48 0.30],'MarkerFaceColor','g')
set(h(1),'MarkerEdgeColor','g','MarkerFaceColor','g')
% set(h(2),'MarkerEdgeColor','none','MarkerFaceColor','g')
box off
set(gca,'color','none');
axis off
axis image
rotate3d on
hold on
view([-90,90])
% view ([180 90])

% figure,
% ft_plot_vol(vol, 'facecolor', 'none', 'edgecolor', [0.7,0.7,0.7]); alpha 0.5; camlight;
% hold on
% h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.k');
% set(h(1),'MarkerEdgeColor',[1 0.48 0.30],'MarkerFaceColor','k')
% view ([-196 56])

if in2 == 1
    [L, Lidx] = sort(network.degrees(network.inside),'descend');
elseif in2 == 2
    [L, Lidx] = sort(network.eigenvector_cent(network.inside),'descend');
elseif in2 == 3
    [L, Lidx] = sort(network.betweenness(network.inside),'descend');
end

% tedge = verbs_wpli.cohspctrm;
nROI = 50;
for i = 1:nROI
    plot3(ROI(Lidx(i),1),ROI(Lidx(i),2),ROI(Lidx(i),3),'g.', 'MarkerSize',32);
end
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0.5
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
        end
    end
end

% cfg               = [];
% cfg.funcolormap   = 'jet';
% cfg.location      = 'max';
% cfg.funparameter  = gtm;
% cfg.method        = 'surface';
% cfg.surfinflated = vol.bnd;
% ft_sourceplot(cfg, network_int);
% alpha 0.4
% colorbar off
% % view([-90,90]);
% % view ([-196 56])
% view ([-263,20])

% vol2 = ft_convert_units(vol, 'mm');
figure,
ft_plot_vol(vol, 'facecolor', 'cortex', 'edgecolor', 'none'); alpha 0.5; camlight;
view ([-196 56])
% view ([27 39])

hold on
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0.5
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
        end
    end
end
hold on
for i = 1:nROI
    plot3(ROI(Lidx(i),1),ROI(Lidx(i),2),ROI(Lidx(i),3),'g.', 'MarkerSize',32);
    id(i) = Lidx(i);
end
set(gca,'color','none')
% title(['sum gtm for sub ', num2str(sub), ' is ', num2str(sum(L(1:nROI)))]);
title('');

%%
% addpath('E:\My Matlab\Connectivity\Conn\conn');
% h = get(0, 'Children');
% if isempty(findobj(h,'tag','CONN functional connectivity toolbox'));
%     conn
% end
% addpath(genpath('E:\My Matlab\SPM\spm12_3\spm12'));
% addpath(genpath('E:\My Matlab\My codes\My GitHub\My tool\MEG_ft_pipline\Time domain'));
% %  ROI2 = cor2mni(10.*ROI, network_int.transform);
%  ROI2 = cor2mni(10.*ROI, network_int.initial );
% %  ROI2 = cor2mni(ROI, network_int.initial );
% conn_mesh_display('', '', '', ROI2(id,:), tedge(id,id), .1);
% %%
% filenameSURF = 'output1.nii';
% filenameVOL = [];
% FSfolder = 'E:\My Matlab\Connectivity\Conn\conn17a\conn\utils\surf';
% sphplots = [];
% connplots = [];
% facealpha = 1;
% position = [-1 0  0];
% conn_mesh_display(filenameSURF,[],FSfolder)
% %%
% filenameSURF = 'output2.nii';
% conn_mesh_display(filenameSURF,[],FSfolder)
% %%
% filenameSURF = 'output3.nii';
% conn_mesh_display(filenameSURF,[],FSfolder)
%% BrainNet viewer
% we increase the threshold here to highlight the dominant links
% edge = source_conn_diff.cohspctrm >.9;
% dlmwrite('.\data\edge.edge',edge,'\t');
%
% [L, Lidx] = sort(network.degrees,'descend');
% nROI = 100;
% edge2 = [];
% edge = source_conn_diff.cohspctrm;
% edge2 = edge(Lidx(1:nROI),:);
%
%
% node = zeros(length(network.degrees(network.inside)),6);
% node(:,1:3) = network.pos(network.inside,:);
% node(:,4)   = 4;
% node(:,5)   = network.degrees(network.inside);
% node(:,6)   = 0;
% node(:,1:3) = node(:,1:3)*10;
% dlmwrite('.\data\node.node',node,' ');
%
% BrainNet_MapCfg('.\data\mesh.nv','.\data\node.node','.\data\edge.edge');
% view([0 -90 0])


%% Clustering
% ROI_corr = [];
% nROI = length(Lidx);
% nROI = 50;
% cl  = [ROI(Lidx(1:nROI),1),ROI(Lidx(1:nROI),2),ROI(Lidx(1:nROI),3)];
% clv = L(1:nROI);
% % clusters = 30;
% clusters = 3;
% obj = fitgmdist(cl,clusters);
% idx = cluster(obj,cl);
% 
% % threshold = 0.8.*max(L);
% % idx = cluster(cl,'cutoff', threshold);
% 
% cfg               = [];
% cfg.funcolormap   = 'jet';
% cfg.location      = 'max';
% cfg.funparameter  = gtm;
% cfg.method        = 'surface';
% cfg.surfinflated  = vol1.bnd;
% ft_sourceplot(cfg, network_int);
% alpha 0.4
% colorbar off
% view ([-263,20])
% hold on
% clr = [];
% scclv = [];
% for i=1:clusters
%     clr(i,:) = rand(3,1);
%     plot3(cl(idx == i,1),cl(idx == i,2),cl(idx == i,3),'*','color', clr(i,:));
%     mcoor = mean(cl(idx == i,:));
%     plot3(mcoor(1),mcoor(2),mcoor(3),'r.', 'MarkerSize',32);
%     % mean gtm
%     cclv = clv(idx == i);
%     %     text(mcoor(1),mcoor(2),mcoor(3), ['c',num2str(i),' (',num2str(sum(cclv)),')']);
%     %     ROI_coor(i,:) = round(mean([cl(idx == i,1),cl(idx == i,2),cl(idx == i,3)],1));
%     scclv(i) = sum(cclv);
% end
% box off
% set(gca,'color','none');
% axis off
% % set(gcf, 'Position', [900   400   800  800]);
% axis image
% rotate3d on
% view ([-90,90])
% 
% [m, midx] = sort(scclv,'descend');
% 
% % in = input('cluster of interest?');
% in  = midx(1);
% cfg               = [];
% cfg.funcolormap   = 'jet';
% cfg.location      = 'max';
% cfg.funparameter  = gtm;
% cfg.method        = 'surface';
% cfg.surfinflated  = vol.bnd;
% ft_sourceplot(cfg, network_int);
% alpha 0.4
% colorbar off
% % view ([-196 56])
% hold on
% for i=in
%     clr = rand(3,1);
%     plot3(cl(idx == i,1),cl(idx == i,2),cl(idx == i,3),'*','color', clr);
%     mcoor = mean(cl(idx == i,:));
%     plot3(mcoor(1),mcoor(2),mcoor(3),'r.', 'MarkerSize',62);
%     % mean gtm
%     cclv = clv(idx == i);
%     %     text(mcoor(1),mcoor(2)+1,mcoor(3), ['c',num2str(i),' (',num2str(sum(cclv)),')']);
%     %     ROI_coor(i,:) = round(mean([cl(idx == i,1),cl(idx == i,2),cl(idx == i,3)],1));
% end
% box off
% set(gca,'color','none');
% axis off
% % set(gcf, 'Position', [900   400   800  800]);
% axis image
% rotate3d on
% % view ([-263,20])
% view ([270,90])

%% saving data
% output.conn.network_int       = network_int;
% output.conn.network           = network;
% output.conn.source_conn_diff  = source_conn_diff;
% save(['.\data\',p], 'output');

%%
% load('atlas_MMP1.0_4k.mat');
% % atlas.pos = source_conn_full.pos; % otherwise the parcellation won't work
%
% cfg = [];
% cfg.parameter = 'plvspctrm';
% source_diff_int_par = ft_sourceinterpolate(cfg, source_conn_full, atlas);
%
%
% %
% cfg = [];
% cfg.parcellation = 'parcellation';
% cfg.parameter    = 'plvspctrm';
% parc_conn = ft_sourceparcellate(cfg, source_diff_int_par, atlas);
% parc_conn.dimord  = 'pos_pos';
% parc_conn.pos = source_conn_full.pos; % otherwise the parcellation won't work
%
% parc_conn = ft_source2sparse(parc_conn);
%
% A = parc_conn.plvspctrm;
% A(isnan(A(:,1)),:)=[];
%
% A(isnan(A)) = [];
% figure;imagesc(parc_conn.plvspctrm);
%
% %%
%
% cfg = [];
% cfg.method    = gtm;
% % cfg.parameter = 'cohspctrm';
% cfg.parameter = 'plvspctrm';
% cfg.threshold = .5;
% network = ft_networkanalysis(cfg,parc_conn);
%
%
%
% cfg           = [];
% cfg.method    = 'degrees';
% cfg.parameter = 'plvspctrm';
% cfg.threshold = .1;
% % network_full = ft_networkanalysis(cfg,source_conn);
% network_parc = ft_networkanalysis(cfg,source_conn_diff);
%
% %
% % figure;imagesc(parc_conn.plvspctrm);
%
%

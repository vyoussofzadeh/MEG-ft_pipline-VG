cfg = [];
cfg.funparameter = 'eigenvector_cent';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest'; 
ft_sourceplot(cfg, output_par_conn);
view ([-70 20 50])
light ('Position',[-70 20 50])
title('Network degrees')
colormap jet


cfg = [];
cfg.funparameter = 'eigenvector_cent';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest'; 
ft_sourceplot(cfg, output_par_conn);
view ([70 20 50])
light ('Position',[-70 20 50])
title('Network degrees')
colormap jet

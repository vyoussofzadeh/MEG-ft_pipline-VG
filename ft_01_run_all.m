clc; clear; close all;

ft_defaults

clc; clear; close all; warning('off');

ft_defaults

dataset = 'Y:\controls\CTL01\MEG\CTL01_verbs.ds';
% dataset  = 'Y:\controls\CTL02\MEG\CTL02_R-kadis-verbsER_20150403_01.ds';
% dataset  = 'Y:\controls\CTL03\MEG\CTL03_R-kadis-verbsER_20150417_01.ds';
% dataset  = 'Y:\controls\CTL06\MEG\CTL0620150430_R-kadis-verbsER_20150430_01.ds';
% dataset = 'Y:\controls\CTL08\MEG\CTL08_R-kadis-verbsER_20150518_01.ds';
% dataset  = 'Y:\controls\CTL09\MEG\CTL09_R-kadis-verbsER_20150514_01.ds';
% dataset = 'Y:\controls\CTL10\MEG\CTL1020150520_R-kadis-verbsER_20150520_01.ds';
% dataset = 'Y:\controls\CTL13\MEG\CTL1320150710_R-kadis-verbsER_20150710_01.ds';
% dataset = 'Y:\controls\CTL14\MEG\CTL1420150618_R-kadis-verbsER_20150618_01.ds';
% dataset = 'Y:\controls\CTL15\MEG\CTL15_20160219\CTL15_R-kadis-verbsER_20160219_01.ds';
% dataset = 'Y:\controls\CTL16\MEG\CTL1620150702_R-kadis-verbsER_20150702_02.ds';
% dataset = 'Y:\controls\CTL17\MEG\CTL1707242015_R-kadis-verbsER_20150724_01.ds';
% dataset = 'Y:\controls\CTL19\MEG\CTL19_R-kadis-verbsER_20150910_02.ds';
% dataset = 'Y:\controls\CTL21\MEG\CTL2120150803_R-kadis-verbsER_20150803_01.ds';
% dataset = 'Y:\controls\CTL22\MEG\CTL22_R-kadis-verbsER_20150807_01.ds';
% dataset = 'Y:\controls\CTL23\MEG\CTL23_R-kadis-verbsER_20150904_01.ds';
% dataset = 'Y:\controls\CTL25\MEG\CTL25_R-kadis-verbsER_20151002_01.ds';
% dataset = 'Y:\controls\CTL28\MEG\CTL28_20151205\CTL28_R-kadis-verbsER_20151205_01.ds';
% dataset = 'Y:\controls\CTL29\MEG\CTL29_20151106\CTL29_R-kadis-verbsER_20151106_01.ds';
% dataset = 'Y:\controls\CTL30\MEG\20160422\CTL30_R-kadis-verbsER_20160422_01.ds';
% dataset = 'Y:\controls\CTL31\MEG\CTL31_20151111\CTL31_R-kadis-verbsER_20151111_01.ds';
% dataset = 'Y:\controls\CTL32\MEG\CTL32_R-kadis-verbsER_20151117_01.ds';
% dataset = 'Y:\controls\CTL33\MEG\CTL33_R-kadis-verbsER_20151210_01.ds';
% dataset = 'Y:\controls\CTL36\MEG\CTL36_R-kadis-verbsER_20160121_01.ds';
% dataset = 'Y:\controls\CTL37\MEG\MRGCTL37_R-kadis-verbsER_20151228_01.ds';
% dataset = 'Y:\controls\CTL38\MEG\CTL38_R-kadis-verbsER_20160225_01.ds';

ev  = ft_read_event(dataset);

mripath = 'Y:\controls\CTL01\nMR\T1';

% sub = input('subject number (e.g = 1)?');
disp(str2double(mripath(16:17)));

%%
disp('preprocessing ...')
ft_1_preprocess
disp('preprocessing was completed');

%%
disp('Timelock ...')
ft_2_timelock
disp('Timelock was completed');

%%
disp('MRI ...')
% output.mri = output.mri;
% output.mri.vol     = output.mri.vol;
% output.mri.mri_aligned = output.mri.mri_aligned;
ft_3_mri
disp('MRI was completed');

%%
% cd '.\Time domain\Noun condition\';
disp('Forward model ...')
ft_4_forwardmodel
disp('Forward model was completed');

%%
% pause
disp('Source analysis ...')
ft_5_sourceanalysis
disp('Source analysis was completed');

%%
% disp('Source stat ...')
ft_6_stat
% disp('Source stat was completed');

%%
disp('Conn analysis ...')
ft_7_connanalysis
disp('Conn analysis was completed');

%%
disp('Parcellation ...')
ft_8_parcellation
% disp('Parcellation was completed');






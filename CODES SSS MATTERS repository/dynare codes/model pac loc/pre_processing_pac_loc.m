clear;

%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
addpath('C:\dynare\5.5\matlab'); %change path accordingly
dynare prepare_model_pac_loc.mod;
save('my_pac_loc_model.mat');

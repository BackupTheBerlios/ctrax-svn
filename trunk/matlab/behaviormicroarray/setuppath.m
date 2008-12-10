% set up the paths
if isempty(which('get_readframe_fcn')),
  tmp = which('fixerrorsgui');
  tmp = strrep(tmp,'fixerrorsgui.m','');
  if ~isempty(tmp),
    tmp = [tmp,'../matlab'];
    if exist(tmp,'dir'),
      addpath(genpath(tmp));
    end
  end
  % now test for one of the functions we need
  if isempty(which('get_readframe_fcn')),
    % ask for help
    fprintf('Where is the "matlab" directory of your Ctrax code?\n');
    dirname = uigetdir('.','Select Ctrax/matlab');
    if isnumeric(dirname) && dirname == 0,
      fprintf('Error setting up path.\n');
      return;
    end
    addpath(genpath(dirname));
    if isempty(which('get_readframe_fcn'))
      fprintf('Error setting up path.\n');
      return;
    end
  end
end
% set up the paths
isdone = ~isempty(which('get_readframe_fcn'));
if isdone, return; end

dirnamestry = {'.','matlab','../matlab','../Ctrax/matlab'};

% try last saved location
rcfile = which('setuppath');
rcfile = strrep(rcfile,'setuppath.m','.setuppathrc.mat');
if exist(rcfile,'file')
  load(rcfile);
  dirnamestry{1} = dirname;
end

for i = 1:length(dirnamestry),
  dirname = dirnamestry{i};
  if exist(dirname,'dir') && exist([dirname,'/filehandling'],'dir') && ...
        exist([dirname,'/misc'],'dir') && ...
        exist([dirname,'/filehandling/get_readframe_fcn.m'],'file'),
    addpath(genpath(dirname));
  end
  isdone = ~isempty(which('get_readframe_fcn'));
  if isdone, 
    fprintf('Found Ctrax/matlab at %s\n',dirname);
    save(rcfile,'dirname');
    return; 
  end
end

dirname = dirnamestry{1};

% ask for help
fprintf('Where is the "matlab" directory of your Ctrax code?\n');
dirname = uigetdir(dirname,'Select Ctrax/matlab');
if isnumeric(dirname) && dirname == 0,
  error('Error setting up path.\n');
  return;
end
addpath(genpath(dirname));
if isempty(which('get_readframe_fcn'))
  error('Error setting up path.\n');
end

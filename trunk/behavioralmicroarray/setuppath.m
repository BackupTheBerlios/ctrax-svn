% set up the paths
isdone = ~isempty(which('ctrax_matlab_misc_check')) &&  ~isempty(which('ctrax_matlab_filehandling_check'));
if isdone, return; end

dirnamestry = {'.','../matlab','..','matlab','../Ctrax/matlab'};

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
    p = genpath(dirname);
    % remove svn
    tmp = textscan(p,'%s','delimiter',':');
    tmp = tmp{1};
    matches = strfind(tmp,'.svn');
    keep = cellfun(@isempty,matches);
    tmp(~keep) = [];
    p1 = sprintf('%s:',tmp{:});
    addpath(p1);
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
p = genpath(dirname);
% remove svn
tmp = textscan(p,'%s','delimiter',':');
tmp = tmp{1};
matches = strfind(tmp,'.svn');
keep = cellfun(@isempty,matches);
tmp(~keep) = [];
p1 = sprintf('%s:',tmp{:});
addpath(p1);
if isempty(which('get_readframe_fcn'))
  error('Error setting up path.\n');
end

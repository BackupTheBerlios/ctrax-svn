% [trx,matname,succeeded] = load_tracks(matname,[moviename])

function [trx,matname,succeeded] = load_tracks(matname,moviename,varargin)

[dosave,savename,annname] = myparse(varargin,'dosave',false,'savename','','annname','');

succeeded = false;
trx = [];

isinteractive = ~exist('matname','var');

if isinteractive,
  helpmsg = 'Choose mat file containing trajectories to load';
  [matname,matpath] = uigetfilehelp('*.mat','Choose mat file containing trajectories','','helpmsg',helpmsg);
  if ~ischar(matname),
    return;
  end
  matname = [matpath,matname];
end

tmp = load(matname);
fprintf('loaded %s\n',matname);
if isfield(tmp,'pairtrx'),
  tmp.trx = tmp.pairtrx;
end
if ~isfield(tmp,'trx'),
  fprintf('no trx variable\n');
  if isfield(tmp,'ntargets'),
    fprintf('Ctrax output file; converting to trx file\n');
    if ~exist('moviename','var'),
      moviename = '?';
    end
    %ds = datestr(now,30);
    fprintf('Calling cleanup_ctrax_data\n');
    [trx,matname] = cleanup_ctrax_data(matname,moviename,tmp,'','dosave',dosave,'savename',savename,'annname',annname);
  else
    msgbox('Could not load data from %s, exiting',matname);
    return;
  end
else
  fprintf('trx variable found\n');
  trx = tmp.trx;
  if exist('moviename','var') && ~isfield(trx,'moviename'),
    for i = 1:length(trx),
      trx(i).moviename = moviename;
    end
  end
end

% member functions can be weird
fprintf('Adding off\n');
for i = 1:length(trx),
  trx(i).off = -trx(i).firstframe + 1;
  trx(i).matname = matname;
end

succeeded = true;
fprintf('Done. returning from load_tracks\n');

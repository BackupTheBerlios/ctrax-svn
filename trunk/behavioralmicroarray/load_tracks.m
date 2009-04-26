% [trx,matname,succeeded] = load_tracks(matname,[moviename])

function [trx,matname,succeeded] = load_tracks(matname,moviename)

succeeded = false;
trx = [];

isinteractive = ~exist('matname','var');

if isinteractive,
  [matname,matpath] = uigetfile('*.mat','Choose mat file containing trajectories','');
  if ~ischar(matname),
    return;
  end
  matname = [matpath,matname];
end

tmp = load(matname);
if isfield(tmp,'pairtrx'),
  tmp.trx = tmp.pairtrx;
end
if ~isfield(tmp,'trx'),
  if isfield(tmp,'ntargets'),
    if ~exist('moviename','var'),
      moviename = '?';
    end
    ds = datestr(now,30);
    [trx,matname] = cleanup_ctrax_data(matname,moviename,tmp,ds);
  else
    msgbox('Could not load data from %s, exiting',matname);
    return;
  end
else
  trx = tmp.trx;
  if exist('moviename','var') && ~isfield(trx,'moviename'),
    for i = 1:length(trx),
      trx(i).moviename = moviename;
    end
  end
end

% member functions can be weird
for i = 1:length(trx),
  trx(i).f2i = @(f) f - trx(i).firstframe + 1;
  trx(i).matname = matname;
end

succeeded = true;
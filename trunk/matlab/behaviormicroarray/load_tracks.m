function [trx,matname] = load_tracks(matname,moviename)

tmp = load(matname);
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
end
% succeeded = save_tracks(trx,matname,...)
% optional arguments:
% 'doappend' = false
% 'varname' = trx
function succeeded = save_tracks(trx,matname,varargin)

[doappend,varname] = myparse(varargin,'doappend',false,'varname','trx');

succeeded = false;

if isfield(trx,'f2i'),
  trx = rmfield(trx,'f2i');
end

try
  if ~strcmp(varname,'trx'),
    eval(sprintf('%s = trx;',varname));
  end
  if doappend,
    save('-append',matname,varname);
  else
    save(matname,varname);
  end
catch ME
  fprintf('Error saving %s to %s. Error info:\n',varname,matname);
  disp(ME);
  return;
end
succeeded = true;
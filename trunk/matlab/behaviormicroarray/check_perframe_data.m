function isok = check_perframe_data(trx)

if ~isfield(trx,'units'),
  fprintf('Error: Units not defined\n');
  isok = false;
  return;
end
  
fns = fieldnames(trx);
% always ignore these fields
ignorefns = {'id','moviename','firstframe','arena','f2i','nframes','endframe',...
            'pxpermm','fps','matname','annname','sex','type','units'};
fns = setdiff(fns,ignorefns);
ignore = [];

for i = 1:length(fns),
  for j = 1:length(trx),
    expectedlength = trx(j).nframes-1;
    if length(trx(j).(fns{i})) < expectedlength,
      ignore(end+1) = i;
      fprintf('Ignoring field %s\n',fns{i});
      break;
    end
  end
end
fns(ignore) = [];

fnsunits = fieldnames(trx(1).units);
nounits = setdiff(fns,fnsunits);
if ~isempty(nounits),
  fprintf('No units for: ');
  fprintf('%s ',nounits{:});
  fprintf('\n');
  isok = false;
  return;
end

isok = true;
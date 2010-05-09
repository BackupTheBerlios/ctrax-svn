function v = isperframe(in)

v = false;

if isstruct(in),
  trx = in;
elseif ischar(in),
  [trx,matname,succeeded] = load_tracks(in);
  if ~succeeded,
    return;
  end
else
  return;
end
if ~isfield(trx,'units'),
  return;
end

processfcns = perframeprop2processfcn(fieldnames(trx));
v = ~isempty(processfcns);
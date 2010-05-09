function fns = get_perframe_propnames(trx)

fns = fieldnames(trx);
ignorefn = false(1,length(fns));
for i = 1:length(fns),
  for fly = 1:length(trx),
    if length(trx(fly).(fns{i})) < trx(fly).nframes-3,
      ignorefn(i) = true;
      break;
    end
  end
end
fns(ignorefn) = [];

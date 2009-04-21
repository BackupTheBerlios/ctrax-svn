function trx = mindist2fly(trx)

if isfield(trx,'dcenter'),
  for i = 1:length(trx),
    trx(i).mindist2fly = trx(i).dcenter;
  end
  return;
end

nflies = length(trx);
for fly = 1:nflies,
  trx(fly).mindist2fly = nan(1,trx(fly).nframes);
end
for fly1 = 1:nflies-1,
  for fly2 = fly1+1:nflies,
    t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
    t1 = min(trx(fly1).endframe,trx(fly2).endframe);
    if t0 > t1, continue; end
    i0 = trx(fly1).f2i(t0);
    i1 = trx(fly1).f2i(t1);
    j0 = trx(fly2).f2i(t0);
    j1 = trx(fly2).f2i(t1);
    dx = trx(fly1).x(i0:i1) - trx(fly2).x(j0:j1);
    dy = trx(fly1).y(i0:i1) - trx(fly2).y(j0:j1);
    d = sqrt(dx.^2 + dy.^2);
    trx(fly1).mindist2fly(i0:i1) = min(trx(fly1).mindist2fly(i0:i1),d);
    trx(fly2).mindist2fly(j0:j1) = min(trx(fly2).mindist2fly(j0:j1),d);
  end
end
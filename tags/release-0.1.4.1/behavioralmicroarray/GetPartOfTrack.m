function trk1 = GetPartOfTrack(trk0,t0,t1)

fns = fieldnames(trk0);
dontsplit = {'id','moviename','firstframe','arena','f2i','nframes','endframe','offset_start'};
if t0 < trk0.firstframe,
  fprintf('firstframe of trx = %d, requesting frame %d\n',trk0.firstframe,t0);
end
if t1 > trk0.endframe,
  fprintf('lastframe of trx = %d, requesting frame %d\n',trk0.endframe,t1);
end
t0 = min(trk0.endframe,max(t0,trk0.firstframe));
t1 = max(trk0.firstframe,min(t1,trk0.endframe));
i0 = trk0.f2i(t0);
i1 = trk0.f2i(t1);
for i = 1:length(fns)
  fn = fns{i};
  if any(strcmpi(fn,dontsplit)),
    trk1.(fn) = trk0.(fn);
  elseif max(size(trk0.(fn))) < trk0.nframes-1,
    trk1.(fn) = trk0.(fn);
  else
    sz = size(trk0.(fn));
    nd = length(sz);
    if i0 > i1,
      trk1.(fn) = zeros([sz(1:end-1),0]);
      continue;
    end
    [tmp,framedim] = min(abs(sz-trk0.nframes));
    nshort = max(0,trk0.nframes-size(trk0.(fn),framedim));
    args = mat2cell(repmat(':',[1,nd]),1,ones(1,nd));
    args{framedim} = i0:i1-nshort;
    trk1.(fn) = trk0.(fn)(args{:});
  end
end
trk1.firstframe = t0;
trk1.endframe = t1;
trk1.nframes = t1 - t0 + 1;
trk1.f2i = @(f) f-t0+1;

function s = printparams(minparams,maxparams,units)

minfns = fieldnames(minparams);
maxfns = fieldnames(maxparams);
fns = union(minfns,maxfns);
s = cell(1,length(fns));
for i = 1:length(fns),
  fn = fns{i};
  unitscurr = units.(fn);
  
  % minimum value
  if ismember(fn,minfns),
    minv = minparams.(fn);
  else
    minv = -inf;
  end
  if ismember(fn,maxfns),
    maxv = maxparams.(fn);
  else
    maxv = inf;
  end
  
  % convert radians to degrees
  israd = strcmpi('rad',unitscurr.num);
  for j = find(israd),
    unitscurr.num{j} = 'deg';
  end
  n = nnz(israd);
  if n >= 1,
    minv = minv*(180/pi)^n;
    maxv = maxv*(180/pi)^n;
  end
  israd = strcmpi('rad',unitscurr.den);
  n = nnz(israd);
  for j = find(israd),
    unitscurr.num{j} = 'deg';
  end
  if n >= 1,
    minv = minv/(180/pi)^n;
    maxv = maxv/(180/pi)^n;
  end
  
  s{i} = sprintf('%f <= %s <= %f [%s]\n',minv,fn,maxv,unitsstring(unitscurr));
  
end

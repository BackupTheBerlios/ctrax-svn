% i'm not sure if this is exactly correct, but it is close enough

function samples = samplepnoreplace(n,k,w,maxiter)

if ~exist('maxiter','var'),
  maxiter = k;
end;
maxiter = min(maxiter,k)

if isscalar(n),
  n = 1:n;
end;

% if we want more than half of the samples, choose samples that are
% not in the set
isreverse = false;
if k > n/2,
  isreverse = true;
  k = n - k;
end;

samples = [];
nleft = k;
samplesleft = n;
wleft = w;

for iter = 1:maxiter,
  newsamples = unique(randsample(samplesleft,nleft,true,wleft));
  samples = [samples,newsamples];
  nleft = nleft - length(newsamples);
  inds = ismember(samplesleft,newsamples);
  samplesleft(inds) = [];
  wleft(inds) = [];
  if nleft == 0,
    break;
  end;
end;

if nleft > 0,
  newsamples = randsample(samplesleft,nleft,false);
  samples = [samples,newsamples];
end;

if isreverse,
  samples = setdiff(n,samples);
end;
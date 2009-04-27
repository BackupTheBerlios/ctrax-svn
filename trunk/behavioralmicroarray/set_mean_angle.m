function mu = set_mean_angle(x,nsamp)

mu = nan;
if ~exist('nsamp','var'),
  nsamp = min(length(x),20);
end
        
minstd = inf;
for samp = 1:nsamp,
  % choose a random point as the base
  if samp == 1, 
    base = 0; 
  else
    base = randsample(x,1);
  end
  d = modrange(x-base,-pi,pi);
  s = std(d,1);
  if s < minstd,
    minstd = s;
    mu = modrange(base + mean(d),-pi,pi);
  end
end

function out = interval_median(in,t0,t1)

nprops = size(in,2);
nintervals = length(t0);
out = zeros(nintervals,nprops);
for i = 1:nintervals,
  out(i,:) = mean(in(t0(i):t1(i),:),1);
end
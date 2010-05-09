function out = interval_median_angle(in,t0,t1)

nprops = size(in,2);
nintervals = length(t0);
out = zeros(nintervals,nprops);
for i = 1:nintervals,
  tmid = round(mean(t0(i),t1(i)));
  xmid = in(tmid,:);
  n = t1(i)-t0(i)+1;
  dmid = modrange(in(t0(i):t1(i),:) - repmat(xmid,[n,1]),-pi,pi);
  dmidmu = median(dmid,1);
  out(i,:) = modrange(dmidmu+xmid,-pi,pi);
end
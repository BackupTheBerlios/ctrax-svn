function out = invertseg(in,n)

if isempty(in.t1),
  out.t1 = 1;
  out.t2 = n;
  return;
end

out = struct('t1',[],'t2',[]);
out.t1 = in.t2(1:end-1);
out.t2 = in.t1(2:end);
if in.t1(1) > 1,
  out.t1 = [1,out.t1];
  out.t2 = [in.t1(1),out.t2];
end

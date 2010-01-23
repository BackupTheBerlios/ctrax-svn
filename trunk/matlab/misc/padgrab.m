function [y,npadl,npadu] = padgrab(x,padv,varargin)

% parse arguments
usagestring = 'Usage: padgrab(x,padv,l1,u1,l2,u2,...)';
nd1 = ndims(x);
nd = length(varargin)/2;  
isrowvec = nd == 1 && nd1 == 2 && size(x,1) == 1;
if isrowvec,
  varargin = [{1,1},varargin];
  nd = 2;
end
if mod(nd,1) ~= 0 || nd > nd1,
  error(usagestring);
end
l = cell2mat(varargin(1:2:end));
u = cell2mat(varargin(2:2:end));
if any(u < l),
  error(usagestring);
end
sz = size(x);
if nd < nd1,
  x = reshape(x,[sz(1:nd-1),prod(sz(nd:end)),1]);
  sz = size(x);
end

l1 = min(max(l,1),sz);
u1 = min(max(u,1),sz);
a = cell(1,nd);
for i = 1:nd,
  a{i} = l1(i):u1(i);
end
y = x(a{:});
npadl = l1-l;
npadu = u-u1;

if any(npadl > 0),
  y = padarray(y,npadl,padv,'pre');
end
if any(npadu > 0),
  y = padarray(y,npadu,padv,'post');
end

% if dotranspose,
%   y = y';
% end
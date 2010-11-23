% s = structappend(sbase,snew)
% s = structappend(sbase,snew,dim)
function s = structappend(sbase,snew,dim)

if nargin < 3,
  if length(sbase) > 1,
    dim = argmax(size(sbase));
  elseif length(snew) > 1,
    dim = argmax(size(snew));
  else
    dim = 2;
  end
end

if isempty(sbase),
  s = snew;
else
  s = cat(dim,sbase,snew);
end
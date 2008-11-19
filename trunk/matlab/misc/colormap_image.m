% Irgb = colormap_image(I,c)

function Irgb = colormap_image(I,c)

if nargin < 2,
  c = colormap;
end
if ndims(I) ~= 2,
  error('input image must be N x M');
end;
a = min(I(:));
b = max(I(:));
d = b-a;
n = rows(c);
J = round((I - a)/d*(n-1)+1);
Irgb = reshape(c(J(:),:),[rows(I),cols(I),3]);
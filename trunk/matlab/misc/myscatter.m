% h = myscatter(x,y,scale,color,[ncolors],[cm],['filled'],marker,...)
%
%  MYSCATTER Scatter/bubble plot.
%     MYSCATTER(X,Y,S,C) displays colored circles at the locations specified
%     by the vectors X and Y (which must be the same size).  This function
%     is similar to the SCATTER function. Instead of plotting objects as
%     patches, it instead plots them as lines, with a different line for
%     each color instead of a different patch for each point. The colors of
%     the lines will not change when the colormap of the axes is changed.
%     Instead, the colormap must be specified by the input. Some of the
%     differences between MYSCATTER and SCATTER are within the **'s below.
%  
%     S determines the area of each marker (in points^2). S can be a
%     vector the same length a X and Y or a scalar. If S is a scalar, 
%     MATLAB draws all the markers the same size. If S is empty, the
%     default size *12* is used.
%     
%     C determines the colors of the markers. When C is a vector the
%     same length as X and Y, the values in C are linearly mapped
%     to the colors in the *jet* colormap. When C is a 
%     length(X)-by-3 matrix, it directly specifies the colors of the  
%     markers as RGB values. C can also be a color string. See ColorSpec.
%  
%     MYSCATTER(X,Y) draws the markers in the default size and color.
%     MYSCATTER(X,Y,S) draws the markers at the specified sizes (S)
%     with a single color. This type of graph is also known as
%     a bubble plot.
%     *MYSCATTER(X,Y,S,C,NCOLORS) uses NCOLORS different colors in the
%     scatter.*
%     *MYSCATTER(X,Y,S,C,NCOLORS,COLORMAP) specifies the colors to use with
%     the COLORMAP input. If COLORMAP is a function handle, then NCOLORS
%     colors are generated by calling COLORMAP(NCOLORS). Otherwise,
%     COLORMAP must be a NCOLORS x 3 matrix, where COLORMAP(i,:) is the i
%     th color to use.*
%     MYSCATTER(...,M) uses the marker M instead of 'o'.
%     MYSCATTER(...,'filled') fills the markers.
%     *All extra arguments are fed into the plot command, and should come in
%     pairs. *
%  
%     H = MYSCATTER(...) returns handles to the scatter objects created.
%  
%     Use PLOT for single color, single marker size scatter plots.
%  
%     Example
%       load seamount
%       scatter(x,y,5,z)
function varargout = myscatter(varargin)

[x,y,s,c,ncolors,cm,cmisfun,isfilled,marker,rest] = parseargs(varargin{:});

[centers,tmp,idx] = unique(c);
if length(centers) > ncolors,
  [counts,centers,idx] = myhist(c,ncolors);
end

ncenters = length(centers);
if ~cmisfun,
  if ncenters < ncolors,
    idx = round(linspace(1,ncolors,ncenters));
    colors = cm(idx,:);
  else
    colors = cm;
  end
else
  colors = cm(ncenters);
end

h = zeros(1,ncenters);
holdstate = ishold;
hold on;
j = 1;
for i = 1:ncenters,
  idxcurr = idx == i;
  if ~any(idxcurr), continue; end;
  h(j) = plot(x(idx==i),y(idx==i),marker,'color',colors(j,:),rest{:});
  if isfilled,
    set(h(j),'markerfacecolor',colors(j,:));
  end
  j = j + 1;
end

if ~holdstate,
  hold off;
end

if cmisfun
  colormap(cm());
else
  colormap(cm);
end
clim = zeros(1,2);
clim(1) = min(c);
clim(2) = max(c);
if clim(2) <= clim(1),
  clim(2) = clim(1) + .001;
end
caxis(clim);

if nargout > 0,
  varargout{1} = h;
end

function [x,y,s,c,ncolors,cm,cmisfun,isfilled,marker,rest] = parseargs(x,y,s,c,ncolors,cm,varargin)

n = length(x);
if ~exist('s','var') || isempty(s),
  s = 12;
end
if ~exist('c','var') || isempty(c),
  c = ones(1,n);
end
if exist('ncolors','var') && ischar(ncolors),
  varargin{end+1} = ncolors;
  ncolors = [];
end
if exist('cm','var') && ischar(cm),
  varargin{end+1} = cm;
  cm = [];
end
if ~exist('ncolors','var') || isempty(ncolors),
  ncolors = 64;
end
if ~exist('cm','var') || isempty(cm),
  cm = @jet;
end
cmisfun = strcmpi(class(cm),'function_handle');

marker = 'o';
markers = {'.','o','x','+','*','s','d','v','^','<','>','p','h'};
rest = {};
isfilled = false;
i = 1;
while true,
  if i > length(varargin), break; end
  if ischar(varargin{i}) && ismember(varargin{i},markers),
    marker = varargin{i};
    i = i + 1;
  elseif ischar(varargin{i}) && strcmpi(varargin{i},'filled'),
    isfilled = true;
    i = i + 1;
  else
    if length(varargin) == i,
      warning('Error parsing arguments, generic plot args should come in pairs. Ignoring last argument.');
      return;
    end
    rest(end+1:end+2) = varargin(i:i+1);
    i = i + 2;
  end
end

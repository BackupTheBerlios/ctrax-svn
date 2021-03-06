function img = getframe_invisible(img,sz)

if nargin < 1,
  img = get(0,'currentfigure');
  if isempty(img),
    error('No window open');
  end
end
issize = nargin == 2;

if ~ishandle(img)
  error('Input is not a handle');
end

GRAYBORDER = 204;

inputType = getInputType(img);

if strcmp(inputType,'figure'),
  
  hfig = img;
  haxes = get(hfig,'children');
  if length(haxes) > 1,
    error('Figure has multiple children, please input handle of axes instead of figure');
  end

else
  
  haxes = img;
  hfig = get(haxes,'parent');

end

pixelsperinch = get(0,'screenpixelsperInch');
pos =  get(hfig,'position');
% If units are normalized, then re-scale paperposition calculation
% to use pixel values based on the screen size.
if strcmp(get(hfig, 'units'), 'normalized')
  screensize = get(0, 'screensize');
  screenwidth = screensize(3);
  screenheight = screensize(4);
  pos = [(pos(1) * screenwidth + 1) (pos(2) * screenheight + 1) (pos(3) * screenwidth) (pos(4) * screenheight)];
end
set(hfig, 'paperposition', pos./pixelsperinch);
renderer = get(hfig,'renderer');
if strcmp(renderer,'painters')
  renderer = 'opengl';
end
%Turn off warning in case opengl is not supported and
%hardcopy needs to use zbuffer
warnstate = warning('off');
noanimate('save',get(haxes,'parent'));
img = hardcopy(haxes, ['-d' renderer], ['-r' num2str(round(pixelsperinch))]);
if numel(img) == 1,
  fprintf('Could not grab invisible figure. Making visible temporarily.\n');
  set(hfig,'visible','on');
  pause(.1);
  img = hardcopy(haxes, ['-d' renderer], ['-r' num2str(round(pixelsperinch))]);
  pause(.1);
  set(hfig,'visible','off');
end

% crop off the extra border
isallgray = all(img==GRAYBORDER,3);
if issize,
  [nr,nc,tmp] = size(img);
  maxfirstcol = nc - sz(2) + 1;
  maxfirstrow = nr - sz(1) + 1;
  if maxfirstcol <= 0 || maxfirstrow <= 0,
    error('image has size smaller than input sz');
  end
  firstcol = min(find(~all(isallgray,1),1),maxfirstcol);
  firstrow = min(find(~all(isallgray,2),1),maxfirstrow);
  lastcol = firstcol + sz(2) - 1;
  lastrow = firstrow + sz(1) - 1;
else
  firstcol = find(~all(isallgray,1),1);
  firstrow = find(~all(isallgray,2),1);
  lastcol = find(~all(isallgray,1),1,'last');
  lastrow = find(~all(isallgray,2),1,'last');
end
img = img(firstrow:lastrow,firstcol:lastcol,:);

noanimate('restore',get(haxes,'parent'));
warning(warnstate);


function inputType = getInputType(frame)
if isscalar(frame) && ishandle(frame) && (frame > 0)
  inputType = get(frame,'type');
elseif isstruct(frame) & isfield(frame,'cdata')
  inputType = 'movie';
elseif isa(frame,'numeric')
  inputType = 'data';
else
  error('Invalid input argument.  Each frame must be a numeric matrix, a MATLAB movie structure, or a handle to a figure or axis.');
end

function varargout = fixerrorsgui(varargin)
% TRX = FIXERRORSGUI(SEQS,MOVIENAME,TRX,ANNNAME,PARAMS,[LOADNAME])
% SEQS: array of suspicious sequences
% MOVIENAME: name of movie
% TRX: positions of flies over time
% ANNNAME: name of annotation file
% PARAMS: parameters used to determine suspicious sequences
% LOADNAME: previously saved results to restart with

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fixerrorsgui_OpeningFcn, ...
                   'gui_OutputFcn',  @fixerrorsgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fixerrorsgui is made visible.
function fixerrorsgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fixerrorsgui (see VARARGIN)

% Choose default command line output for fixerrorsgui
handles.output = hObject;

% read inputs
handles.seqs = varargin{1};
handles.moviename = varargin{2};
handles.trx = varargin{3};
handles.annname = varargin{4};
handles.params = varargin{5};
if length(varargin) > 5,
  handles.savename = varargin{6};
  if exist(handles.savename,'file')
    load(handles.savename);
    handles.doneseqs = doneseqs;
    handles.moviename = moviename;
    handles.trx = trx;
    handles.seqs = seqs;
    didload = true;
  else
    didload = false;
  end
else
  didload = false;
end

% set up to read movie
[handles.readframe,handles.nframes,handles.fid] = get_readframe_fcn(handles.moviename);

% initialize parameters

% initialize state
isseqleft = false;
for i = 1:length(handles.seqs),
  if ~strcmpi(handles.seqs(i).type,'dummy'),
    isseqleft = true;
    break;
  end
end
if ~isseqleft,
  if ~didload,
    handles.doneseqs = [];
  end
  guidata(hObject,handles);
  msgbox('No suspicious sequences to be corrected. Exiting. ','All Done');
  uiresume(handles.figure1);
  return;
end
handles = SetSeq(handles,i,true);
handles.nselect = 0;
handles.selected = [];
handles.motionobj = [];
handles.plotpath = 'All Flies';
handles.nframesplot = 101;
handles.zoommode = 'sequence';
handles.undolist = {};

handles.bgthresh = 10;
%handles.lighterthanbg = 1;
handles.bgcolor = nan;
[handles.ang_dist_wt,handles.maxjump,bgtype,bgmed,bgmean,...
  tmp,handles.bgthresh] = ...
  read_ann(handles.annname,'ang_dist_wt','max_jump',...
  'bg_algorithm','background_median','background_mean','bg_type',...
  'n_bg_std_thresh_low');
if tmp == 0,
  handles.lighterthanbg = 1;
elseif tmp == 1,
  handles.lighterthanbg = -1;
else
  handles.lighterthanbg = 0;
end
if strcmpi(bgtype,'median'),
  handles.bgmed = bgmed;
else
  handles.bgmed = bgmean;
end

% initialize data structures

% initialize structure to hold seqs that are done
if ~didload,
  handles.doneseqs = [];
end
handles.nflies = length(handles.trx);

% initialize gui

InitializeFrameSlider(handles);
SetFrameNumber(handles);
handles = PlotFirstFrame(handles);
InitializeDisplayPanel(handles);
SetErrorTypes(handles);
handles.bgmed = reshape(handles.bgmed,[handles.nc,handles.nr])';
% Update handles structure
guidata(hObject, handles);

Play(handles,handles.figure1);

% UIWAIT makes fixerrorsgui wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = fixerrorsgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
deletetrx = [];
for fly = 1:length(handles.trx),
  if isdummytrk(handles.trx(fly))
    deletetrx(end+1) = fly;
  end
end
handles.trx(deletetrx) = [];
varargout{1} = handles.trx;
delete(handles.figure1);

function handles = PlotFirstFrame(handles)

axes(handles.mainaxes);
im = handles.readframe(handles.f);
[handles.nr,handles.nc] = size(im);
handles.him = imagesc(im);
colormap gray; axis image; hold on;
zoom reset;

handles.colors = SetFlyColors(handles);

handles.hellipse = zeros(1,handles.nflies);
handles.hcenter = handles.hellipse;
handles.hhead = handles.hellipse;
handles.htail = handles.hellipse;
handles.hleft = handles.hellipse;
handles.hright = handles.hellipse;
handles.htailmarker = handles.hellipse;
handles.hpath = handles.hellipse;
for fly = 1:handles.nflies,
  [handles.hellipse(fly),handles.hcenter(fly),handles.hhead(fly),...
    handles.htail(fly),handles.hleft(fly),handles.hright(fly),...
    handles.htailmarker(fly),handles.hpath(fly)] = ...
    InitFly(handles.colors(fly,:));
  UpdateFlyPathVisible(handles);
  FixUpdateFly(handles,fly);
end

ZoomInOnSeq(handles);

function PlotFrame(handles)

im = handles.readframe(handles.f);
set(handles.him,'cdata',im);
for fly = 1:handles.nflies,
  FixUpdateFly(handles,fly);
  if ~isdummytrk(handles.trx(fly))
    if length(handles.trx(fly).x) ~= handles.trx(fly).nframes || ...
        1 + handles.trx(fly).endframe - handles.trx(fly).firstframe ~= handles.trx(fly).nframes,
      keyboard;
    end
  end
end

function ZoomInOnSeq(handles,seq)

if strcmpi(handles.zoommode,'whole arena'),
  return;
end

if exist('seq','var'),
  handles.seq = seq;
end

BORDER = round(min(handles.nr,handles.nc)/30);

frames = max(min(handles.seq.frames)-10,1):max(handles.seq.frames)+10;
nfliesseq = length(handles.seq.flies);
nframesseq = length(frames);
x0 = nan(nfliesseq,nframesseq);
x1 = x0; y0 = x0; y1 = x0; 
for flyi = 1:nfliesseq,
  fly = handles.seq.flies(flyi);
  for fi = 1:nframesseq,
    f = frames(fi);
    i = handles.trx(fly).f2i(f);
    if isalive(handles.trx(fly),f)
      [x0(flyi,fi),x1(flyi,fi),y0(flyi,fi),y1(flyi,fi)] = ...
        ellipse_to_bounding_box(handles.trx(fly).x(i),handles.trx(fly).y(i),...
        handles.trx(fly).a(i)*2,handles.trx(fly).b(i)*2,handles.trx(fly).theta(i));
    end
  end
end
badidx = isnan(x0);
x0(badidx) = []; y0(badidx) = []; x1(badidx) = []; y1(badidx) = [];

xlim = [min(x0(:))-BORDER,max(x1(:))+BORDER];
xlim = max(min(xlim,handles.nc),1);
ylim = [min(y0(:))-BORDER,max(y1(:))+BORDER];
ylim = max(min(ylim,handles.nr),1);
set(handles.mainaxes,'xlim',xlim,'ylim',ylim);

function v = isdummytrk(trk)

v = any(isnan(trk.x));

function SetFlyVisible(handles,fly,v)

if isdummytrk(handles.trx(fly))
  return;
end

set(handles.hellipse(fly),'visible',v);
set(handles.hcenter(fly),'visible',v);
set(handles.hleft(fly),'visible',v);
set(handles.hright(fly),'visible',v);
set(handles.hhead(fly),'visible',v);
set(handles.htail(fly),'visible',v);
set(handles.htailmarker(fly),'visible',v);
%set(handles.hpath(fly),'visible',v);

function FixUpdateFly(handles,fly)

if isdummytrk(handles.trx(fly)),
  return;
end

ii = handles.trx(fly).f2i(handles.f);

if isalive(handles.trx(fly),handles.f)
  SetFlyVisible(handles,fly,'on');
  i = ii;
else
  SetFlyVisible(handles,fly,'off');
  i = 1;
end

x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
ellipseupdate(handles.hellipse(fly),a,b,x,y,theta);

xleft = x - b*cos(theta+pi/2);
yleft = y - b*sin(theta+pi/2);
xright = x + b*cos(theta+pi/2);
yright = y + b*sin(theta+pi/2);
xhead = x + a*cos(theta);
yhead = y + a*sin(theta);
xtail = x - a*cos(theta);
ytail = y - a*sin(theta);

set(handles.htailmarker(fly),'xdata',[xtail,x],'ydata',[ytail,y]);
set(handles.hleft(fly),'xdata',xleft,'ydata',yleft);
set(handles.hright(fly),'xdata',xright,'ydata',yright);
set(handles.hhead(fly),'xdata',xhead,'ydata',yhead);
set(handles.htail(fly),'xdata',xtail,'ydata',ytail);
set(handles.hcenter(fly),'xdata',x,'ydata',y);

i0 = ii - floor((handles.nframesplot-1)/2);
i1 = ii + handles.nframesplot - 1;
i0 = max(i0,1);
i1 = min(i1,handles.trx(fly).nframes);
set(handles.hpath(fly),'xdata',handles.trx(fly).x(i0:i1),...
  'ydata',handles.trx(fly).y(i0:i1));


function [hellipse,hcenter,hhead,htail,hleft,hright,htailmarker,hpath] = InitFly(color)

hpath = plot(0,0,'.-','color',color,'hittest','off');
htailmarker = plot([0,0],[0,0],'-','color',color,'hittest','off');
hellipse = ellipsedraw(10,10,0,0,0);
set(hellipse,'color',color,'linewidth',2);
set(hellipse,'buttondownfcn','fixerrorsgui(''ellipse_buttondown'',gcbo,[],guidata(gcbo))');
hleft = plot(0,0,'o','markersize',6,'color',color,'markerfacecolor','w');
set(hleft,'buttondownfcn','fixerrorsgui(''left_buttondown'',gcbo,[],guidata(gcbo))');
hright = plot(0,0,'o','markersize',6,'color',color,'markerfacecolor','w');
set(hright,'buttondownfcn','fixerrorsgui(''right_buttondown'',gcbo,[],guidata(gcbo))');
hhead = plot(0,0,'o','markersize',6,'color',color,'markerfacecolor','w');
set(hhead,'buttondownfcn','fixerrorsgui(''head_buttondown'',gcbo,[],guidata(gcbo))');
htail = plot(0,0,'o','markersize',6,'color',color,'markerfacecolor','w');
set(htail,'buttondownfcn','fixerrorsgui(''tail_buttondown'',gcbo,[],guidata(gcbo))');
hcenter = plot(0,0,'o','markersize',6,'color',color,'markerfacecolor','w');
set(hcenter,'buttondownfcn','fixerrorsgui(''center_buttondown'',gcbo,[],guidata(gcbo))');

function UpdateFlyPathVisible(handles)

hObject = handles.plotpathmenu;
contents = get(hObject,'String');
s = contents{get(hObject,'Value')};
handles.plotpath = s;

for fly = 1:handles.nflies,
  if isdummytrk(handles.trx(fly))
    if ishandle(handles.hpath(fly)) && handles.hpath(fly) > 0,
      delete(handles.hpath(fly));
    end
    continue;
  end
  if strcmpi(handles.plotpath,'all flies') || ...
      (strcmpi(handles.plotpath,'seq flies') && ismember(fly,handles.seq.flies)),
    set(handles.hpath(fly),'visible','on');
  else
    set(handles.hpath(fly),'visible','off');
  end
end

function tail_buttondown(hObject,eventdata,handles)

fly = find(handles.htail==hObject);
if isempty(fly), return; end
handles.motionobj = {'tail',fly};
guidata(hObject,handles);

function head_buttondown(hObject,eventdata,handles)

fly = find(handles.hhead==hObject);
if isempty(fly), return; end
handles.motionobj = {'head',fly};
guidata(hObject,handles);

function right_buttondown(hObject,eventdata,handles)

fly = find(handles.hright==hObject);
if isempty(fly), return; end
handles.motionobj = {'right',fly};
guidata(hObject,handles);

function left_buttondown(hObject,eventdata,handles)

fly = find(handles.hleft==hObject);
if isempty(fly), return; end
handles.motionobj = {'left',fly};
guidata(hObject,handles);

function ellipse_buttondown(hObject,eventdata,handles)

fly = find(handles.hellipse==hObject,1);
if isempty(fly), return; end

set(handles.selectedflytext,'string',sprintf('Selected Fly %d',fly));

% are we selecting flies?
if handles.nselect == 0, return; end;

handles = SelectFly(handles,fly);
guidata(hObject,handles);

function center_buttondown(hObject,eventdata,handles)

fly = find(handles.hcenter==hObject);
if isempty(fly), return; end
handles.motionobj = {'center',fly};
guidata(hObject,handles);

function handles = move_center(fly,handles)

tmp = get(handles.mainaxes,'CurrentPoint');

% outside of the axis
if tmp(1,3) ~= 1,
  return;
end

i = handles.trx(fly).f2i(handles.f);
handles.trx(fly).x(i) = tmp(1,1);
handles.trx(fly).y(i) = tmp(1,2);
FixUpdateFly(handles,fly);

function handles = move_head(fly,handles)

tmp = get(handles.mainaxes,'CurrentPoint');
% outside of the axis
if tmp(1,3) ~= 1,
  return;
end
x1 = tmp(1,1);
y1 = tmp(1,2);
i = handles.trx(fly).f2i(handles.f);

x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
theta = handles.trx(fly).theta(i);

x2 = x - a*cos(theta);
y2 = y - a*sin(theta);
x = (x1+x2)/2;
y = (y1+y2)/2;
theta = atan2(y1-y2,x1-x2);
a = sqrt( (x1-x)^2 + (y1-y)^2 )/2;

handles.trx(fly).x(i) = x;
handles.trx(fly).y(i) = y;
handles.trx(fly).a(i) = a;
handles.trx(fly).theta(i) = theta;

FixUpdateFly(handles,fly);

function handles = move_left(fly,handles)

tmp = get(handles.mainaxes,'CurrentPoint');
% outside of the axis
if tmp(1,3) ~= 1,
  return;
end
x3 = tmp(1,1);
y3 = tmp(1,2);
i = handles.trx(fly).f2i(handles.f);

x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);

% compute the distance from this point to the major axis
d = -sin(theta)*(x3 - x) + cos(theta)*(y3 - y);
% compute projection onto minor axis
x3 = x - d * sin(theta);
y3 = y + d * cos(theta);

x4 = x + b*cos(theta+pi/2);
y4 = y + b*sin(theta+pi/2);

x = (x3+x4)/2;
y = (y3+y4)/2;
b = sqrt((x3-x)^2 + (y3-y)^2)/2;

handles.trx(fly).x(i) = x;
handles.trx(fly).y(i) = y;
handles.trx(fly).b(i) = b;

FixUpdateFly(handles,fly);

function handles = move_right(fly,handles)

tmp = get(handles.mainaxes,'CurrentPoint');
% outside of the axis
if tmp(1,3) ~= 1,
  return;
end
x4 = tmp(1,1);
y4 = tmp(1,2);
i = handles.trx(fly).f2i(handles.f);

x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);

% compute the distance from this point to the major axis
d = -sin(theta)*(x4 - x) + cos(theta)*(y4 - y);
% compute projection onto minor axis
x4 = x - d * sin(theta);
y4 = y + d * cos(theta);

x3 = x - b*cos(theta+pi/2);
y3 = y - b*sin(theta+pi/2);

x = (x3+x4)/2;
y = (y3+y4)/2;
b = sqrt((x3-x)^2 + (y3-y)^2)/2;

handles.trx(fly).x(i) = x;
handles.trx(fly).y(i) = y;
handles.trx(fly).b(i) = b;

FixUpdateFly(handles,fly);

function handles = move_tail(fly,handles)

tmp = get(handles.mainaxes,'CurrentPoint');
% outside of the axis
if tmp(1,3) ~= 1,
  return;
end
x2 = tmp(1,1);
y2 = tmp(1,2);
i = handles.trx(fly).f2i(handles.f);

x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
theta = handles.trx(fly).theta(i);

x1 = x + a*cos(theta);
y1 = y + a*sin(theta);
x = (x1+x2)/2;
y = (y1+y2)/2;
theta = atan2(y1-y2,x1-x2);
a = sqrt( (x1-x)^2 + (y1-y)^2 )/2;

handles.trx(fly).x(i) = x;
handles.trx(fly).y(i) = y;
handles.trx(fly).a(i) = a;
handles.trx(fly).theta(i) = theta;

FixUpdateFly(handles,fly);


function handles = SelectFly(handles,fly)

if ismember(fly,handles.selected),
  % set the current fly as unselected
  SetFlySelected(handles,fly,false);
  i = find(handles.selected==fly,1);
  handles.selected(i) = [];
  handles.selected(end+1) = 0;
else
  % set the current fly as selected
  SetFlySelected(handles,fly,true);
  % unselect another fly if necessary
  unselect = handles.selected(end);
  if unselect > 0,
    SetFlySelected(handles,unselect,false);
  end
  % store selected
  handles.selected = [fly,handles.selected(1:end-1)];
end

fprintf('selected = %d\n',handles.selected);

function SetFlySelected(handles,fly,v)

if v,
  set(handles.hellipse(fly),'color',handles.colors(fly,:)*.5+.5,'linewidth',3);
  set(handles.hcenter(fly),'visible','off');
  set(handles.hleft(fly),'visible','off');
  set(handles.hright(fly),'visible','off');
  set(handles.hhead(fly),'visible','off');
  set(handles.htail(fly),'visible','off');
  set(handles.hpath(fly),'linewidth',2);
else
  set(handles.hellipse(fly),'color',handles.colors(fly,:),'linewidth',2);
  set(handles.hcenter(fly),'visible','on');
  set(handles.hleft(fly),'visible','on');
  set(handles.hright(fly),'visible','on');
  set(handles.hhead(fly),'visible','on');
  set(handles.htail(fly),'visible','on');
  set(handles.hpath(fly),'linewidth',1);
end

function colors = SetFlyColors(handles)

colors = jet(handles.nflies);
idx = round(linspace(1,handles.nflies,length(handles.seq.flies)));
for i = 1:length(handles.seq.flies),
  tmp = colors(handles.seq.flies(i),:);
  colors(handles.seq.flies,:) = colors(idx(i),:);
  colors(idx(i),:) = tmp;
end

function InitializeDisplayPanel(handles)

i = find(strcmpi(get(handles.plotpathmenu,'string'),handles.plotpath),1);
set(handles.plotpathmenu,'value',i);
i = find(strcmpi(get(handles.zoommenu,'string'),handles.zoommode),1);
set(handles.zoommenu,'value',i);
set(handles.nframesplotedit,'string',num2str(handles.nframesplot));

function InitializeFrameSlider(handles)

set(handles.frameslider,'max',handles.nframes,'min',1,'sliderstep',[1,20]/(handles.nframes-1));

function SetFrameNumber(handles,hObject)

if nargin < 2,
  hObject = -1;
end

if hObject ~= handles.frameslider,
  set(handles.frameslider,'Value',handles.f);
end
if hObject ~= handles.frameedit,
  set(handles.frameedit,'string',num2str(handles.f));
end
if handles.f < handles.seq.frames(1),
  set(handles.frameofseqtext,'string','Before Sequence','backgroundcolor',[1,0,0],...
    'foregroundcolor',[1,1,1]);
elseif handles.f > handles.seq.frames(end),
  set(handles.frameofseqtext,'string','After Sequence','backgroundcolor',[1,0,0],...
    'foregroundcolor',[1,1,1]);
elseif handles.f == handles.seq.frames(1),
  set(handles.frameofseqtext,'string','Frame of Seq: 1','backgroundcolor',[0,0,1],...
    'foregroundcolor',[1,1,1]);
elseif handles.f == handles.seq.frames(end),
  set(handles.frameofseqtext,'string',...
    sprintf('Frame of Seq: %d',handles.f-handles.seq.frames(1)+1),...
    'backgroundcolor',[1,1,0]/2,'foregroundcolor',[1,1,1]);
else
  set(handles.frameofseqtext,'string',...
    sprintf('Frame of Seq: %d',handles.f-handles.seq.frames(1)+1),...
    'backgroundcolor',[.7,.7,.7],'foregroundcolor',[0,0,0]);
end
i = find(handles.seq.frames == handles.f);
if isempty(i),
  set(handles.suspframetext,'string','Susp: --');
else
  set(handles.suspframetext,'string',sprintf('Susp: %f',handles.seq.suspiciousness(i)));
end

if ~isempty(handles.motionobj),
  if ~isalive(handles.trx(handles,motionobj{2}),handles.f),
    handles.motionobj = [];
  end
end

function v = isalive(track,f)

v = ~isdummytrk(track) && track.firstframe <= f && track.endframe >= f;

% --- Executes on slider movement.
function frameslider_Callback(hObject, eventdata, handles)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.f = round(get(hObject,'value'));
SetFrameNumber(handles,hObject);
PlotFrame(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function frameslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function frameedit_Callback(hObject, eventdata, handles)
% hObject    handle to frameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameedit as text
%        str2double(get(hObject,'String')) returns contents of frameedit as a double
f = str2double(get(hObject,'String'));
if isnan(f),
  set(hObject,'string',num2str(handles.f));
  return;
end
handles.f = round(f);
handles.f = max(f,1);
handles.f = min(f,handles.nframes);
if handles.f ~= f,
  set(hObject,'string',num2str(handles.f));
end
SetFrameNumber(handles,handles.f);
PlotFrame(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function frameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nexterrortypemenu.
function nexterrortypemenu_Callback(hObject, eventdata, handles)
% hObject    handle to nexterrortypemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nexterrortypemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nexterrortypemenu


% --- Executes during object creation, after setting all properties.
function nexterrortypemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nexterrortypemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sortbymenu.
function sortbymenu_Callback(hObject, eventdata, handles)
% hObject    handle to sortbymenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sortbymenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sortbymenu


% --- Executes during object creation, after setting all properties.
function sortbymenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sortbymenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in correctbutton.
function correctbutton_Callback(hObject, eventdata, handles)
% hObject    handle to correctbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% add to undolist
handles.undolist{end+1} = {'correct',handles.seqi,handles.seq};

% move from seqs to doneseqs
if isempty(handles.doneseqs),
  handles.doneseqs = handles.seq;
else
  handles.doneseqs(end+1) = handles.seq;
end
%handles.seqs(handles.seqi) = [];
handles.seqs(handles.seqi).type = 'dummy';

handles.seqs = check_suspicious_sequences(handles.trx,handles.annname,handles.seqs,handles.params{:});

SetErrorTypes(handles);

% quit if this is the last sequence
if strcmpi(get(handles.correctbutton,'string'),'finish')
  msgbox('All suspicious sequences have been corrected. Saving and exiting. ','All Done');
  savebutton_Callback(hObject, [], handles);
  uiresume(handles.figure1);
  return;
end


nexterrortype_type = {'Track Birth','birth',
  'Track Death','death',
  'Match Cost Ambiguity','swap',
  'Large Jump','jump',
  'Large Change in Orientation','orientchange',
  'Velocity & Orient. Mismatch','orientvelmismatch',
  'Large Major Axis','largemajor'};
try
  contents = get(handles.nexterrortypemenu,'string');
  if length(contents) == 1,
    s = contents;
  else
    v = get(handles.nexterrortypemenu,'value');
    if v > length(contents),
      set(handles.nexterrortypemenu,'value',length(contents));
      v = length(contents);
    end
    s = contents{v};
  end
catch
  keyboard;
end

% what is the next type of error
nexttype = nexterrortype_type{strcmpi(nexterrortype_type(:,1),s),2};
flies = [];
frames = [];
susp = [];
idx = [];

% find an error of type nexttype
for i = 1:length(handles.seqs),
  
  % if this is the right type of error
  if strcmpi(handles.seqs(i).type,nexttype),
    
    % store frames, flies, suspiciousness for this seq
    if strcmpi(nexttype,'swap'),
      flies(end+1) = handles.seqs(i).flies(1)*handles.nflies+handles.seqs(i).flies(2);
    else
      flies(end+1) = handles.seqs(i).flies;
    end
    frames(end+1) = handles.seqs(i).frames(1);
    susp(end+1) = max(handles.seqs(i).suspiciousness);
    idx(end+1) = i;
    
  end
  
end

if isempty(flies), keyboard; end

% choose error of this type if there are more than one
contents = get(handles.sortbymenu,'string');
sortby = contents{get(handles.sortbymenu,'value')};
if strcmpi(sortby,'suspiciousness'),
  j = argmax(susp);
  handles = SetSeq(handles,idx(j));
elseif strcmpi(sortby,'frame number'),
  j = argmin(frames);
  handles = SetSeq(handles,idx(j));
elseif strcmpi(sortby,'fly'),
  if strcmpi(handles.seq.type,'swap'),
    currfly = handles.seq.flies(1)*handles.nflies +handles.seq.flies(2);
  else
    currfly = handles.seq.flies;
  end
  issamefly = flies == currfly;
  if any(issamefly),
    nextfly = currfly;
  else
    nextfly = min(flies);
  end
  nextflies = find(flies == nextfly);
  j = nextflies(argmin(frames(nextflies)));
  handles = SetSeq(handles,idx(j));
end

guidata(hObject,handles);

Play(handles,hObject);

function handles = SetSeq(handles,seqi,isfirstframe)

handles.seqi = seqi;
handles.seq = handles.seqs(seqi);
handles.f = handles.seq.frames(1);
handles.nselect = 0;
handles.selected = [];
set(handles.errnumbertext,'string',sprintf('Error: %d/%d',seqi,length(handles.seqs)));
set(handles.seqframestext,'string',sprintf('Frames: %d:%d',handles.seq.frames(1),handles.seq.frames(end)));
set(handles.seqfliestext,'string',['Flies: [',num2str(handles.seq.flies),']']);
set(handles.seqtypetext,'string',sprintf('Type: %s',handles.seq.type));
set(handles.seqsusptext,'string',sprintf('Susp: %f',max(handles.seq.suspiciousness)));
if nargin < 3 || ~isfirstframe,
  SetFrameNumber(handles);
  PlotFrame(handles);
  ZoomInOnSeq(handles);
end


% --- Executes on button press in backbutton.
function backbutton_Callback(hObject, eventdata, handles)
% hObject    handle to backbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global defaultpath;

if ~isfield(handles,'savename') || isempty(handles.savename),
  [path,defaultfilename] = filenamesplit(handles.moviename);
  defaultfilename = splitext(defaultfilename);
  defaultfilename = ['fixed_',defaultfilename,'.mat'];
  if ~isempty(defaultpath),
    defaultfilename = [defaultpath,defaultfilename];
  end
  [filename,defaultpath] = uiputfile('*.mat','Save Fixed Tracks As',defaultfilename);
  handles.savename = [defaultpath,filename];
end

trx = handles.trx;
seqs = handles.seqs;
doneseqs = handles.doneseqs;
moviename = handles.moviename;
seqi = handles.seqi;
params = handles.params;
save(handles.savename,'trx','seqs','doneseqs','moviename','seqi','params');

guidata(hObject,handles);

% --- Executes on button press in quitbutton.
function quitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = questdlg('Save before quitting?','Save?','Yes','No','Yes');
if strcmpi(v,'yes'),
  savebutton_Callback(handles.savebutton, eventdata, handles);
end
uiresume(handles.figure1);

% --- Executes on button press in undobutton.
function undobutton_Callback(hObject, eventdata, handles)
% hObject    handle to undobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in deletedoitbutton.
function deletedoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletedoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('You must first select a fly track to delete. See Delete Track Instructions Panel',...
    'No Fly Selected');
  return;
end

fly = handles.selected;
if handles.f <= handles.trx(fly).firstframe,
  handles.undolist{end+1} = {'delete',handles.f,fly,...
    handles.trx(fly)};
  handles = DeleteFly(handles,fly);
  % remove events involving this fly
  handles = RemoveFlyEvent(handles,fly,-inf,inf);
else  
  handles.undolist{end+1} = {'delete',handles.f,fly,...
    GetPartOfTrack(handles.trx(fly),handles.f,inf)};
  handles.trx(fly) = GetPartOfTrack(handles.trx(fly),1,handles.f-1);
  % remove events involving this fly in the deleted interval
  handles = RemoveFlyEvent(handles,fly,handles.f,inf);
  SetFlySelected(handles,fly,false);
  FixUpdateFly(handles,fly);
end

handles.nselect = 0;
handles.selected = [];
EnablePanel(handles.editpanel,'on');
set(handles.deletepanel,'visible','off');
guidata(hObject,handles);

function trk = GetPartOfTrack(trk,f0,f1)

i0 = trk.f2i(f0);
i1 = trk.f2i(f1);
i0 = max(1,i0);
i1 = min(i1,trk.nframes);
trk.x = trk.x(i0:i1);
trk.y = trk.y(i0:i1);
trk.a = trk.a(i0:i1);
trk.b = trk.b(i0:i1);
trk.theta = trk.theta(i0:i1);
trk.nframes = max(0,i1-i0+1);
trk.firstframe = max(f0,trk.firstframe);
trk.endframe = min(trk.endframe,f1);
trk.f2i = @(f) f - trk.firstframe + 1;

function trk1 = CatTracks(trk1,trk2)

n = trk2.nframes;
trk1.x(end+1:end+n) = trk2.x;
trk1.y(end+1:end+n) = trk2.y;
trk1.a(end+1:end+n) = trk2.a;
trk1.b(end+1:end+n) = trk2.b;
trk1.theta(end+1:end+n) = trk2.theta;
trk1.nframes = trk1.nframes + n;
trk1.endframe = trk1.endframe+n;

% --- Executes on button press in deletecancelbutton.
function deletecancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletecancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected > 0,
  SetFlySelected(handles,handles.selected,false);
end
handles.nselect = 0;
handles.selected = [];
set(handles.deletepanel,'visible','off');
EnablePanel(handles.editpanel,'on');
guidata(hObject,handles);

function EnablePanel(h,v)

children = get(h,'children');
for hchild = children,
  try
    set(hchild,'enable',v);
  catch
  end
end

% --- Executes on button press in renamedoitbutton.
function renamedoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to renamedoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nnz(handles.selected) ~= 2, 
    errordlg('You must first select the two flies two swap. See Swap Identities Instructions Panel',...
    'Bad Selection');
  return;
end

fly1 = handles.selected(1);
fly2 = handles.selected(2);
f = handles.f;

if ~isalive(handles.trx(fly1),f) || ~isalive(handles.trx(fly2),f),
  errordlg('Both flies must be alive in the selected frame.',...
    'Bad Selection');
  return;
end

i1 = handles.trx(fly1).f2i(f);
i2 = handles.trx(fly2).f2i(f);

handles.undolist{end+1} = {'swap',f,[fly1,fly2]};
trk1 = GetPartOfTrack(handles.trx(fly1),f,inf);
trk2 = GetPartOfTrack(handles.trx(fly2),f,inf);
handles.trx(fly1) = GetPartOfTrack(handles.trx(fly1),1,f-1);
handles.trx(fly2) = GetPartOfTrack(handles.trx(fly2),1,f-1);
handles.trx(fly1) = CatTracks(handles.trx(fly1),trk2);
handles.trx(fly2) = CatTracks(handles.trx(fly2),trk1);

handles = FixDeathEvent(handles,fly1);
handles = FixDeathEvent(handles,fly2);

handles = SwapEvents(handles,fly1,fly2,f,inf);

FixUpdateFly(handles,fly1);
FixUpdateFly(handles,fly2);

SetFlySelected(handles,fly1,false);
SetFlySelected(handles,fly2,false);
handles.nselect = 0;
handles.selected = [];
set(handles.swappanel,'visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

function handles = SwapEvents(handles,fly1,fly2,f0,f1)

for i = 1:length(handles.seqs)
  if min(handles.seqs(i).frames) < f0 || max(handles.seqs(i).frames) > f1,
    continue;
  end
  if ismember(fly1,handles.seqs(i).flies) && ~ismember(fly2,handles.seqs(i).flies)
    handles.seqs(i).flies = union(setdiff(handles.seqs(i).flies,fly1),fly2);
  end
  if ismember(fly2,handles.seqs(i).flies) && ~ismember(fly1,handles.seqs(i).flies)
    handles.seqs(i).flies = union(setdiff(handles.seqs(i).flies,fly2),fly1);
  end
end

% --- Executes on button press in renamecancelbutton.
function renamecancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to renamecancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for fly = handles.selected,
  if fly > 0,
    SetFlySelected(handles,fly,false);
  end
end
handles.nselect = 0;
handles.selected = [];
set(handles.swappanel,'visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

% --- Executes on mouse press over axes background.
function mainaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);

% Hint: delete(hObject) closes the figure
%delete(hObject);


% --- Executes on button press in debugbutton.
function debugbutton_Callback(hObject, eventdata, handles)
% hObject    handle to debugbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard;


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'motionobj') || isempty(handles.motionobj), return; end

if strcmpi(handles.motionobj{1},'center'),
  handles = move_center(handles.motionobj{2},handles);
elseif strcmpi(handles.motionobj{1},'head'),
  handles = move_head(handles.motionobj{2},handles);
elseif strcmpi(handles.motionobj{1},'tail'),
  handles = move_tail(handles.motionobj{2},handles);
elseif strcmpi(handles.motionobj{1},'left'),
  handles = move_left(handles.motionobj{2},handles);
elseif strcmpi(handles.motionobj{1},'right'),
  handles = move_right(handles.motionobj{2},handles);
end

guidata(hObject,handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.motionobj = [];
guidata(hObject,handles);

% --- Executes on selection change in editmenu.
function editmenu_Callback(hObject, eventdata, handles)
% hObject    handle to editmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns editmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from editmenu


% --- Executes during object creation, after setting all properties.
function editmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gobutton.
function gobutton_Callback(hObject, eventdata, handles)
% hObject    handle to gobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% what are we doing?
contents = get(handles.editmenu,'string');
s = contents{get(handles.editmenu,'value')};
if strcmpi(s,'delete track...'),
  set(handles.deletepanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 1;
  handles.selected = 0;
elseif strcmpi(s,'interpolate...'),
  set(handles.interpolatepanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 1;
  handles.selected = 0;
  handles.interpolatefirstframe = -1;
  set(handles.interpolatedoitbutton,'enable','off');
elseif strcmpi(s,'connect tracks...'),
  set(handles.connectpanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 1;
  handles.selected = 0;
  handles.connectfirstframe = -1;
  handles.connectfirstfly = -1;
  set(handles.connectdoitbutton,'enable','off');
elseif strcmpi(s,'swap identities...')
  set(handles.swappanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 2;
  handles.selected = zeros(1,2);
elseif strcmpi(s,'extend track...'),
  set(handles.extendpanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 1;
  handles.selected = 0;
  handles.extendfirstframe = -1;
  set(handles.extenddoitbutton,'enable','off');
elseif strcmpi(s,'auto track...'),
  set(handles.autotrackpanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 1;
  handles.selected = 0;
  handles.autotrackframe = -1;
  set(handles.autotrackdoitbutton,'enable','off');
  set(handles.autotracksettingsbutton,'enable','off');
elseif strcmpi(s,'flip orientation...'),
  set(handles.flippanel,'visible','on');
  EnablePanel(handles.editpanel,'off');
  handles.nselect = 1;
  handles.selected = 0;
  handles.flipframe = -1;
  set(handles.flipdoitbutton,'enable','off');
end
guidata(hObject,handles);

% --- Executes on selection change in seekmenu.
function seekmenu_Callback(hObject, eventdata, handles)
% hObject    handle to seekmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns seekmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from seekmenu


% --- Executes during object creation, after setting all properties.
function seekmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seekmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previousbutton.
function previousbutton_Callback(hObject, eventdata, handles)
% hObject    handle to previousbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = get(handles.seekmenu,'value');
contents = get(handles.seekmenu,'string');
s = contents{value};

if strcmpi(s,'birth nearby'),
  
  nextnearbirth = -1;
  nextnearframe = -inf;
  xlim = get(handles.mainaxes,'xlim');
  ylim = get(handles.mainaxes,'ylim');
  for i = 1:length(handles.seqs),
    if ~strcmpi(handles.seqs(i).type,'birth'),
      continue;
    end
    f = handles.seqs(i).frames;
    if f >= handles.f, 
      continue;
    end
    fly = handles.seqs(i).flies;
    j = handles.trx(fly).f2i(f);
    x = handles.trx(fly).x(j);
    y = handles.trx(fly).y(j);
    if x >= xlim(1) && x <= xlim(2) && y >= ylim(1) && y <= ylim(2),
      if nextnearframe < f,
        nextnearbirth = i;
        nextnearframe = f;
      end
    end
  end
  
  if nextnearbirth == -1,
    msgbox('Sorry! There are no fly births in the current axes before the current frame.',...
      'Could Not Find Birth');
    return;
  end
  
  handles.lastframe = handles.f;
  handles.f = nextnearframe;
  SetFrameNumber(handles,hObject);
  PlotFrame(handles);
  
  guidata(hObject,handles);
  
elseif strcmpi(s,'death nearby'),
  
  nextneardeath = -1;
  nextnearframe = -inf;
  xlim = get(handles.mainaxes,'xlim');
  ylim = get(handles.mainaxes,'ylim');
  for i = 1:length(handles.seqs),
    if ~strcmpi(handles.seqs(i).type,'death'),
      continue;
    end
    f = handles.seqs(i).frames;
    if f >= handles.f, 
      continue;
    end
    fly = handles.seqs(i).flies;
    j = handles.trx(fly).f2i(f);
    x = handles.trx(fly).x(j);
    y = handles.trx(fly).y(j);
    if x >= xlim(1) && x <= xlim(2) && y >= ylim(1) && y <= ylim(2),
      if nextnearframe < f,
        nextneardeath = i;
        nextnearframe = f;
      end
    end
  end
  
  if nextneardeath == -1,
    msgbox('Sorry! There are no fly deaths in the current axes before the current frame.',...
      'Could Not Find Death');
    return;
  end
  
  handles.lastframe = handles.f;
  handles.f = nextnearframe;
  SetFrameNumber(handles,hObject);
  PlotFrame(handles);
  
  guidata(hObject,handles);
  
end

% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = get(handles.seekmenu,'value');
contents = get(handles.seekmenu,'string');
s = contents{value};

if strcmpi(s,'birth nearby'),
  
  nextnearbirth = -1;
  nextnearframe = inf;
  xlim = get(handles.mainaxes,'xlim');
  ylim = get(handles.mainaxes,'ylim');
  for i = 1:length(handles.seqs),
    if ~strcmpi(handles.seqs(i).type,'birth'),
      continue;
    end
    f = handles.seqs(i).frames;
    if f <= handles.f, 
      continue;
    end
    fly = handles.seqs(i).flies;
    j = handles.trx(fly).f2i(f);
    x = handles.trx(fly).x(j);
    y = handles.trx(fly).y(j);
    if x >= xlim(1) && x <= xlim(2) && y >= ylim(1) && y <= ylim(2),
      if nextnearframe > f,
        nextnearbirth = i;
        nextnearframe = f;
      end
    end
  end
  
  if nextnearbirth == -1,
    msgbox('Sorry! There are no fly births in the current axes after the current frame.',...
      'Could Not Find Birth');
    return;
  end
  
  handles.lastframe = handles.f;
  handles.f = nextnearframe;
  SetFrameNumber(handles,hObject);
  PlotFrame(handles);
  
  guidata(hObject,handles);
  
elseif strcmpi(s,'death nearby'),
  
  nextneardeath = -1;
  nextnearframe = inf;
  xlim = get(handles.mainaxes,'xlim');
  ylim = get(handles.mainaxes,'ylim');
  for i = 1:length(handles.seqs),
    if ~strcmpi(handles.seqs(i).type,'death'),
      continue;
    end
    f = handles.seqs(i).frames;
    if f <= handles.f, 
      continue;
    end
    fly = handles.seqs(i).flies;
    j = handles.trx(fly).f2i(f);
    x = handles.trx(fly).x(j);
    y = handles.trx(fly).y(j);
    if x >= xlim(1) && x <= xlim(2) && y >= ylim(1) && y <= ylim(2),
      if nextnearframe > f,
        nextneardeath = i;
        nextnearframe = f;
      end
    end
  end
  
  if nextneardeath == -1,
    msgbox('Sorry! There are no fly deaths in the current axes after the current frame.',...
      'Could Not Find Death');
    return;
  end
  
  handles.lastframe = handles.f;
  handles.f = nextnearframe;
  SetFrameNumber(handles,hObject);
  PlotFrame(handles);
  
  guidata(hObject,handles);
  
end


% --- Executes on selection change in plotpathmenu.
function plotpathmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotpathmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotpathmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotpathmenu
UpdateFlyPathVisible(handles);

% --- Executes during object creation, after setting all properties.
function plotpathmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotpathmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nframesplotedit_Callback(hObject, eventdata, handles)
% hObject    handle to nframesplotedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nframesplotedit as text
%        str2double(get(hObject,'String')) returns contents of nframesplotedit as a double
v = str2double(get(hObject,'string'));
if isempty(v), 
  set(hObject,'string',num2str(handles.f));
else
  handles.nframesplot = v;
  for fly = 1:handles.nflies,
    FixUpdateFly(handles,fly);
  end
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nframesplotedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nframesplotedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in zoommenu.
function zoommenu_Callback(hObject, eventdata, handles)
% hObject    handle to zoommenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns zoommenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zoommenu
contents = get(hObject,'String');
s = contents{get(hObject,'Value')};
if strcmpi(handles.zoommode,s),
  return;
end
handles.zoommode = s;
if strcmpi(s,'whole arena'),
  set(handles.mainaxes,'xlim',[1,handles.nc],'ylim',[1,handles.nr]);
else
  ZoomInOnSeq(handles);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function zoommenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoommenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interpolatedoitbutton.
function interpolatedoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to interpolatedoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isalive(handles.trx(handles.interpolatefly),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end

SetFlySelected(handles,handles.interpolatefly,false);
handles.selected = [];

f0 = handles.interpolatefirstframe;
f1 = handles.f;
if f0 > f1,
  tmp = f0; f0 = f1; f1 = tmp;
end
fly = handles.interpolatefly;

% save to undo list
handles.undolist{end+1} = {'interpolate',[f0,f1],fly,...
  GetPartOfTrack(handles.trx(fly),f0,f1)};

% interpolate between f0 and f1
i0 = handles.trx(fly).f2i(f0);
i1 = handles.trx(fly).f2i(f1);
x0 = handles.trx(fly).x(i0);
y0 = handles.trx(fly).y(i0);
a0 = handles.trx(fly).a(i0);
b0 = handles.trx(fly).b(i0);
theta0 = handles.trx(fly).theta(i0);
x1 = handles.trx(fly).x(i1);
y1 = handles.trx(fly).y(i1);
a1 = handles.trx(fly).a(i1);
b1 = handles.trx(fly).b(i1);
theta1 = handles.trx(fly).theta(i1);
nframesinterp = f1-f0+1;
handles.trx(fly).x(i0:i1) = linspace(x0,x1,nframesinterp);
handles.trx(fly).y(i0:i1) = linspace(y0,y1,nframesinterp);
handles.trx(fly).a(i0:i1) = linspace(a0,a1,nframesinterp);
handles.trx(fly).b(i0:i1) = linspace(b0,b1,nframesinterp);

dtheta = modrange(theta1-theta0,-pi,pi);
thetainterp = linspace(0,dtheta,nframesinterp)+theta0;
handles.trx(fly).theta(i0:i1) = modrange(thetainterp,-pi,pi);

delete(handles.hinterpolate);
set(handles.interpolatefirstframebutton,'string','First Frame','Enable','on');
set(handles.interpolatedoitbutton,'enable','off');
set(handles.interpolatepanel','visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

FixUpdateFly(handles,fly);

% --- Executes on button press in cancelinterpolatebutton.
function cancelinterpolatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelinterpolatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hinterpolate') && ishandle(handles.hinterpolate),
  delete(handles.hinterpolate);
end
handles.nselect = 0;
if isfield(handles,'interpolatefly'),
  SetFlySelected(handles,handles.interpolatefly,false);
end
handles.selected = [];
set(handles.interpolatefirstframebutton,'string','First Frame','Enable','on');
set(handles.interpolatedoitbutton,'enable','off');
set(handles.interpolatepanel','visible','off');
EnablePanel(handles.editpanel,'on');

% --- Executes on button press in interpolatefirstframebutton.
function interpolatefirstframebutton_Callback(hObject, eventdata, handles)
% hObject    handle to interpolatefirstframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('Please select fly track to interpolate first.','No Fly Selected');
  return;
end
if ~isalive(handles.trx(handles.selected),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end
handles.interpolatefly = handles.selected;
handles.nselect = 0;
handles.selected = [];
handles.interpolatefirstframe = handles.f;
set(handles.interpolatedoitbutton,'enable','on');
set(handles.interpolatefirstframebutton,'enable','off');
set(handles.interpolatefirstframebutton,'string',sprintf('First = %d',handles.f));

% draw the fly
fly = handles.interpolatefly;
i = handles.trx(fly).f2i(handles.f);
x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
handles.hinterpolate = ellipsedraw(a,b,x,y,theta);
color = handles.colors(fly,:);
set(handles.hinterpolate,'color',color*.75,'linewidth',3,'linestyle','--',...
  'hittest','off');

guidata(hObject,handles);

% --- Executes on button press in interpolatefirstframebutton.
function extendfirstflybutton_Callback(hObject, eventdata, handles)
% hObject    handle to interpolatefirstframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('Please select fly track to extend first.','No Fly Selected');
  return;
end
if ~isalive(handles.trx(handles.selected),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end
handles.extendfly = handles.selected;
handles.nselect = 0;
handles.selected = [];
set(handles.extenddoitbutton,'enable','on');
set(handles.extendfirstflybutton,'enable','off');

% draw the fly
fly = handles.extendfly;
i = handles.trx(fly).f2i(handles.f);
x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
handles.hextend = ellipsedraw(a,b,x,y,theta);
color = handles.colors(fly,:);
set(handles.hextend,'color',color*.75,'linewidth',3,'linestyle','--',...
  'hittest','off');

guidata(hObject,handles);

function UpdateInterpolateFly(handles)

fly = handles.interpolatefly;
i = handles.trx(fly).f2i(handles.interpolatefirstframe);
x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
ellipseupdate(handles.hinterpolate(fly),a,b,x,y,theta);


% --- Executes on button press in connectdoitbutton.
function connectdoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to connectdoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('Please select fly track to connect first.','No Fly Selected');
  return;
end

fly2 = handles.selected;

if ~isalive(handles.trx(fly2),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end

SetFlySelected(handles,handles.connectfirstfly,false);
SetFlySelected(handles,fly2,false);
handles.selected = [];
handles.nselect = 0;

f1 = handles.connectfirstframe;
f2 = handles.f;
fly1 = handles.connectfirstfly;

if f1 > f2,
  tmp = f1; f1 = f2; f2 = tmp;
  tmp = fly1; fly1 = fly2; fly2 = tmp;
end

% save to undo list
handles.undolist{end+1} = {'connect',[f1,f2],[fly1,fly2],...
  [GetPartOfTrack(handles.trx(fly1),f1+1,inf),...
  GetPartOfTrack(handles.trx(fly2),1,f2-1)]};

% interpolate between f1 and f2
i1 = handles.trx(fly1).f2i(f1);
i2 = handles.trx(fly2).f2i(f2);
x1 = handles.trx(fly1).x(i1);
y1 = handles.trx(fly1).y(i1);
a1 = handles.trx(fly1).a(i1);
b1 = handles.trx(fly1).b(i1);
theta1 = handles.trx(fly1).theta(i1);
x2 = handles.trx(fly2).x(i2);
y2 = handles.trx(fly2).y(i2);
a2 = handles.trx(fly2).a(i2);
b2 = handles.trx(fly2).b(i2);
theta2 = handles.trx(fly2).theta(i2);
nframesinterp = f2-f1+1;

xinterp = linspace(x1,x2,nframesinterp);
yinterp = linspace(y1,y2,nframesinterp);
ainterp = linspace(a1,a2,nframesinterp);
binterp = linspace(b1,b2,nframesinterp);
dtheta = modrange(theta2-theta1,-pi,pi);
thetainterp = modrange(linspace(0,dtheta,nframesinterp)+theta1,-pi,pi);

% will we need to cut?
f3 = handles.trx(fly2).endframe;
if f3 < handles.trx(fly1).endframe,
  % if fly1 outlives fly2, then delete all of fly1 after death of fly2
  handles.trx(fly1) = GetPartOfTrack(handles.trx(fly1),1,f3);
  % delete events involving fly1 in frames f3 and after
  handles = RemoveFlyEvent(handles,fly1,f3+1,inf);
elseif f3 > handles.trx(fly1).endframe,
  % we will need to append track
  nappend = f3 - handles.trx(fly1).endframe;
  handles.trx(fly1).x(end+1:end+nappend) = 0;
  handles.trx(fly1).y(end+1:end+nappend) = 0;
  handles.trx(fly1).a(end+1:end+nappend) = 0;
  handles.trx(fly1).b(end+1:end+nappend) = 0;
  handles.trx(fly1).theta(end+1:end+nappend) = 0;
  handles.trx(fly1).nframes = handles.trx(fly1).nframes+nappend;
  handles.trx(fly1).endframe = f3;
end

% copy over the interpolation
idx = i1:handles.trx(fly1).f2i(f2);
handles.trx(fly1).x(idx) = xinterp;
handles.trx(fly1).y(idx) = yinterp;
handles.trx(fly1).a(idx) = ainterp;
handles.trx(fly1).b(idx) = binterp;
handles.trx(fly1).theta(idx) = thetainterp;

% copy over fly2
idx1 = handles.trx(fly1).f2i(f2):handles.trx(fly1).f2i(f3);
idx2 = handles.trx(fly2).f2i(f2):handles.trx(fly2).f2i(f3);
handles.trx(fly1).x(idx1) = handles.trx(fly2).x(idx2);
handles.trx(fly1).y(idx1) = handles.trx(fly2).y(idx2);
handles.trx(fly1).a(idx1) = handles.trx(fly2).a(idx2);
handles.trx(fly1).b(idx1) = handles.trx(fly2).b(idx2);
handles.trx(fly1).theta(idx1) = handles.trx(fly2).theta(idx2);

% delete fly
handles = DeleteFly(handles,fly2);
% replace fly2 with fly1 for frames f2 thru f3
handles = ReplaceFlyEvent(handles,fly2,fly1,f2,f3);
handles = RemoveFlyEvent(handles,fly2,-inf,inf);
handles = FixDeathEvent(handles,fly1);

delete(handles.hconnect);
set(handles.connectfirstflybutton,'string','First Fly','Enable','on');
set(handles.connectdoitbutton,'enable','off');
set(handles.connectpanel','visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

FixUpdateFly(handles,fly1);

function handles = RemoveFlyEvent(handles,fly,f0,f1)

for i = 1:length(handles.seqs)
  if ismember(fly,handles.seqs(i).flies) && f0 <= min(handles.seqs(i).frames) && ...
      f1 >= max(handles.seqs(i).frames)
    handles.seqs(i).type = 'dummy';
  end
end

function handles = ReplaceFlyEvent(handles,fly0,fly1,f0,f1)

%for i = 1:length(handles.seqs)
%  if ismember(fly0,handles.seqs(i).flies) && f0 <= min(handles.seqs(i).frames) && ...
%      f1 >= max(handles.seqs(i).frames)
%    handles.seqs(i).flies = union(setdiff(handles.seqs(i).flies,fly0),fly1);
%  end
%end

function handles = RemoveBirthEvent(handles,fly)

%if handles.trx(fly).firstframe > 1,
%  for i = 1:length(handles.seqs)
%    if strcmpi(handles.seqs(i).type,'birth'),
%      if fly ~= handles.seqs(i).flies,
%        continue;
%      end
%      if isempty(handles.doneseqs),
%        handles.doneseqs = handles.seqs(i);
%      else
%        handles.doneseqs(end+1) = handles.seqs(i);
%      end
%      handles.seqs(i).type = 'dummy';
%    end
%  end
%end

function handles = RemoveDeathEvent(handles,fly)

for i = 1:length(handles.seqs)
  if strcmpi(handles.seqs(i).type,'death'),
    if fly ~= handles.seqs(i).flies,
      continue;
    end
    if isempty(handles.doneseqs),
      handles.doneseqs = handles.seqs(i);
    else
      handles.doneseqs(end+1) = handles.seqs(i);
    end
    handles.seqs(i).type = 'dummy';
  end
end  

function handles = FixDeathEvent(handles,fly)

f = handles.trx(fly).endframe;
if f == handles.nframes,
  handles = RemoveDeathEvent(handles,fly);
else
  for i = 1:length(handles.seqs)
    if ~strcmpi(handles.seqs(i).type,'death'),
      continue;
    end
    if fly ~= handles.seqs(i).flies,
      continue;
    end
    handles.seqs(i).frames = f;
  end
end  

function handles = FixBirthEvent(handles,fly)

f = handles.trx(fly).firstframe;
if f == 1,
  handles = RemoveBirthEvent(handles,fly);
else
  for i = 1:length(handles.seqs)
    if ~strcmpi(handles.seqs(i).type,'birth'),
      continue;
    end
    if fly ~= handles.seqs(i).flies,
      continue;
    end
    handles.seqs(i).frames = f;
  end
end

function handles = DeleteFly(handles,fly)

% find birth and death event for this fly, if it exists
handles = RemoveBirthEvent(handles,fly);
handles = RemoveDeathEvent(handles,fly);

fns = fieldnames(handles.trx(fly));
for i = 1:length(fns),
  fn = fns{i};
  handles.trx(fly).(fn) = nan;
end
delete(handles.hellipse(fly));
delete(handles.hcenter(fly));
delete(handles.hhead(fly));
delete(handles.htail(fly));
delete(handles.hleft(fly));
delete(handles.hright(fly));
delete(handles.htailmarker(fly));
delete(handles.hpath(fly));

% --- Executes on button press in connectcancelbutton.
function connectcancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to connectcancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected > 0,
  SetFlySelected(handles,handles.selected,false);
end
if isfield(handles,'connectfirstfly') && handles.connectfirstfly > 0,
  SetFlySelected(handles,handles.connectfirstfly,false);
end
if isfield(handles,'hconnect') && ishandle(handles.hconnect),
  delete(handles.hconnect);
end
set(handles.connectfirstflybutton,'enable','on','string','First Fly');
handles.nselect = 0;
handles.selected = [];
set(handles.connectpanel,'visible','off');
EnablePanel(handles.editpanel,'on');
guidata(hObject,handles);

% --- Executes on button press in connectfirstflybutton.
function connectfirstflybutton_Callback(hObject, eventdata, handles)
% hObject    handle to connectfirstflybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('Please select fly track to connect first.','No Fly Selected');
  return;
end
if ~isalive(handles.trx(handles.selected),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end
handles.connectfirstfly = handles.selected;
handles.nselect = 1;
handles.selected = 0;
handles.connectfirstframe = handles.f;
set(handles.connectdoitbutton,'enable','on');
set(handles.connectfirstflybutton,'enable','off');
set(handles.connectfirstflybutton,'string',sprintf('First = %d',handles.f));

% draw the fly
fly = handles.connectfirstfly;
i = handles.trx(fly).f2i(handles.f);
x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
handles.hconnect = ellipsedraw(a,b,x,y,theta);
color = handles.colors(fly,:);
set(handles.hconnect,'color',color*.75,'linewidth',3,'linestyle','--',...
  'hittest','off');

guidata(hObject,handles);

function SetErrorTypes(handles)

isbirth = false; isdeath = false;
isswap = false; isjump = false;
isorientchange = false; isorientvelmismatch = false;
islargemajor = false;
for i = 1:length(handles.seqs),
  if strcmpi(handles.seqs(i).type,'dummy'),
    continue;
  end
  eval(sprintf('is%s = true;',handles.seqs(i).type));
end
s = {};
if isbirth,
  s{end+1} = 'Track Birth';
end
if isdeath
  s{end+1} = 'Track Death';
end
if isswap,
  s{end+1} = 'Match Cost Ambiguity';
end
if isjump,
  s{end+1} = 'Large Jump';
end
if isorientchange,
  s{end+1} = 'Large Change in Orientation';
end
if isorientvelmismatch,
  s{end+1} = 'Velocity & Orient. Mismatch';
end
if islargemajor,
 s{end+1} = 'Large Major Axis';
end
content = get(handles.nexterrortypemenu,'string');
v = get(handles.nexterrortypemenu,'value');
if v > length(content),
  set(handles.nexterrortypemenu,'value',length(content));
  v = length(content);
end
sel = content{v};
if isempty(s),
  set(handles.nexterrortypemenu,'string','No more errors','value',1);
  set(handles.correctbutton,'string','Finish');
else
  set(handles.nexterrortypemenu,'string',s);
  set(handles.correctbutton,'string','Correct');
  i = find(strcmpi(sel,s));
  if ~isempty(i),
    set(handles.nexterrortypemenu,'value',i);
  else
    if length(s) >= v,
      set(handles.nexterrortypemenu,'value',min(v,length(s)));
    end
  end
end


% --- Executes on button press in playstopbutton.
function playstopbutton_Callback(hObject, eventdata, handles)
% hObject    handle to playstopbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(hObject,'string'),'play'),
  Play(handles,hObject);
else
  handles.isplaying = false;
  guidata(hObject,handles);
end

function Play(handles,hObject)

handles.isplaying = true;
set(handles.playstopbutton,'string','Stop','backgroundcolor',[.5,0,0]);
guidata(hObject,handles);
f0 = max(1,handles.seq.frames(1)-10);
f1 = min(handles.nframes,handles.seq.frames(end)+10);

for f = f0:f1,
  
  handles.f = f;
  SetFrameNumber(handles);
  PlotFrame(handles);
  drawnow;
  handles = guidata(hObject);

  if ~handles.isplaying,
    break;
  end
  
end

handles.f = f;

if handles.isplaying,
  handles.f = handles.seq.frames(1);
  SetFrameNumber(handles);
  PlotFrame(handles);  
end

handles.isplaying = false;
set(handles.playstopbutton,'string','Play','backgroundcolor',[0,.5,0]);
guidata(hObject,handles);


% --- Executes on button press in extenddoitbutton.
function extenddoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to extenddoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isalive(handles.trx(handles.extendfly),handles.f),
  errordlg('Selected fly is alive in current frame!','Bad Selection');
  return;
end

SetFlySelected(handles,handles.extendfly,false);
handles.selected = [];

f = handles.f;
fly = handles.extendfly;

% save to undo list
handles.undolist{end+1} = {'interpolate',f,fly};

% extend
if f < handles.trx(fly).firstframe,
  n = handles.trx(fly).firstframe - f;
  handles.trx(fly).x = [zeros(1,n),handles.trx(fly).x];
  handles.trx(fly).y = [zeros(1,n),handles.trx(fly).y];
  handles.trx(fly).a = [zeros(1,n),handles.trx(fly).a];
  handles.trx(fly).b = [zeros(1,n),handles.trx(fly).b];
  handles.trx(fly).theta = [zeros(1,n),handles.trx(fly).theta];
  handles.trx(fly).x(1:n) = handles.trx(fly).x(n+1);
  handles.trx(fly).y(1:n) = handles.trx(fly).y(n+1);
  handles.trx(fly).a(1:n) = handles.trx(fly).a(n+1);
  handles.trx(fly).b(1:n) = handles.trx(fly).b(n+1);
  handles.trx(fly).theta(1:n) = handles.trx(fly).theta(n+1);
  handles.trx(fly).firstframe = f;
  handles.trx(fly).f2i = @(f) f - handles.trx(fly).firstframe + 1;
  handles.trx(fly).nframes = length(handles.trx(fly).x);
  % move the death event
  handles = FixDeathEvent(handles,fly);
else
  n = f - handles.trx(fly).endframe;
  handles.trx(fly).x = [handles.trx(fly).x,zeros(1,n)];
  handles.trx(fly).y = [handles.trx(fly).y,zeros(1,n)];
  handles.trx(fly).a = [handles.trx(fly).a,zeros(1,n)];
  handles.trx(fly).b = [handles.trx(fly).b,zeros(1,n)];
  handles.trx(fly).theta = [handles.trx(fly).theta,zeros(1,n)];
  handles.trx(fly).x(end-n+1:end) = handles.trx(fly).x(end-n);
  handles.trx(fly).y(end-n+1:end) = handles.trx(fly).y(end-n);
  handles.trx(fly).a(end-n+1:end) = handles.trx(fly).a(end-n);
  handles.trx(fly).b(end-n+1:end) = handles.trx(fly).b(end-n);
  handles.trx(fly).theta(end-n+1:end) = handles.trx(fly).theta(end-n);
  handles.trx(fly).endframe = f;
  handles.trx(fly).nframes = length(handles.trx(fly).x);
  % move the death event
  handles = FixDeathEvent(handles,fly);
end

delete(handles.hextend);
set(handles.extendfirstflybutton,'Enable','on');
set(handles.extenddoitbutton,'enable','off');
set(handles.extendpanel','visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

FixUpdateFly(handles,fly);

% --- Executes on button press in extendcancelbutton.
function extendcancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to extendcancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hextend') && ishandle(handles.hextend),
  delete(handles.hextend);
end
handles.nselect = 0;
if isfield(handles,'extendfly'),
  SetFlySelected(handles,handles.extendfly,false);
end
handles.selected = [];
set(handles.extendfirstflybutton,'Enable','on');
set(handles.extenddoitbutton,'enable','off');
set(handles.extendpanel','visible','off');
EnablePanel(handles.editpanel,'on');
guidata(hObject,handles);

% --- Executes on button press in autotrackdoitbutton.
function autotrackdoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to autotrackdoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f0 = min(handles.f,handles.autotrackframe);
f1 = max(handles.f,handles.autotrackframe);

SetFlySelected(handles,handles.autotrackfly,false);
handles.selected = [];

fly = handles.autotrackfly;

% save to undo list
handles.undolist{end+1} = {'autotrack',[f0,f1],fly,GetPartOfTrack(handles.trx(fly),f0,f1)};

set(handles.autotrackcancelbutton,'string','Stop');
set(handles.autotrackdoitbutton,'enable','off');
handles.stoptracking = false;

% track
seq.flies = fly;
seq.frames = f0:min(f1,handles.trx(fly).endframe);
if get(handles.showtrackingbutton,'value')
  ZoomInOnSeq(handles,seq);
end
handles.stoptracking = false;
handles = FixTrackFly(fly,f0,f1,handles);
handles = FixDeathEvent(handles,fly);

delete(handles.hautotrack);
set(handles.autotrackcancelbutton,'string','Cancel');
set(handles.autotrackfirstframebutton,'Enable','on');
set(handles.autotrackdoitbutton,'enable','off');
set(handles.autotrackpanel','visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

FixUpdateFly(handles,fly);


% --- Executes on button press in autotrackcancelbutton.
function autotrackcancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to autotrackcancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(handles.autotrackcancelbutton,'string'),'stop')
  handles.stoptracking = true;
else
  if isfield(handles,'hautotrack') && ishandle(handles.hautotrack),
    delete(handles.hautotrack);
  end
  handles.nselect = 0;
  if isfield(handles,'autotrackfly'),
    SetFlySelected(handles,handles.autotrackfly,false);
  end
  handles.selected = [];
  set(handles.autotrackfirstframebutton,'Enable','on');
  set(handles.autotrackdoitbutton,'enable','off');
  set(handles.autotrackpanel','visible','off');
  EnablePanel(handles.editpanel,'on');
end
guidata(hObject,handles);

% --- Executes on button press in autotrackfirstframebutton.
function autotrackfirstframebutton_Callback(hObject, eventdata, handles)
% hObject    handle to autotrackfirstframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('Please select fly track to track first.','No Fly Selected');
  return;
end
if ~isalive(handles.trx(handles.selected),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end
handles.autotrackfly = handles.selected;
handles.autotrackframe = handles.f;

handles.nselect = 0;
handles.selected = [];
set(handles.autotrackdoitbutton,'enable','on');
set(handles.autotrackfirstframebutton,'enable','off');
set(handles.autotracksettingsbutton,'enable','on');
% draw the fly
fly = handles.autotrackfly;
i = handles.trx(fly).f2i(handles.f);
x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
handles.hautotrack = ellipsedraw(a,b,x,y,theta);
color = handles.colors(fly,:);
set(handles.hautotrack,'color',color*.75,'linewidth',3,'linestyle','--',...
  'hittest','off');
handles.bgcurr = handles.bgmed;

guidata(hObject,handles);

% --- Executes on button press in autotracksettingsbutton.
function autotracksettingsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to autotracksettingsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = retrack_settings(handles);
guidata(hObject,handles);

function handles = FixTrackFly(fly,f0,f1,handles)

trk = handles.trx(fly);
se = strel('disk',2);
for f = f0+1:f1

  drawnow;
  handles = guidata(handles.figure1);
  if isfield(handles,'stoptracking') && handles.stoptracking
    break;
  end
  
  i = trk.f2i(f);
  [isfore,dfore,xpred,ypred,thetapred,r0,r1,c0,c1,im] = FixBgSub(fly,f,handles);

  [cc,ncc] = bwlabel(isfore);
  isdeleted = [];
  for fly2 = 1:handles.nflies,
    if fly2 == fly, continue; end
    if ~isalive(handles.trx(fly2),f), continue; end
    i2 = handles.trx(fly2).f2i(f);
    if handles.trx(fly2).x(i2)-(2*handles.trx(fly2).a(i2)+5) > c1 || ...
        handles.trx(fly2).x(i2) + (2*handles.trx(fly2).a(i2)+5)< c0 || ...
        handles.trx(fly2).y(i2) + (2*handles.trx(fly2).a(i2)+5)< r0 || ...
        handles.trx(fly2).y(i2) - (2*handles.trx(fly2).a(i2)+5)> r1,
      continue;
    end
    bw = ellipsepixels([handles.trx(fly2).x(i2),handles.trx(fly2).y(i2),...
      handles.trx(fly2).a(i2)*4,handles.trx(fly2).b(i2)*4,handles.trx(fly2).theta(i2)],...
      [r0,r1,c0,c1]);
    j = 1;
    while true,
      if j > ncc,
        break;
      end

      if ismember(j,isdeleted), 
        j = j + 1;
        continue; 
      end
      fracoverlap = sum(dfore((cc==j) & bw)) / sum(dfore(cc==j));
      if fracoverlap > .75
        isfore(cc==j) = false;
        isdeleted(end+1) = j;
        cc(cc==j) = 0;
      elseif fracoverlap > 0
        bw = imdilate(bw,se);
        isfore(bw) = false;
        cc(bw) = 0;
        tmp = cc == j;
        tmp = imopen(tmp,se);
        %[cctmp,ncctmp] = bwlabel(tmp);
        %if ncctmp > 1
        %  areas = regionprops(cctmp,'area');
        %  areas = getstructarrayfield(areas,'Area');
        %  k = argmax(areas);
        %else
        %  k = 1;
        %end
        %tmp = cctmp==k;
        isfore(cc==j) = false;
        cc(cc==j) = 0;
        cc(tmp) = j;
        isfore(tmp) = true;
        [cctmp,ncctmp] = bwlabel(tmp);
        for k = 2:ncctmp,
          ncc = ncc+1;
          cc(cctmp==k) = ncc;
        end
      end
      j = j + 1;
    end
  end
  % choose the closest connected component
  if ~any(isfore(:)),
    msgbox(sprintf('Frame %d: Could not find the selected fly. Quitting',f));
    return;
  end
  %[cc,ncc] = bwlabel(isfore);
  [tmp1,tmp2,cc] = unique(cc);
  cc = reshape(cc,size(isfore))-1;
  if tmp1(1) == 0
    ncc = length(tmp1)-1;
  end
  xfit = zeros(1,ncc);
  yfit = zeros(1,ncc);
  thetafit = zeros(1,ncc);
  afit = zeros(1,ncc);
  bfit = zeros(1,ncc);
  for j = 1:ncc,
    [y,x] = find(cc==j);
    w = dfore(cc==j);
    [mu,S] = weighted_mean_cov([x,y],w(:));
    xfit(j) = mu(1);
    yfit(j) = mu(2);
    [afit(j),bfit(j),thetafit(j)] = cov2ell(S);
  end
  afit = afit / 2;
  bfit = bfit / 2;
  xfit = xfit + c0 - 1;
  yfit = yfit + r0 - 1;
  if ncc == 1,
    j = 1;
  else
    err = (xpred - xfit).^2 + (ypred - yfit).^2 + handles.ang_dist_wt*(modrange(thetapred - thetafit,-pi/2,pi/2)).^2;
    j = argmin(err);
  end
  trk.x(i) = xfit(j);
  trk.y(i) = yfit(j);
  trk.theta(i) = thetafit(j);
  trk.a(i) = afit(j);
  trk.b(i) = bfit(j);
  

  dtheta = modrange(trk.theta(i)-trk.theta(i-1),-pi/2,pi/2);
  trk.theta(i) = trk.theta(i-1)+dtheta;
  
  handles.trx(fly).x(i) = trk.x(i);
  handles.trx(fly).y(i) = trk.y(i);
  handles.trx(fly).a(i) = trk.a(i);
  handles.trx(fly).b(i) = trk.b(i);
  handles.trx(fly).theta(i) = trk.theta(i);
  handles.f = f;
  if handles.trx(fly).endframe < handles.f
    handles.trx(fly).endframe = f;
  end
  handles.trx(fly).nframes = length(handles.trx(fly).x);
  guidata(handles.figure1,handles);

  if get(handles.showtrackingbutton,'value')
    PlotFrame(handles);
    xlim = get(handles.mainaxes,'xlim');
    ylim = get(handles.mainaxes,'ylim');
    if trk.x(i) < xlim(1) || trk.x(i) > xlim(2) || trk.y(i) < ylim(1) || trk.y(i) > ylim(2)
      seq.frames = [max(f0,f-50),min(f1,f+50)];
      seq.flies = fly;
      ZoomInOnSeq(handles,seq);
    end
  else
    set(handles.frameedit,'string',sprintf('%05d',f));
  end

end


% --- Executes on button press in showtrackingbutton.
function showtrackingbutton_Callback(hObject, eventdata, handles)
% hObject    handle to showtrackingbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showtrackingbutton


% --- Executes on button press in flipdoitbutton.
function flipdoitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to flipdoitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isalive(handles.trx(handles.flipfly),handles.f),
  errordlg('Selected fly is not alive in current frame','Bad Selection');
  return;
end

SetFlySelected(handles,handles.flipfly,false);
handles.selected = [];

f = handles.f;
fly = handles.flipfly;

% save to undo list
handles.undolist{end+1} = {'flip',handles.flipframe,f,fly};

% flip
for f = handles.flipframe:f,
  i = handles.trx(fly).f2i(f);
  handles.trx(fly).theta(i) = modrange(handles.trx(fly).theta(i)+pi,-pi,pi);
end

delete(handles.hflip);
set(handles.flipfirstframebutton,'Enable','on');
set(handles.flipdoitbutton,'enable','off');
set(handles.flippanel','visible','off');
EnablePanel(handles.editpanel,'on');

guidata(hObject,handles);

FixUpdateFly(handles,fly);


% --- Executes on button press in flipcancelbutton.
function flipcancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to flipcancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hflip') && ishandle(handles.hflip),
  delete(handles.hflip);
end
handles.nselect = 0;
if isfield(handles,'flipfly'),
  SetFlySelected(handles,handles.flipfly,false);
end
handles.selected = [];
set(handles.flipfirstframebutton,'Enable','on');
set(handles.flipdoitbutton,'enable','off');
set(handles.flippanel','visible','off');
EnablePanel(handles.editpanel,'on');
guidata(hObject,handles);

% --- Executes on button press in flipfirstframebutton.
function flipfirstframebutton_Callback(hObject, eventdata, handles)
% hObject    handle to flipfirstframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.selected == 0,
  errordlg('Please select fly track to flip first.','No Fly Selected');
  return;
end
if ~isalive(handles.trx(handles.selected),handles.f),
  errordlg('Selected fly is not alive in current frame!','Bad Selection');
  return;
end
handles.flipfly = handles.selected;
handles.flipframe = handles.f;
handles.nselect = 0;
handles.selected = [];
set(handles.flipdoitbutton,'enable','on');
set(handles.flipfirstframebutton,'enable','off');

% draw the fly
fly = handles.flipfly;
i = handles.trx(fly).f2i(handles.f);
x = handles.trx(fly).x(i);
y = handles.trx(fly).y(i);
a = 2*handles.trx(fly).a(i);
b = 2*handles.trx(fly).b(i);
theta = handles.trx(fly).theta(i);
handles.hflip = ellipsedraw(a,b,x,y,theta);
color = handles.colors(fly,:);
set(handles.hflip,'color',color*.75,'linewidth',3,'linestyle','--',...
  'hittest','off');

guidata(hObject,handles);


% --- Executes on button press in printbutton.
function printbutton_Callback(hObject, eventdata, handles)
% hObject    handle to printbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for fly = 1:length(handles.trx)
  fprintf('Track %d: firstframe = %d, endframe = %d, nframes = %d, length(x) = %d\n',...
    fly,handles.trx(fly).firstframe,handles.trx(fly).endframe,handles.trx(fly).nframes,...
    length(handles.trx(fly).x));
end

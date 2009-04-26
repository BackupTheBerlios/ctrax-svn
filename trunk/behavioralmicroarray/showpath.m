function varargout = showpath(varargin)
% SHOWPATH M-file for showpath.fig
%      SHOWPATH, by itself, creates a new SHOWPATH or raises the existing
%      singleton*.
%
%      H = SHOWPATH returns the handle to a new SHOWPATH or the handle to
%      the existing singleton*.
%
%      SHOWPATH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWPATH.M with the given input arguments.
%
%      SHOWPATH('Property','Value',...) creates a new SHOWPATH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showpath_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showpath_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showpath

% Last Modified by GUIDE v2.5 18-Jan-2008 17:03:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @showpath_OpeningFcn, ...
  'gui_OutputFcn',  @showpath_OutputFcn, ...
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


% --- Executes just before showpath is made visible.
function showpath_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showpath (see VARARGIN)

% Choose default command line output for showpath
handles.output = hObject;

% set data
handles.flies = varargin{1};
handles.midframe = varargin{2};
handles.readframe = varargin{3};
handles.trx = varargin{4};
handles.savename = varargin{5};
if ischar(handles.readframe),
  handles.moviename = handles.readframe;
  [handles.readframe,handles.nframes,handles.fid] = get_readframe_fcn(handles.moviename);
end
if ischar(handles.trx),
  handles.matname = handles.trx;
  handles.trx = createdata_perfile('',handles.matname,'',false);
end
handles.nframes = max(getstructarrayfield(handles.trx,'endframe'));
handles.nflies = length(handles.flies);
handles.nfliestot = length(handles.trx);

% set parameters
global pathlength;
if isempty(pathlength),
  pathlength = 2001;
end
pathlength = max(min(pathlength,handles.nframes),1);
handles.startframe = max(1,handles.midframe - ceil((pathlength-1)/2));
handles.endframe = handles.startframe + pathlength - 1;
if handles.endframe > handles.nframes,
  handles.endframe = handles.nframes;
  handles.startframe = handles.endframe - pathlength + 1;
end
pathlength = handles.endframe - handles.startframe + 1;
handles.midframe = round((handles.startframe + handles.endframe)/2);
handles.f = handles.midframe;
handles.zoomwidth = 100;

% setup gui
InitializeFrameSlider(handles);
SetFrameNumber(handles);
InitializeDisplayPanel(handles);

% plot first frame
handles = PlotFirstFrame(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showpath wait for user response (see UIRESUME)
%uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = showpath_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;

function handles = PlotFirstFrame(handles)

global pathlength;

if handles.nflies == 1,
  handles.colors = [1,0,0];
elseif handles.nflies == 2,
  handles.colors = [1,0,0;0,1,1];
else
  handles.colors = lines(handles.nflies);
end
handles.colorsother = jet(handles.nfliestot);
% do we need to swap colors of other flies?
for i = 1:handles.nflies,
  fly = handles.flies(i);
  j = find(all(handles.colorsother == repmat(handles.colors(i,:),[handles.nfliestot,1]),2));
  if j == fly,
    continue;
  end
  if isempty(j)
    handles.colorsother(fly,:) = handles.colors(i,:);
    continue;
  end
  % choose a fly to swap with
  k = randsample(setdiff(1:handles.nfliestot,fly),1);
  handles.colorsother(k,:) = handles.colorsother(fly,:);
  %handles.colorsother(fly,:) = [1,1,1];
end
handles.colorsother(handles.flies,:) = 1;

axes(handles.mainaxes);
im = handles.readframe(handles.f);
[handles.nr,handles.nc] = size(im);
handles.him = imagesc(im);
axis image; colormap gray; hold on;
zoom reset;
ax = axis;
handles.hmarker = zeros(1,handles.nfliestot);
handles.hpathother = zeros(1,handles.nfliestot);
handles.hpathdir = cell(1,handles.nflies);
handles.hpath = zeros(1,handles.nflies);
handles.hpathline = zeros(1,handles.nflies);
handles.htext = cell(1,handles.nflies);

frametext = 20:20:pathlength-1;

for fly = 1:handles.nfliestot,
  if ismember(fly,handles.flies),
    continue;
  else
    handles.hpathother(fly) = plot([0,0],[0,0],'.-');
    set(handles.hpathother(fly),'color',handles.colorsother(fly,:),...
      'hittest','off','visible','off');
  end
end
for i = 1:handles.nflies,
  fly = handles.flies(i);
  x = handles.trx(fly).x(handles.startframe:handles.endframe);
  y = handles.trx(fly).y(handles.startframe:handles.endframe);
  xframe = x(frametext);
  yframe = y(frametext);
  theta = handles.trx(fly).theta(handles.startframe:handles.endframe);
  xtail = x + 2*cos(theta+pi);
  ytail = y + 2*sin(theta+pi);
  x(frametext) = nan;
  y(frametext) = nan;
  xtail(frametext) = nan;
  ytail(frametext) = nan;
  handles.hpathline(i) = plot(x,y,'w-','hittest','off');
  handles.hpathdir{i} = plot([x;xtail],[y;ytail],'-','color',handles.colors(i,:));
  handles.htext{i} = text(xframe,yframe,num2str(frametext'));
  set(handles.htext{i},'color',handles.colors(i,:),'horizontalalignment','center',...
    'clipping','on');
  handles.hpath(i) = plot(x,y,'.','color',handles.colors(i,:),'hittest','off');
  set(handles.hpathdir{i},'buttondownfcn','showpath(''pathclick_callback'',gcbo,[],guidata(gcbo))');
  set(handles.htext{i},'buttondownfcn','showpath(''pathclick_callback'',gcbo,[],guidata(gcbo))');
end

for fly = 1:handles.nfliestot,
  handles.hmarker(fly) = drawflyo(0,0,0,10,10);
  set(handles.hmarker(fly),'color',handles.colorsother(fly,:),'linewidth',2,'hittest','off');
  if ~ismember(fly,handles.flies),
    set(handles.hmarker(fly),'visible','off');
  end
end
axes(handles.zoomlocaxes);
imagesc(im); axis image; colormap gray; axis off;
hold on;
handles.hzoombox = plot(ax([1,1,2,2,1]),ax([3,4,4,3,3]),'r','linewidth',2);

axes(handles.mainaxes);
handles.hzoom = zoom(handles.figure1);
set(handles.hzoom,'ActionPostCallback',@UpdateZoomBox);

UpdatePlot(handles);

function UpdatePlot(handles)

im = handles.readframe(handles.f);
set(handles.him,'cdata',im);
for fly = 1:handles.nfliestot,
  if ~isalive(handles.trx(fly),handles.f)
    set(handles.hmarker(fly),'visible','off');
    if ~ismember(fly,handles.flies),
      set(handles.hpathother(fly),'visible','off');
    end
  else
    i = handles.trx(fly).f2i(handles.f);
    idx = max(i-20,1):min(i+20,handles.trx(fly).nframes);
    if ~ismember(fly,handles.flies)
      set(handles.hpathother(fly),'xdata',handles.trx(fly).x(idx),...
        'ydata',handles.trx(fly).y(idx));
    end
    updatefly(handles.hmarker(fly),handles.trx(fly).x(i),handles.trx(fly).y(i),...
      handles.trx(fly).theta(i),handles.trx(fly).a(i),handles.trx(fly).b(i));
  end
end

if get(handles.zoombutton,'value'),
  ZoomInOnFlies(handles);
end

function ZoomInOnFlies(handles)
  
x0 = inf;
x1 = 0;
y0 = inf;
y1 = 0;
for fly = handles.flies
  if isalive(handles.trx(fly),handles.f)
    i = handles.trx(fly).f2i(handles.f);
    [xa,xb,ya,yb] = ellipse_to_bounding_box(handles.trx(fly).x(i),...
      handles.trx(fly).y(i),handles.trx(fly).a(i)*2,...
      handles.trx(fly).b(i)*2,handles.trx(fly).theta(i));
    x0 = min(x0,xa);
    x1 = max(x1,xb);
    y0 = min(y0,ya);
    y1 = max(y1,yb);
  end
  if ~isinf(x0),
    BORDER = handles.zoomwidth/2;
    x0 = max(x0-BORDER,1); x1 = min(x1+BORDER,handles.nc);
    y0 = max(y0-BORDER,1); y1 = min(y1+BORDER,handles.nr);
    set(handles.mainaxes,'xlim',[x0,x1],'ylim',[y0,y1]);
    ax = [x0,x1,y0,y1];
    set(handles.hzoombox,'xdata',ax([1,1,2,2,1]),'ydata',ax([3,4,4,3,3]));
  end
end

function UpdateZoomBox(obj,event_obj)

handles = guidata(obj);
axes(handles.mainaxes);
set(handles.zoombutton,'value',false);
ax = axis;
set(handles.hzoombox,'xdata',ax([1,1,2,2,1]),'ydata',ax([3,4,4,3,3]));

function InitializeFrameSlider(handles)

set(handles.frameslider,'max',handles.endframe,'min',handles.startframe,...
  'sliderstep',[1,20]/(handles.endframe-handles.startframe-1));

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

function v = isalive(track,f)

v = track.firstframe <= f && track.endframe >= f;

function InitializeDisplayPanel(handles)

set(handles.startframeedit,'string',num2str(handles.startframe));
set(handles.endframeedit,'string',num2str(handles.endframe));
set(handles.plotothersbutton,'value',false);

function startframeedit_Callback(hObject, eventdata, handles)
% hObject    handle to startframeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startframeedit as text
%        str2double(get(hObject,'String')) returns contents of startframeedit as a double
f = str2double(get(hObject,'String'));
if f == handles.startframe,
  return;
end
if isnan(f)
  set(hObject,'string',num2str(handles.f));
  return;
end
if f < 1, f = 1; end;
if f > handles.endframe, f = handles.endframe; end;
handles.startframe = round(f);
if handles.startframe ~= str2double(get(hObject,'String')),
  set(hObject,'string',num2str(handles.f));
end

handles = UpdatePaths(handles);

guidata(hObject,handles);

function handles = UpdatePaths(handles)

global pathlength;
pathlength = handles.endframe - handles.startframe + 1;

for i = 1:handles.nflies,
  delete(handles.hpathdir{i});
  delete(handles.htext{i});
end

frametext = 20:20:pathlength-1;

for i = 1:handles.nflies,
  fly = handles.flies(i);
  x = handles.trx(fly).x(handles.startframe:handles.endframe);
  y = handles.trx(fly).y(handles.startframe:handles.endframe);
  xframe = x(frametext);
  yframe = y(frametext);
  theta = handles.trx(fly).theta(handles.startframe:handles.endframe);
  xtail = x + 2*cos(theta+pi);
  ytail = y + 2*sin(theta+pi);
  x(frametext) = nan;
  y(frametext) = nan;
  xtail(frametext) = nan;
  ytail(frametext) = nan;
  set(handles.hpathline(i),'xdata',x,'ydata',y);
  handles.hpathdir{i} = plot([x;xtail],[y;ytail],'-','color',handles.colors(i,:));
  handles.htext{i} = text(xframe,yframe,num2str(frametext'));
  set(handles.htext{i},'color',handles.colors(i,:),'horizontalalignment','center',...
    'clipping','on');
  set(handles.hpath(i),'xdata',x,'ydata',y);
  set(handles.hpathdir{i},'buttondownfcn','showpath(''pathclick_callback'',gcbo,[],guidata(gcbo))');
  set(handles.htext{i},'buttondownfcn','showpath(''pathclick_callback'',gcbo,[],guidata(gcbo))');

end
InitializeFrameSlider(handles);

if handles.f < handles.startframe
  handles.f = handles.startframe;
  SetFrameNumber(handles);
  UpdatePlot(handles);
end
if handles.f > handles.endframe
  handles.f = handles.endframe;
  SetFrameNumber(handles);
  UpdatePlot(handles);
end

% --- Executes during object creation, after setting all properties.
function startframeedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startframeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endframeedit_Callback(hObject, eventdata, handles)
% hObject    handle to endframeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endframeedit as text
%        str2double(get(hObject,'String')) returns contents of endframeedit as a double
f = str2double(get(hObject,'String'));
if f == handles.endframe,
  return;
end
if isnan(f)
  set(hObject,'string',num2str(handles.f));
  return;
end
if f < handles.startframe, f = handles.startframe; end;
if f > handles.nframes, f = handles.nframes; end;
handles.endframe = round(f);
if handles.endframe ~= str2double(get(hObject,'String')),
  set(hObject,'string',num2str(handles.endframe));
end

handles = UpdatePaths(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function endframeedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endframeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoombutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoombutton
v = get(hObject,'value');
axes(handles.mainaxes);
if v,
  handles.lastzoom = axis;
  guidata(hObject,handles);
  ZoomInOnFlies(handles);
else
  set(handles.mainaxes,'xlim',handles.lastzoom(1:2),'ylim',handles.lastzoom(3:4));
  set(handles.hzoombox,'xdata',handles.lastzoom([1,1,2,2,1]),...
    'ydata',handles.lastzoom([3,4,4,3,3]));
end

% --- Executes on slider movement.
function frameslider_Callback(hObject, eventdata, handles)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
f = get(hObject,'Value');
handles.f = min(handles.nframes,max(1,round(f)));
if f ~= handles.f,
  set(hObject,'value',handles.f);
end

guidata(hObject,handles);
SetFrameNumber(handles,hObject);
UpdatePlot(handles);


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
f = str2double(get(hObject,'string'));
if isnan(f)
  set(hObject,'string',num2str(handles.f));
  return;
end
if f < handles.startframe, f = handles.startframe; end;
if f > handles.endframe, f = handles.endframe; end;
handles.f = round(f);
if handles.f ~= f,
  set(hObject,'string',num2str(handles.f));
end
SetFrameNumber(handles,hObject);
UpdatePlot(handles);
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


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newtrk.flies = handles.flies;
newtrk.startframe = handles.startframe;
newtrk.endframe = handles.endframe;
newtrk.f = handles.f;
axes(handles.mainaxes);
newtrk.ax = axis;

if exist(handles.savename,'file'),
  load(handles.savename);
  chosentrx(end+1) = newtrk;
else
  chosentrx = newtrk;
end
save(handles.savename,'chosentrx');

% --- Executes on button press in discardbutton.
function discardbutton_Callback(hObject, eventdata, handles)
% hObject    handle to discardbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);

% --- Executes on button press in plotothersbutton.
function plotothersbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotothersbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotothersbutton

v = get(hObject,'Value');
if v,
  for fly = 1:handles.nfliestot,
    if isalive(handles.trx(fly),handles.f)
      set(handles.hpathother(fly),'visible','on');
      set(handles.hmarker(fly),'visible','on');
    end
  end
else
  for fly = 1:handles.nfliestot,
    if ~ismember(fly,handles.flies),
      set(handles.hpathother(fly),'visible','off');
      set(handles.hmarker(fly),'visible','off');
    end
  end
end

function pathclick_callback(hObject,eventdata,handles)

j = -1;
for i = 1:handles.nflies,
  if any(handles.hpathdir{i} == hObject),
    j = find(handles.hpathdir{i} == hObject);
    break;
  elseif any(handles.htext{i} == hObject),
    j = find(handles.htext{i} == hObject);
    j = str2num(get(handles.htext{i}(j),'string'));
    break;
  end
end
if j <= 0,
  return;
end

fly = handles.flies(i);
handles.f = handles.startframe+j-1;
SetFrameNumber(handles);
UpdatePlot(handles);
guidata(hObject,handles);

% --- Executes on button press in playpausebutton.
function playpausebutton_Callback(hObject, eventdata, handles)
% hObject    handle to playpausebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(hObject,'string'),'Play'),
  handles.isplaying = true;
  set(handles.playpausebutton,'string','Pause','backgroundcolor',[.5,0,0]);
  guidata(hObject,handles);
  for f = handles.f+1:handles.nframes,
    tic;
    handles.f = f;    
    UpdatePlot(handles);
    SetFrameNumber(handles);
    drawnow;
    handles = guidata(hObject);    
    if ~handles.isplaying,
      break;
    end
  end
  handles.f = f;
end
handles.isplaying = false;
set(handles.playpausebutton,'string','Play','backgroundcolor',[0,.5,0]);
guidata(hObject,handles);

% --- Executes on button press in upzoomwidthbutton.
function upzoomwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to upzoomwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoomwidth = max(round(handles.zoomwidth*1.1),1);
set(handles.zoomwidthtext,'string',sprintf('Zoom Width: %d px',handles.zoomwidth));
guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes on button press in downzoomwidthbutton.
function downzoomwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to downzoomwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoomwidth = max(round(handles.zoomwidth/1.1),1);
set(handles.zoomwidthtext,'string',sprintf('Zoom Width: %d px',handles.zoomwidth));
guidata(hObject,handles);
UpdatePlot(handles);


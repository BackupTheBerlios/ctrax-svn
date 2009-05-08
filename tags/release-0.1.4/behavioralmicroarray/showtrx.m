function varargout = showtrx(varargin)
% paths = showtrx(moviename,matname)
% 
% if moviename/matname are not provided, the file selection dialog will be
% brought up. 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showtrx_OpeningFcn, ...
                   'gui_OutputFcn',  @showtrx_OutputFcn, ...
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


% --- Executes just before showtrx is made visible.
function showtrx_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showtrx (see VARARGIN)

global showtrxlastmovie;
setuppath;

% Choose default command line output for showtrx
handles.output = hObject;

% get movie, matfile names
if length(varargin) >= 1,
  handles.moviename = varargin{1};
else
  helpmsg = 'Choose movie to show';
  if ~isempty(showtrxlastmovie),
    [handles.moviename, handles.pathname] = ...
      uigetfilehelp({'*.fmf';'*.sbfmf';'*.*'},'Choose Movie Name',...
      showtrxlastmovie,'helpmsg',helpmsg);
  else
    [handles.moviename, handles.pathname] = ...
      uigetfilehelp({'*.fmf';'*.sbfmf';'*.*'},'Choose Movie Name','helpmsg',helpmsg);
  end
  if ~ischar(handles.moviename),
    return;
  end
  handles.moviename = [handles.pathname,handles.moviename];
  showtrxlastmovie = handles.moviename;
end
if length(varargin) >= 2,
  handles.matname = varargin{2};
else
  [handles.matname,tmp] = splitext(handles.moviename);
  handles.matname = [handles.matname,'.mat'];
  helpmsg = sprintf('Choose trx file with which to annotate movie %s',handles.moviename);
  [handles.matname,pathname] = uigetfilehelp({'*.mat';'*.*'},'Choose Mat File Name',handles.matname,'helpmsg',helpmsg);
  if ~ischar(handles.matname),
    return;
  end
  handles.matname = [pathname,handles.matname];
end

% read in data
[handles.readframe,handles.nframes,handles.fid] = get_readframe_fcn(handles.moviename);
handles.trx = createdata_perfile('',handles.matname,'',false);
handles.nflies = length(handles.trx);

% set parameters
handles.fps = 10;
handles.pathlength = 50;
handles.isplaying = false;
handles.f = 1;
handles.isselected = false(1,handles.nflies);
handles.iszooming = false;
handles.zoomwidth = 100;
handles.lastflyselected = 1;
handles.propnames = setpropnames(handles.trx);

% setup gui
InitializeFrameSlider(handles);
SetFrameNumber(handles);
InitializeDisplayPanel(handles);

% plot first frame
handles = PlotFirstFrame(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showtrx wait for user response (see UIRESUME)
uiwait(handles.figure1);

function propnames = setpropnames(trx)

propnames = fieldnames(trx);
removeprops = [];
for i = 1:length(propnames),
  if ~isnumeric(trx(1).(propnames{i})) && ...
      ~ischar(trx(1).(propnames{i})),
    removeprops(end+1) = i;
  end
end
propnames(removeprops) = [];


% --- Outputs from this function are returned to the command line.
function varargout = showtrx_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'fid'),
  try
    fclose(handles.fid);
  catch
  end
end
delete(handles.figure1);

function handles = PlotFirstFrame(handles)

handles.colors = jet(handles.nflies);

axes(handles.mainaxes);
im = handles.readframe(handles.f);
[handles.nr,handles.nc] = size(im);
handles.him = imagesc(im);
axis image; colormap gray; hold on;
zoom reset;
ax = axis;
handles.hmarker = zeros(1,handles.nflies);
handles.hpath = zeros(1,handles.nflies);
handles.hselected = zeros(1,handles.nflies);
for fly = 1:handles.nflies,
  handles.hpath(fly) = plot([0,0],[0,0],'.-');
  set(handles.hpath(fly),'color',handles.colors(fly,:));
  set(handles.hpath(fly),'buttondownfcn','showtrx(''selectfly_callback'',gcbo,[],guidata(gcbo))');
end
for fly = 1:handles.nflies,
  handles.hmarker(fly) = drawflyo(0,0,0,10,10);
  set(handles.hmarker(fly),'color',handles.colors(fly,:),'linewidth',2);
  set(handles.hmarker(fly),'buttondownfcn','showtrx(''selectfly_callback'',gcbo,[],guidata(gcbo))');
end
for fly = 1:handles.nflies,
  handles.hselected(fly) = drawflyo(0,0,0,10,10);
  set(handles.hselected(fly),'color','w','linestyle','--','linewidth',2,'hittest','off');
  set(handles.hselected(fly),'visible','off');
end
axes(handles.zoomlocaxes);
imagesc(im); axis image; colormap gray; axis off;
hold on;
handles.hzoombox = plot(ax([1,1,2,2,1]),ax([3,4,4,3,3]),'r','linewidth',2);

axes(handles.mainaxes);
handles.hzoom = zoom(handles.figure1);
set(handles.hzoom,'ActionPostCallback',@UpdateZoomBox);

UpdatePlot(handles,handles.figure1);

function UpdatePlot(handles,hObject)

try
  oldv = get(hObject,'interruptible');
  set(hObject,'interruptible','on');
  im = handles.readframe(handles.f);
  set(hObject,'interruptible',oldv);
  set(handles.him,'cdata',im);
catch
  fprintf('Failed to read frame %d, not updating image\n',handles.f);
end
for fly = 1:handles.nflies,
  if ~isalive(handles.trx(fly),handles.f)
    set(handles.hpath(fly),'visible','off');
    set(handles.hmarker(fly),'visible','off');
    set(handles.hselected(fly),'visible','off');
  else
    i = handles.trx(fly).f2i(handles.f);
    set(handles.hpath(fly),'xdata',handles.trx(fly).x(max(i-handles.pathlength,1):i),...
      'ydata',handles.trx(fly).y(max(i-handles.pathlength,1):i));
    updatefly(handles.hmarker(fly),handles.trx(fly).x(i),handles.trx(fly).y(i),...
      handles.trx(fly).theta(i),handles.trx(fly).a(i),handles.trx(fly).b(i));
    updatefly(handles.hselected(fly),handles.trx(fly).x(i),handles.trx(fly).y(i),...
      handles.trx(fly).theta(i),handles.trx(fly).a(i),handles.trx(fly).b(i));
    % are we only plotting the selected flies?
    if ~handles.isselected(fly) && get(handles.plotselectedbutton,'value'),
      set(handles.hpath(fly),'visible','off');
      set(handles.hmarker(fly),'visible','off');
    else
      set(handles.hpath(fly),'visible','on');
      set(handles.hmarker(fly),'visible','on');
    end
  end
end
updatepropvalue(handles);
if get(handles.zoombutton,'value'),
  ZoomInOnFlies(handles);
end

function ZoomInOnFlies(handles)

if ~any(handles.isselected),
  return;
else
  x0 = inf;
  x1 = 0;
  y0 = inf;
  y1 = 0;
  for fly = find(handles.isselected)
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

function v = isalive(track,f)

v = track.firstframe <= f && track.endframe >= f;

function InitializeDisplayPanel(handles)

set(handles.speedtext,'string',sprintf('Playback Speed: %.1f fps',handles.fps));
set(handles.taillengthedit,'string',num2str(handles.pathlength));
set(handles.fliesselectedtext,'string','[ ]');
set(handles.lastflyselectedtext,'string',num2str(handles.lastflyselected));
set(handles.lastclickedprop,'string',handles.propnames);
updatepropvalue(handles);

function updatepropvalue(handles)

v = get(handles.lastclickedprop,'value');
fn = handles.propnames{v};
if handles.f < handles.trx(handles.lastflyselected).firstframe,
  x = nan;
else
  i = handles.trx(handles.lastflyselected).f2i(handles.f);
  j = min(i,length(handles.trx(handles.lastflyselected).(fn)));
  x = handles.trx(handles.lastflyselected).(fn)(j);
end
if ischar(x),
  set(handles.lastclickedpropvalue,'string',x);
else
  set(handles.lastclickedpropvalue,'string',num2str(x));
end

% --- Executes on button press in playpausebutton.
function playpausebutton_Callback(hObject, eventdata, handles)
% hObject    handle to playpausebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.isplaying,
  handles.isplaying = true;
  set(handles.playpausebutton,'string','Pause','backgroundcolor',[.5,0,0]);
  guidata(hObject,handles);
  for f = handles.f+1:handles.nframes,
    tic;
    handles.f = f;    
    UpdatePlot(handles,hObject);
    SetFrameNumber(handles);
    spf = 1 / handles.fps;
    t = toc;
    pause(max(spf-t,.05));
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
UpdatePlot(handles,hObject);

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
if isnan(f)
  set(hObject,'string',num2str(handles.f));
  return;
end
if f < 1, f = 1; end;
if f > handles.nframes, f = handles.nframes; end;
handles.f = round(f);
if handles.f ~= f,
  set(hObject,'string',num2str(handles.f));
end
SetFrameNumber(handles,hObject);
UpdatePlot(handles,hObject);
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

% --- Executes on button press in increasespeedbutton.
function increasespeedbutton_Callback(hObject, eventdata, handles)
% hObject    handle to increasespeedbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fps = handles.fps * 1.5;
set(handles.speedtext,'string',sprintf('Playback Speed: %.1f fps',handles.fps));
guidata(hObject,handles);

% --- Executes on button press in decreasespeedbutton.
function decreasespeedbutton_Callback(hObject, eventdata, handles)
% hObject    handle to decreasespeedbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fps = handles.fps / 1.5;
set(handles.speedtext,'string',sprintf('Playback Speed: %.1f fps',handles.fps));
guidata(hObject,handles);


function taillengthedit_Callback(hObject, eventdata, handles)
% hObject    handle to taillengthedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of taillengthedit as text
%        str2double(get(hObject,'String')) returns contents of taillengthedit as a double
v = str2double(get(hObject,'String'));
if isnan(v)
  set(hObject,'string',num2str(handles.pathlength));
  return;
end
if v < 0, v = 0; end;
if v > handles.nframes-1, v = handles.nframes-1; end;
handles.pathlength = round(v);
if v ~= handles.pathlength
  set(hObject,'string',num2str(handles.pathlength));
end
UpdatePlot(handles,hObject);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function taillengthedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to taillengthedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function mainaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on mouse press over axes background.
function selectfly_callback(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'interruptible','off');

fly = find(handles.hmarker == hObject);
if isempty(fly)
  fly = find(handles.hpath == hObject);
end
if handles.isselected(fly)
  handles.isselected(fly) = false;
  set(handles.hselected(fly),'visible','off');
else
  handles.isselected(fly) = true;
  set(handles.hselected(fly),'visible','on');
end

if get(handles.zoombutton,'value'),
  ZoomInOnFlies(handles);
end

s = ['[ ',sprintf('%d ',find(handles.isselected)),']'];
set(handles.fliesselectedtext,'string',s);
set(handles.lastflyselectedtext,'string',num2str(fly));
handles.lastflyselected = fly;

guidata(hObject,handles);

updatepropvalue(handles);

set(hObject,'interruptible','on');

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);
uiresume(handles.figure1);


% --- Executes on button press in debugbutton.
function debugbutton_Callback(hObject, eventdata, handles)
% hObject    handle to debugbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of debugbutton
keyboard;


% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoombutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.zoombutton,'value'),
  UpdatePlot(handles,hObject);
end
guidata(hObject,handles);


% --- Executes on button press in plotselectedbutton.
function plotselectedbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotselectedbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotselectedbutton

UpdatePlot(handles,hObject);


% --- Executes on button press in upzoomwidthbutton.
function upzoomwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to upzoomwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoomwidth = max(round(handles.zoomwidth*1.1),1);
set(handles.zoomwidthtext,'string',sprintf('Zoom Width: %d px',handles.zoomwidth));
guidata(hObject,handles);
UpdatePlot(handles,hObject);

% --- Executes on button press in downzoomwidthbutton.
function downzoomwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to downzoomwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoomwidth = max(round(handles.zoomwidth/1.1),1);
set(handles.zoomwidthtext,'string',sprintf('Zoom Width: %d px',handles.zoomwidth));
guidata(hObject,handles);
UpdatePlot(handles,hObject);


% --- Executes on button press in quitbutton.
function quitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);


% --- Executes on selection change in lastclickedprop.
function lastclickedprop_Callback(hObject, eventdata, handles)
% hObject    handle to lastclickedprop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lastclickedprop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lastclickedprop
updatepropvalue(handles);

% --- Executes during object creation, after setting all properties.
function lastclickedprop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lastclickedprop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in computeperframebutton.
function computeperframebutton_Callback(hObject, eventdata, handles)
% hObject    handle to computeperframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[succeeded,savename,trx] = compute_perframe_stats_f('matname',handles.matname);
if succeeded,
  handles.matname = savename;
  handles.trx = trx;
  v = handles.propnames{get(handles.lastclickedprop,'value')};
  handles.propnames = setpropnames(handles.trx);
  i = find(strcmpi(v,handles.propnames),1);
  if isempty(i),
    i = 1;
  end
  set(handles.lastclickedprop,'string',handles.propnames);
  set(handles.lastclickedprop,'value',i);
  updatepropvalue(handles);
end
guidata(hObject,handles);


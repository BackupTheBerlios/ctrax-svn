function varargout = labelmanybehaviors(varargin)
% [starts,ends] = labelmanybehaviors(trx,fly,moviename)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @labelmanybehaviors_OpeningFcn, ...
                   'gui_OutputFcn',  @labelmanybehaviors_OutputFcn, ...
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


% --- Executes just before labelmanybehaviors is made visible.
function labelmanybehaviors_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to labelmanybehaviors (see VARARGIN)

% Choose default command line output for labelmanybehaviors
if length(varargin) < 4,
  error('Usage: [starts,ends] = labelmanybehaviors(trx,fly,moviename,behaviors,...');
end
handles.PATHLENGTH = 500;
handles.output = hObject;
handles.trx = varargin{1};
handles.fly = varargin{2};
handles.moviename = varargin{3};
handles.behaviors = varargin{4};
[handles.segstarts,handles.segends,handles.labels,...
  handles.autosavefilename,...
  handles.scattercommand,...
  handles.behaviorcolors,...
  restart,...
  handles.oldsegstarts,handles.oldsegends,handles.oldlabels,...
  handles.oldbehaviors,handles.oldbehaviorcolors] = ...
  myparse(varargin(5:end),...
  'segstarts',[],'segends',[],'labels',{},...
  'autosavefilename','labeledbehaviors_autosave.mat',...
  'scattercommand','min(trk.velmag,7)',...
  'behaviorcolors',zeros(0,3),...
  'restart',false,...
  'oldsegstarts',[],'oldsegends',[],'oldlabels',{},...
  'oldbehaviors',{},'oldbehaviorcolors',zeros(0,3));

if length(handles.segstarts) ~= length(handles.segends) || ...
    length(handles.segstarts) ~= length(handles.labels),
  error('segstarts, segends, and labels must all be the same size');
end
if restart,
  if ~isempty(handles.segstarts),
    error('segstarts, segends, and labels must be empty if restart set to true');
  end
  tmp = load(handles.autosavefilename);
  if ~tmp.firstframe == handles.trx(handles.fly).firstframe,
    error('Firstframe for autosave data = %d, but for input fly = %d',tmp.firstframe,handles.trx(handles.fly).firstframe);
  end
  if ~tmp.fly == handles.fly,
    error('Fly for autosave data = %d, but for input = %d',tmp.fly,handles.fly);
  end
  handles.segstarts = tmp.starts;
  handles.segends = tmp.ends;
  handles.labels = tmp.labels;
end
if ~isempty(handles.segstarts),
  handles.segstarts = handles.segstarts - handles.trx(handles.fly).firstframe+1;
  handles.segends = handles.segends-handles.trx(handles.fly).firstframe+1;
end
if ~isempty(handles.oldsegstarts),
  handles.oldsegstarts = handles.oldsegstarts - handles.trx(handles.fly).firstframe+1;
  handles.oldsegends = handles.oldsegends-handles.trx(handles.fly).firstframe+1;
end

handles.behaviors = [handles.behaviors,setdiff(unique(handles.labels),handles.behaviors)];
handles.nbehaviors = length(handles.behaviors);
handles.oldbehaviors = [handles.oldbehaviors,setdiff(unique(handles.oldlabels),handles.oldbehaviors)];
handles.noldbehaviors = length(handles.oldbehaviors);
set(handles.scatteredit,'string',handles.scattercommand);
set(handles.edit_pathlength,'String',num2str(handles.PATHLENGTH));

n = handles.nbehaviors - size(handles.behaviorcolors,1);
if n > 0,
  if n <= 6,
    newcolors = lines(n);
  else
    newcolors = jet(n)*.7;
  end
  handles.behaviorcolors = cat(1,handles.behaviorcolors,newcolors);
end
handles.behaviorcurr = handles.behaviors{1};

% store the offset between the left side of the panel position and the
% right side of the figure
set(handles.uipanel_right,'Units','Pixels');
set(handles.uipanel_bottom,'Units','Pixels');
set(handles.figure1,'Units','Pixels');
set(handles.mainaxes,'Units','Pixels');
figpos = get(handles.figure1,'Position');
axespos = get(handles.mainaxes,'Position');
panelpos = get(handles.uipanel_right,'Position');
handles.axesoffright = figpos(3) - (axespos(1)+axespos(3));
handles.axesofftop = figpos(4) - (axespos(2)+axespos(4));
handles.paneloffright = figpos(3) - panelpos(1);
handles.minfigsz = figpos(3:4);

% set parameters
handles.fps = 10;
handles.isplaying = false;
handles.f = handles.trx(handles.fly).firstframe;
handles.iszooming = false;
handles.zoomwidth = 100;
s = strrep(handles.scattercommand,'trk','handles.trx(handles.fly)');
try
  handles.c = eval(s);
catch
  handles.scattercommand = '1';
  handles.c = 1;
end
[handles.readframe,nframes,handles.fid] = get_readframe_fcn(handles.moviename);
handles.nframes = handles.trx(handles.fly).endframe;
handles.isseg = -ones(1,handles.trx(handles.fly).nframes-1);
handles.hstarts = [];
handles.hends = [];
handles.behaviorvalue = 1;
handles.behaviorcurr = handles.behaviors{handles.behaviorvalue};
% setup gui
InitializeFrameSlider(handles);
SetFrameNumber(handles);
InitializeDisplayPanel(handles);

% plot first frame
handles = PlotFirstFrame(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes labelmanybehaviors wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = labelmanybehaviors_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%[starts,ends] = get_interval_ends(handles.isseg>0);
%ends = ends - 1;
starts = handles.segstarts + handles.trx(handles.fly).firstframe-1;
ends = handles.segends + handles.trx(handles.fly).firstframe-1;
varargout{1} = starts;
varargout{2} = ends;
varargout{3} = handles.labels;
if isfield(handles,'fid'),
  try
    fclose(handles.fid);
  catch
  end
end
delete(handles.figure1);

function handles = PlotFirstFrame(handles)

axes(handles.mainaxes);
handles.bd = get(handles.mainaxes,'buttondownfcn');
im = handles.readframe(handles.f);
[handles.nr,handles.nc] = size(im);
handles.him = image(uint8(repmat(im,[1,1,3])));
set(handles.him,'hittest','off');
axis image; hold on;
zoom reset;
ax = axis;

if isempty(handles.behaviorcolors),
  if handles.nbehaviors <= 6,
    handles.behaviorcolors = lines(handles.nbehaviors);
  else
    handles.behaviorcolors = jet(handles.nbehaviors)*.7;
  end
end

if isempty(handles.oldbehaviorcolors),
  handles.oldbehaviorcolors = zeros(handles.noldbehaviors,3);
end

for i = 1:length(handles.oldsegstarts),
  i1 = handles.oldsegstarts(i); i2 = handles.oldsegends(i);
  %i1 = handles.trx(handles.fly).off+(f1);
  %i2 = handles.trx(handles.fly).off+(f2);
  behaviorcurr = handles.oldlabels{i};
  behaviorvalue = find(strcmp(behaviorcurr,handles.oldbehaviors),1);
  plot(handles.trx(handles.fly).x(i1:i2+1),...
    handles.trx(handles.fly).y(i1:i2+1),'-','color',handles.oldbehaviorcolors(behaviorvalue,:),'linewidth',6,'hittest','off');
end


for i = 1:length(handles.segstarts),
  i1 = handles.segstarts(i); i2 = handles.segends(i);
  %i1 = handles.trx(handles.fly).off+(f1);
  %i2 = handles.trx(handles.fly).off+(f2);
  behaviorcurr = handles.labels{i};
  behaviorvalue = find(strcmp(behaviorcurr,handles.behaviors),1);
  handles.isseg(i1:i2) = ...
    plot(handles.trx(handles.fly).x(i1:i2+1),...
    handles.trx(handles.fly).y(i1:i2+1),'-','color',handles.behaviorcolors(behaviorvalue,:),'linewidth',6,'hittest','off');
  handles.hstarts(i) = plot(handles.trx(handles.fly).x(i1),handles.trx(handles.fly).y(i1),'o','color',handles.behaviorcolors(behaviorvalue,:));
  handles.hends(i) = plot(handles.trx(handles.fly).x(i2),handles.trx(handles.fly).y(i2),'s','color',handles.behaviorcolors(behaviorvalue,:));
end
i0 = max(1,handles.f-handles.PATHLENGTH+handles.trx(handles.fly).off);
i1 = min(handles.trx(handles.fly).nframes,handles.f+handles.PATHLENGTH+handles.trx(handles.fly).off);
handles.htrx = plot(handles.trx(handles.fly).x(i0:i1),handles.trx(handles.fly).y(i0:i1),'w-','hittest','off');
[handles.hpath,handles.hpathcenters] = myscatter(handles.trx(handles.fly).x,handles.trx(handles.fly).y,[],handles.c,[],@jet,'.');
caxis([min(handles.c),max(handles.c)]);
colormap jet;
%colorbar;
set(handles.hpath,'hittest','off');
handles.hmarker = drawflyo(0,0,0,10,10);
handles.hcenter = plot(0,0,'o','color','k','hittest','off');
set(handles.hmarker,'color','k','linewidth',2,'hittest','off');
axes(handles.zoomlocaxes);
handles.him_zoomloc = image(uint8(repmat(im,[1,1,3]))); axis image; axis off;
hold on;
handles.hzoombox = plot(ax([1,1,2,2,1]),ax([3,4,4,3,3]),'r','linewidth',2);

axes(handles.mainaxes);
handles.hzoom = zoom(handles.figure1);
set(handles.hzoom,'ActionPostCallback',@UpdateZoomBox);

handles = UpdatePlot(handles);
guidata(handles.figure1,handles);

function handles = UpdatePlot(handles)

im = handles.readframe(handles.f);
set(handles.him,'cdata',uint8(repmat(im,[1,1,3])));
set(handles.him_zoomloc,'cdata',uint8(repmat(im,[1,1,3])));

fly = handles.fly;
i = handles.trx(fly).off+(handles.f);
handles = UpdateBehaviorLabels(handles);

i0 = max(1,handles.f-handles.PATHLENGTH+handles.trx(handles.fly).off);
i1 = min(handles.trx(handles.fly).nframes,handles.f+handles.PATHLENGTH+handles.trx(handles.fly).off);
set(handles.htrx,'xdata',handles.trx(handles.fly).x(i0:i1),...
  'ydata',handles.trx(handles.fly).y(i0:i1));
if i1 > numel(handles.c),
  c = [handles.c(i0:end),repmat(handles.c(end),[1,i1-numel(handles.c)])];
else
  c = handles.c(i0:i1);
end
updatemyscatter(handles.hpath,...
  handles.hpathcenters,...
  handles.trx(handles.fly).x(i0:i1),...
  handles.trx(handles.fly).y(i0:i1),[],...
  c);
for segi = 1:length(handles.segstarts),
  j0 = handles.segstarts(segi); j1 = handles.segends(segi);
  if j1 < i0 || j0 > i1,
    set(handles.isseg(j0),'XData',[],'YData',[]);
  else
    set(handles.isseg(j0),'XData',handles.trx(handles.fly).x(max(i0,j0):min(i1,j1)),...
      'YData',handles.trx(handles.fly).y(max(i0,j0):min(i1,j1)));
  end
  if i0 <= j0 && j0 <= i1,
    set(handles.hstarts(segi),'Visible','on');
  else
    set(handles.hstarts(segi),'Visible','off');
  end
  if i0 <= j1 && j1 <= i1,
    set(handles.hends(segi),'Visible','on');
  else
    set(handles.hends(segi),'Visible','off');
  end

end

updatefly(handles.hmarker,handles.trx(fly).x(i),handles.trx(fly).y(i),...
  handles.trx(fly).theta(i),handles.trx(fly).a(i),handles.trx(fly).b(i));
set(handles.hcenter,'xdata',handles.trx(fly).x(i),'ydata',handles.trx(fly).y(i));

if get(handles.zoombutton,'value'),
  ZoomInOnFlies(handles);
end

set(handles.mainaxes,'buttondownfcn',handles.bd);

function handles = UpdateBehaviorLabels(handles)

fly = handles.fly;
i = handles.trx(fly).off+(handles.f);
if handles.isseg(i) > 0,
  j = find(handles.segends >= i,1);
  if ~strcmp(handles.behaviorcurr,handles.labels{j}),
    handles.behaviorcurr = handles.labels{j};
    handles.behaviorvalue = find(strcmp(handles.behaviors,handles.behaviorcurr));
    set(handles.addlabelbutton,'Value',handles.behaviorvalue,'BackgroundColor',handles.behaviorcolors(handles.behaviorvalue,:));
  end
  set(handles.hmarker,'color',handles.behaviorcolors(handles.behaviorvalue,:));
  set(handles.hcenter,'color',handles.behaviorcolors(handles.behaviorvalue,:));
else
  set(handles.hmarker,'color','k');
  set(handles.hcenter,'color','k');
end


function ZoomInOnFlies(handles)

x0 = inf;
x1 = 0;
y0 = inf;
y1 = 0;
fly = handles.fly;
if isalive(handles.trx(fly),handles.f)
  i = handles.trx(fly).off+(handles.f);
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

function UpdateZoomBox(obj,event_obj)

handles = guidata(obj);
axes(handles.mainaxes);
set(handles.mainaxes,'buttondownfcn',handles.bd);
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
if handles.f >= handles.trx(handles.fly).endframe,
  set(handles.ductr_text,'string',sprintf('du_ctr: nan'));
  set(handles.dvctr_text,'string',sprintf('dv_ctr: nan'));
  set(handles.dtheta_text,'string',sprintf('dtheta: nan'));
  set(handles.ducor_text,'string',sprintf('du_cor: nan'));
  set(handles.dvcor_text,'string',sprintf('dv_cor: nan'));
  set(handles.corfrac_text,'string',sprintf('cor: nan'));
else
  i = handles.trx(handles.fly).off+(handles.f);
  set(handles.scattervaltext,'string',sprintf('%f',handles.c(i)));
  set(handles.ductr_text,'string',sprintf('du_ctr: %f',handles.trx(handles.fly).du_ctr(i)));
  set(handles.dvctr_text,'string',sprintf('dv_ctr: %f',handles.trx(handles.fly).dv_ctr(i)));
  set(handles.dtheta_text,'string',sprintf('dtheta: %f',180/pi*handles.trx(handles.fly).dtheta(i)));
  set(handles.ducor_text,'string',sprintf('du_cor: %f',handles.trx(handles.fly).du_cor(i)));
  set(handles.dvcor_text,'string',sprintf('dv_cor: %f',handles.trx(handles.fly).dv_cor(i)));
  set(handles.corfrac_text,'string',sprintf('cor: %f',handles.trx(handles.fly).corfrac_maj(i)));
end
function v = isalive(track,f)

v = track.firstframe <= f && track.endframe >= f;

function InitializeDisplayPanel(handles)

set(handles.speedtext,'string',sprintf('Playback Speed: %.1f fps',handles.fps));
set(handles.addlabelbutton,'String',handles.behaviors,'Value',handles.behaviorvalue,'BackgroundColor',handles.behaviorcolors(handles.behaviorvalue,:));

% --- Executes on button press in playpausebutton.
function playpausebutton_Callback(hObject, eventdata, handles)
% hObject    handle to playpausebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.isplaying,
  handles.isplaying = true;
  set(handles.playpausebutton,'string','Pause','backgroundcolor',[.5,0,0]);
  guidata(hObject,handles);
  for f = handles.f+1:handles.trx(handles.fly).endframe,
    tic;
    handles.f = f;    
    handles = UpdatePlot(handles);
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
handles.f = min(handles.trx(handles.fly).endframe,max(handles.trx(handles.fly).firstframe,round(f)));
if f ~= handles.f,
  set(hObject,'value',handles.f);
end

SetFrameNumber(handles,hObject);
handles = UpdatePlot(handles);
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
handles = UpdatePlot(handles);
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

% --- Executes on mouse press over axes background.
function mainaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(handles.hzoom,'enable'),'on'), return; end
pt = get(handles.mainaxes,'currentpoint');
x = pt(1,1); y = pt(1,2);
i0 = max(1,handles.f-handles.PATHLENGTH+handles.trx(handles.fly).off);
i1 = min(handles.trx(handles.fly).nframes,handles.f+handles.PATHLENGTH+handles.trx(handles.fly).off);
d = sqrt((handles.trx(handles.fly).x(i0:i1)-x).^2 + ...
  (handles.trx(handles.fly).y(i0:i1)-y).^2);
[mind,t] = min(d);
t = t + i0 - 1;
ax = axis;
maxd = max(ax(2)-ax(1),ax(4)-ax(3))/20;
if mind < maxd,
  handles.f = t+handles.trx(handles.fly).firstframe-1;
  SetFrameNumber(handles,hObject);
  handles = UpdatePlot(handles);
  guidata(hObject,handles);
end

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
  handles = UpdatePlot(handles);
end
guidata(hObject,handles);

% --- Executes on button press in upzoomwidthbutton.
function upzoomwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to upzoomwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoomwidth = max(round(handles.zoomwidth*1.1),1);
set(handles.zoomwidthtext,'string',sprintf('Zoom Width: %d px',handles.zoomwidth));
handles = UpdatePlot(handles);
guidata(hObject,handles);

% --- Executes on button press in downzoomwidthbutton.
function downzoomwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to downzoomwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoomwidth = max(round(handles.zoomwidth/1.1),1);
set(handles.zoomwidthtext,'string',sprintf('Zoom Width: %d px',handles.zoomwidth));
handles = UpdatePlot(handles);
guidata(hObject,handles);


% --- Executes on button press in addstartbutton.
function addstartbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addstartbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f1 = handles.f;
i1 = handles.trx(handles.fly).off+(f1);
if handles.isseg(i1)>0 || f1 > handles.trx(handles.fly).endframe || ...
    f1 < handles.trx(handles.fly).firstframe,
  return;
end
i2 = find(handles.isseg(i1+1:end)>0,1,'first');
if isempty(i2),
  i2 = handles.trx(handles.fly).nframes-1;
else
  i2 = i2 + i1 - 1;
end

handles.currentstarti = i1;
%i1 = handles.trx(handles.fly).off+(f1);
%i2 = handles.trx(handles.fly).off+(f2);
handles.isseg(i1:i2) = plot(handles.trx(handles.fly).x(i1:i2+1),...
  handles.trx(handles.fly).y(i1:i2+1),'-','color',handles.behaviorcolors(handles.behaviorvalue,:),'linewidth',6,'hittest','off');
chil = get(handles.mainaxes,'children');
chil = [chil(2:end-1);chil(1);chil(end)];
set(handles.mainaxes,'children',chil);
j = find(handles.segstarts >= i1,1);
hstart = plot(handles.trx(handles.fly).x(i1),handles.trx(handles.fly).y(i1),'o','color',handles.behaviorcolors(handles.behaviorvalue,:));
hend = plot(handles.trx(handles.fly).x(i2),handles.trx(handles.fly).y(i2),'s','color',handles.behaviorcolors(handles.behaviorvalue,:));
if isempty(j),
  j = length(handles.segstarts)+1;
end
handles.segstarts = [handles.segstarts(1:j-1),i1,handles.segstarts(j:end)];
handles.segends = [handles.segends(1:j-1),i2,handles.segends(j:end)];
handles.labels = [handles.labels(1:j-1),{handles.behaviorcurr},handles.labels(j:end)];
handles.hstarts = [handles.hstarts(1:j-1),hstart,handles.hstarts(j:end)];
handles.hends = [handles.hends(1:j-1),hend,handles.hends(j:end)];

handles = UpdatePlot(handles);

guidata(hObject,handles);

% --- Executes on button press in addendbutton.
function addendbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addendbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% end last interval at f1
f1 = handles.f;
i1 = handles.trx(handles.fly).off+(f1);
if i1 == handles.trx(handles.fly).nframes,
  return;
end
if handles.isseg(i1)<=0 || f1 > handles.trx(handles.fly).endframe || ...
    f1 < handles.trx(handles.fly).firstframe,
  return;
end

j = find(handles.segstarts <= i1 & handles.segends >= i1,1);

% reset until f2 to -1
i2 = handles.segends(j);
handles.isseg(i1+1:i2) = -1;

% start of this interval
i0 = handles.segstarts(j);
%i0 = find(handles.isseg(1:i1-1)<=0,1,'last');
%if isempty(i0),
%  i0 = 1;
%  %f0 = handles.trx(handles.fly).firstframe;
%else
%  i0 = i0 + 1;
%  %f0 = handles.trx(handles.fly).firstframe + i0 - 1;
%end

% set plot data
set(handles.isseg(i1),'xdata',handles.trx(handles.fly).x(i0:i1+1),...
  'ydata',handles.trx(handles.fly).y(i0:i1+1));
set(handles.hends(j),'xdata',handles.trx(handles.fly).x(i1+1),...
  'ydata',handles.trx(handles.fly).y(i1+1));
handles.segends(j) = i1;

handles = UpdateBehaviorLabels(handles);

guidata(hObject,handles);

starts = handles.segstarts + handles.trx(handles.fly).firstframe-1;
ends = handles.segends + handles.trx(handles.fly).firstframe-1;
labels = handles.labels;
fly = handles.fly;
firstframe = handles.trx(handles.fly).firstframe;

save(handles.autosavefilename,'starts','ends','labels','firstframe','fly');
fprintf('Autosaving to file %s\n',handles.autosavefilename);

function scatteredit_Callback(hObject, eventdata, handles)
% hObject    handle to scatteredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scatteredit as text
%        str2double(get(hObject,'String')) returns contents of scatteredit as a double

s = get(hObject,'string');
try
  s = strrep(s,'trk','handles.trx(handles.fly)');
  c = eval(s);
catch
  warning('Could not evaluate scattercommand, reverting.');
  set(handles.scatteredit,'string',handles.scattercommand);
  return;
end
handles.c = c;
caxis([min(handles.c),max(handles.c)]);
delete(handles.hpath(ishandle(handles.hpath)&handles.hpath>0));
handles.hpath = myscatter(handles.trx(handles.fly).x,handles.trx(handles.fly).y,[],handles.c,[],@jet,'.');
set(handles.hpath,'hittest','off');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function scatteredit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scatteredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deletesegbutton.
function deletesegbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletesegbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f = handles.f;
off = handles.trx(handles.fly).firstframe-1;
i = handles.trx(handles.fly).off+(f);
if handles.isseg(i) <= 0,
  fprintf('No label at this frame. Not deleting\n');
  return;
end
fidx = handles.trx(handles.fly).off+(f);
tmp = find(handles.segstarts<=fidx);
[i0,j] = max(handles.segstarts(tmp));
j = tmp(j);
%i0 = find(handles.isseg(1:i-1)<=0,1,'last');
%if isempty(i0),
%  i0 = 1;
%else
%  i0 = i0 - 1;
%end
i1 = handles.segends(j);
%i1 = find(handles.isseg(i+1:end)<=0,1,'first');
%if isempty(i1),
%  i1 = handles.nframes;
%else
%  i1 = i1 + i - 1;
%end
delete(handles.isseg(i));
delete(handles.hstarts(j));
delete(handles.hends(j));
handles.segstarts(j) = [];
handles.segends(j) = [];
handles.labels(j) = [];
handles.hstarts(j) = [];
handles.hends(j) = [];
handles.isseg(i0:i1) = -1;

handles = UpdateBehaviorLabels(handles);

guidata(hObject,handles);


% --- Executes on button press in addlabelbutton.
function addlabelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addlabelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviorvalue = get(hObject,'value');
handles.behaviorcurr = handles.behaviors{handles.behaviorvalue};
set(hObject,'BackgroundColor',handles.behaviorcolors(handles.behaviorvalue,:));

f = handles.f;
i = handles.trx(handles.fly).off+(f);
if handles.isseg(i) > 0,
  % switch the label
  j = find(handles.segends >= i,1);
  handles.labels{j} = handles.behaviorcurr;
  handles.behaviorvalue = find(strcmp(handles.behaviorcurr,handles.behaviors),1);
  set(handles.isseg(i),'color',handles.behaviorcolors(handles.behaviorvalue,:));
end

guidata(hObject,handles);


% --- Executes on button press in prevlabel_button.
function prevlabel_button_Callback(hObject, eventdata, handles)
% hObject    handle to prevlabel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = handles.trx(handles.fly).off+(handles.f);
j = find(handles.segends < i,1,'last');
if isempty(j), return; end
j = handles.segends(j);
handles.f = handles.trx(handles.fly).firstframe + j - 1;
SetFrameNumber(handles,hObject);
handles = UpdatePlot(handles);
guidata(hObject,handles);

% --- Executes on button press in nextlabel_button.
function nextlabel_button_Callback(hObject, eventdata, handles)
% hObject    handle to nextlabel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = handles.trx(handles.fly).off+(handles.f);
j = find(handles.segstarts > i,1);
if isempty(j), return; end
j = handles.segstarts(j);
handles.f = handles.trx(handles.fly).firstframe + j - 1;
SetFrameNumber(handles,hObject);
handles = UpdatePlot(handles);
guidata(hObject,handles);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'minfigsz'),
  return;
end
figpos = get(handles.figure1,'Position');
figpos0 = figpos;
figpos(3) = max(figpos(3),handles.minfigsz(1));
figpos(4) = max(figpos(4),handles.minfigsz(2));
if any(figpos ~= figpos0),
  set(handles.figure1,'Position',figpos);
end

axespos = get(handles.mainaxes,'Position');
axespos(3) = figpos(3) - axespos(1) - handles.axesoffright;
axespos(4) = figpos(4) - axespos(2) - handles.axesofftop;
set(handles.mainaxes,'Position',axespos);

panelpos = get(handles.uipanel_bottom,'Position');
c = axespos(1) + axespos(3)/2;
panelpos(1) = c - panelpos(3)/2;
set(handles.uipanel_bottom,'Position',panelpos);

panelpos = get(handles.uipanel_right,'Position');
panelpos(1) = figpos(3) - handles.paneloffright;
set(handles.uipanel_right,'Position',panelpos);

%keyboard;



function edit_pathlength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathlength as text
%        str2double(get(hObject,'String')) returns contents of edit_pathlength as a double
pathlength = str2double(get(hObject,'String'));
if isnan(pathlength),
  set(hObject,'String',num2str(handles.PATHLENGTH));
else
  handles.PATHLENGTH = pathlength;
  handles = UpdatePlot(handles);
  guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function edit_pathlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
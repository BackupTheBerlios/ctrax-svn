function varargout = labelbehaviorssocial(varargin)
% [starts,ends] = labelbehaviorssocial(trx,fly,moviename)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @labelbehaviorssocial_OpeningFcn, ...
                   'gui_OutputFcn',  @labelbehaviorssocial_OutputFcn, ...
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


% --- Executes just before labelbehaviorssocial is made visible.
function labelbehaviorssocial_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to labelbehaviorssocial (see VARARGIN)

% Choose default command line output for labelbehaviorssocial
if length(varargin) < 3 || length(varargin) == 4 || length(varargin) > 8,
  error('Usage: [starts,ends] = labelbehaviorssocial(maintrk,trx,fly,moviename, <starts,ends,otherfly>,<labels>,<scattercommand>)');
end
handles.output = hObject;
handles.trx = varargin{1};
handles.fly = varargin{2};
handles.moviename = varargin{3};
if length(varargin) >= 6,
  handles.segstarts = varargin{4};
  handles.segends = varargin{5};
  handles.otherfly = varargin{6};
else
  handles.segstarts = [];
  handles.segends = [];
  handles.otherfly = [];
end
if length(varargin) >= 7 && ~isempty(varargin{7}),
  handles.labels = varargin{7};
else
  handles.labels = struct('f',cell(1,0),'s',cell(1,0));
end
if length(varargin) >= 8 && ~isempty(varargin{8}),
  handles.scattercommand = varargin{8};
else
  handles.scattercommand = 'min(trk.velmag,7)';
end
set(handles.scatteredit,'string',handles.scattercommand);
if length(varargin) >= 9 && ~isempty(varargin{9}),
  handles.segcolors = varargin{9};
end

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
handles.hotherflyline = [];
% setup gui
InitializeFrameSlider(handles);
SetFrameNumber(handles);
InitializeDisplayPanel(handles);

% plot first frame
handles = PlotFirstFrame(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes labelbehaviorssocial wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = labelbehaviorssocial_OutputFcn(hObject, eventdata, handles) 
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
varargout{3} = handles.otherfly;
if ~isfield(handles,'labels'),
  varargout{4} = [];
else
  varargout{4} = handles.labels;
end
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
colors = lines(6);
for i = 1:length(handles.segstarts),
  i1 = handles.segstarts(i); i2 = handles.segends(i);
  f1 = handles.segstarts(i) + handles.trx(handles.fly).firstframe - 1;
  f2 = handles.segends(i) + handles.trx(handles.fly).firstframe - 1;
  handles.isseg(i1:i2) = ...
    plot(handles.trx(handles.fly).x(i1:i2),...
    handles.trx(handles.fly).y(i1:i2),'m','linewidth',4,'hittest','off');
  handles.hstarts(i) = plot(handles.trx(handles.fly).x(i1),handles.trx(handles.fly).y(i1),'o','color',colors(mod(i-1,6)+1,:));
  handles.hends(i) = plot(handles.trx(handles.fly).x(i2),handles.trx(handles.fly).y(i2),'s','color',colors(mod(i-1,6)+1,:));
  fmid = ceil((f1+f2)/2);
  imid1 = handles.trx(handles.fly).f2i(fmid);
  otherfly = handles.otherfly(i);
  imid2 = handles.trx(otherfly).f2i(fmid);
  handles.hotherflyline(i) = ...
    plot([handles.trx(handles.fly).x(imid1),handles.trx(otherfly).x(imid2)],...
    [handles.trx(handles.fly).y(imid1),handles.trx(otherfly).y(imid2)],...
    '-','color',colors(mod(i-1,6)+1,:),'hittest','off');
  if isfield(handles,'segcolors'),
    set(handles.isseg(i1),'color',handles.segcolors(i,:));
  end
end
plot(handles.trx(handles.fly).x,handles.trx(handles.fly).y,'-','color',[.7,.7,.7],'hittest','off');
handles.hpath = myscatter(handles.trx(handles.fly).x,handles.trx(handles.fly).y,[],handles.c,[],@jet,'.');
caxis([min(handles.c),max(handles.c)]);
colormap jet;
%colorbar;
set(handles.hpath,'hittest','off');
handles.hmarker = drawflyo(0,0,0,10,10);
handles.nflies = length(handles.trx);
handles.hotherflymarker = nan(1,handles.nflies);
colors = jet(handles.nflies);
[tmp,replacecolor] = min( (1-colors(:,1)) + colors(:,2) + colors(:,3) );
colors = [colors(1:replacecolor-1,:);colors(replacecolor+1:end,:)];
colors = [colors(1:handles.fly-1,:);[1,0,0];colors(handles.fly:end,:)];
handles.otherflycolors = colors;
if handles.fly ~= 1,
  handles.otherflyselected = 1;
else
  handles.otherflyselected = 2;
end
for i = 1:handles.nflies,
  handles.hotherflymarker(i) = drawflyo(0,0,0,10,10);
  set(handles.hotherflymarker(i),'color',colors(i,:),'linewidth',1);
  set(handles.hotherflymarker(i),'buttondownfcn',...
    @(hObject,eventdata,handles)labelbehaviorssocial(...
    'otherfly_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
end
i = handles.otherflyselected;
set(handles.hotherflymarker(i),'color',(1+colors(i,:))/2,'linewidth',4);
handles.hcenter = plot(0,0,'o','color','r','hittest','off');
set(handles.hmarker,'color','r','linewidth',2,'hittest','off');
handles.hotherflyselectedline = plot([0,0],[0,0],'w-','hittest','off');

handles.hlabel = zeros(size(handles.labels));
for i = 1:length(handles.labels),
  j = handles.trx(handles.fly).f2i(handles.labels(i).f);
  handles.hlabel(i) = text(handles.trx(handles.fly).x(j),...
    handles.trx(handles.fly).y(j),handles.labels(i).s,'horizontalalignment','left',...
    'color','r','fontname','times','fontsize',16,'clipping','on');
  if handles.isseg(handles.labels(i).f),
    set(handles.hlabel(i),'color',get(handles.isseg(handles.labels(i).f),'color'));
  end
end
axes(handles.zoomlocaxes);
image(uint8(repmat(im,[1,1,3]))); axis image; axis off;
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

fly = handles.fly;
i = handles.trx(fly).f2i(handles.f);
updatefly(handles.hmarker,handles.trx(fly).x(i),handles.trx(fly).y(i),...
  handles.trx(fly).theta(i),handles.trx(fly).a(i),handles.trx(fly).b(i));
for j = 1:handles.nflies,
  if j == fly,
    continue;
  end
  if handles.f < handles.trx(j).firstframe || handles.f > handles.trx(j).endframe,
    set(handles.hotherflymarker(j),'visible','off');
    continue;
  else
    set(handles.hotherflymarker(j),'visible','on');
  end
  updatefly(handles.hotherflymarker(j),handles.trx(j),handles.trx(j).f2i(handles.f));
end
set(handles.hcenter,'xdata',handles.trx(fly).x(i),'ydata',handles.trx(fly).y(i));

% for j = 1:length(handles.segstarts),
%   if handles.f < handles.segstarts(j) || handles.f > handles.segends(j),
%     set(handles.hotherflyline(j),'visible','off');
%   else
%     fly2 = handles.otherfly(j);
%     k = handles.trx(fly2).f2i(handles.f);
%     set(handles.hotherflyline(j),'visible','on',...
%       'xdata',[handles.trx(fly).x(i),handles.trx(fly2).x(k)],...
%       'ydata',[handles.trx(fly).y(i),handles.trx(fly2).y(k)]);
%   end
% end

otherfly = handles.otherflyselected;
i2 = handles.trx(otherfly).f2i(handles.f);
set(handles.hotherflyselectedline,...
  'xdata',[handles.trx(handles.fly).x(i),handles.trx(otherfly).x(i2)],...
  'ydata',[handles.trx(handles.fly).y(i),handles.trx(otherfly).y(i2)]);

if get(handles.zoombutton,'value'),
  ZoomInOnFlies(handles);
end

set(handles.mainaxes,'buttondownfcn',handles.bd);

function ZoomInOnFlies(handles)

x0 = inf;
x1 = 0;
y0 = inf;
y1 = 0;
fly = handles.fly;
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
  i = handles.trx(handles.fly).f2i(handles.f);
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

function otherfly_ButtonDownFcn(hObject,eventdata,handles)

if strcmpi(get(handles.hzoom,'enable'),'on'), return; end
fly = find(handles.hotherflymarker == hObject,1);
oldfly = handles.otherflyselected;
set(handles.hotherflymarker(oldfly),'color',handles.otherflycolors(oldfly,:),...
  'linewidth',1);
set(handles.hotherflymarker(fly),'color',(1+handles.otherflycolors(fly,:))/2,...
  'linewidth',4);
handles.otherflyselected = fly;
i = handles.trx(handles.fly).f2i(handles.f);
i2 = handles.trx(fly).f2i(handles.f);
set(handles.hotherflyselectedline,...
  'xdata',[handles.trx(handles.fly).x(i),handles.trx(fly).x(i2)],...
  'ydata',[handles.trx(handles.fly).y(i),handles.trx(fly).y(i2)]);

guidata(handles.figure1,handles);

% --- Executes on mouse press over axes background.
function mainaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(handles.hzoom,'enable'),'on'), return; end
pt = get(handles.mainaxes,'currentpoint');
x = pt(1,1); y = pt(1,2);
d = sqrt((handles.trx(handles.fly).x-x).^2 + ...
  (handles.trx(handles.fly).y-y).^2);
[mind,t] = min(d);
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
i1 = handles.trx(handles.fly).f2i(f1);
if handles.isseg(i1)>0 || f1 > handles.trx(handles.fly).endframe || ...
    f1 < handles.trx(handles.fly).firstframe,
  return;
end
i2 = find(handles.isseg(i1+1:end)>0,1,'first');
if isempty(i2),
  i2 = handles.trx(handles.fly).nframes;
else
  i2 = i2 + i1 - 1;
end
handles.currentstarti = i1;
%i1 = handles.trx(handles.fly).f2i(f1);
%i2 = handles.trx(handles.fly).f2i(f2);
colors = lines(6);
i = mod(length(handles.segstarts)-1,6)+1;
handles.isseg(i1:i2) = plot(handles.trx(handles.fly).x(i1:i2),...
  handles.trx(handles.fly).y(i1:i2),'-','color',colors(i,:),'linewidth',4,'hittest','off');
chil = get(handles.mainaxes,'children');
chil = [chil(2:end-1);chil(1);chil(end)];
set(handles.mainaxes,'children',chil);
handles.segstarts(end+1) = i1;
handles.segends(end+1) = i2;
fly2 = handles.otherflyselected;
handles.otherfly(end+1) = fly2;
handles.hstarts(end+1) = plot(handles.trx(handles.fly).x(i1),handles.trx(handles.fly).y(i1),'o','color',colors(i,:));
handles.hends(end+1) = plot(handles.trx(handles.fly).x(i2),handles.trx(handles.fly).y(i2),'s','color',colors(i,:));

j = handles.trx(fly2).f2i(f1);
x = [handles.trx(handles.fly).x(i1),handles.trx(fly2).x(j)];
y = [handles.trx(handles.fly).y(i1),handles.trx(fly2).y(j)];
handles.hotherflyline(end+1) = plot(x,y,'.-','color',colors(i,:),'hittest','off');

guidata(hObject,handles);

% --- Executes on button press in addendbutton.
function addendbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addendbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% end last interval at f1
f1 = handles.f;
i1 = handles.trx(handles.fly).f2i(f1);
if handles.isseg(i1)<=0 || f1 > handles.trx(handles.fly).endframe || ...
    f1 < handles.trx(handles.fly).firstframe,
  return;
end

% reset until f2 to -1
i2 = handles.segends(end);
handles.isseg(i1+1:i2) = -1;

% start of this interval
i0 = handles.segstarts(end);
%i0 = find(handles.isseg(1:i1-1)<=0,1,'last');
%if isempty(i0),
%  i0 = 1;
%  %f0 = handles.trx(handles.fly).firstframe;
%else
%  i0 = i0 + 1;
%  %f0 = handles.trx(handles.fly).firstframe + i0 - 1;
%end

% set plot data
set(handles.isseg(i1),'xdata',handles.trx(handles.fly).x(i0:i1),...
  'ydata',handles.trx(handles.fly).y(i0:i1));
set(handles.hends(end),'xdata',handles.trx(handles.fly).x(i1),...
  'ydata',handles.trx(handles.fly).y(i1));
handles.segends(end) = i1;

guidata(hObject,handles);

starts = handles.segstarts + handles.trx(handles.fly).firstframe-1;
ends = handles.segends + handles.trx(handles.fly).firstframe-1;
if ~isfield(handles,'labels'),
  labels = [];
else
  labels = handles.labels;
end

save labeledbehaviors_autosave.mat starts ends labels

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
i = handles.trx(handles.fly).f2i(f);
if handles.isseg(i) <= 0,
  fprintf('No label at this frame. Not deleting\n');
  return;
end
fidx = handles.trx(handles.fly).f2i(f);
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
delete(handles.hotherflyline(j));
handles.segstarts(j) = [];
handles.segends(j) = [];
handles.otherfly(j) = [];
handles.hstarts(j) = [];
handles.hends(j) = [];
handles.hotherflyline(j) = [];
handles.isseg(i0:i1) = -1;
if isfield(handles,'labels') && ~isempty(handles.labels),
  fs = getstructarrayfield(handles.labels,'f');
  f0 = i0 + off; f1 = i1 + off;
  dodelete = fs >= f0 & fs <= f1;
  if any(dodelete),
    delete(handles.hlabel(dodelete));
    handles.labels(dodelete) = [];
  end
end

guidata(hObject,handles);


% --- Executes on button press in addlabelbutton.
function addlabelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addlabelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f = handles.f;
if handles.isseg(f) <= 0,
  return;
end

s = inputdlg('Enter label for segment: ');
newlabel.f = f;
newlabel.s = s;
i = handles.trx(handles.fly).f2i(f);
newh = text(handles.trx(handles.fly).x(i),...
  handles.trx(handles.fly).y(i),s,'horizontalalignment','left',...
  'color','r','fontname','times','fontsize',16,'clipping','on');
if ~isfield(handles,'labels'),
  handles.labels = newlabel;
  handles.hlabel = newh;
else
  handles.labels(end+1) = newlabel;
  handles.hlabel(end+1) = newh;
end
guidata(hObject,handles);


% --- Executes on button press in prevlabel_button.
function prevlabel_button_Callback(hObject, eventdata, handles)
% hObject    handle to prevlabel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = handles.trx(handles.fly).f2i(handles.f);
j = find(handles.isseg(1:i-2)<=0&handles.isseg(2:i-1)>0,1,'last');
if isempty(j),
  return;
end
j = j + 1;
handles.f = handles.trx(handles.fly).firstframe + j - 1;
SetFrameNumber(handles,hObject);
handles = UpdatePlot(handles);
guidata(hObject,handles);

% --- Executes on button press in nextlabel_button.
function nextlabel_button_Callback(hObject, eventdata, handles)
% hObject    handle to nextlabel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = handles.trx(handles.fly).f2i(handles.f);
j = find(handles.isseg(i:end-1)<=0&handles.isseg(i+1:end)>0,1,'first');
if isempty(j),
  return;
end
j = j + i;
handles.f = handles.trx(handles.fly).firstframe + j - 1;
SetFrameNumber(handles,hObject);
handles = UpdatePlot(handles);
guidata(hObject,handles);


% --- Executes on button press in otherflybutton.
function otherflybutton_Callback(hObject, eventdata, handles)
% hObject    handle to otherflybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



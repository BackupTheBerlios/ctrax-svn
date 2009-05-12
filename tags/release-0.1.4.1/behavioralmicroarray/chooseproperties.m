function varargout = chooseproperties(varargin)
% CHOOSEPROPERTIES M-file for chooseproperties.fig
%      CHOOSEPROPERTIES, by itself, creates a new CHOOSEPROPERTIES or raises the existing
%      singleton*.
%
%      H = CHOOSEPROPERTIES returns the handle to a new CHOOSEPROPERTIES or the handle to
%      the existing singleton*.
%
%      CHOOSEPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOOSEPROPERTIES.M with the given input arguments.
%
%      CHOOSEPROPERTIES('Property','Value',...) creates a new CHOOSEPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chooseproperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chooseproperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chooseproperties

% Last Modified by GUIDE v2.5 04-Dec-2008 09:34:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chooseproperties_OpeningFcn, ...
                   'gui_OutputFcn',  @chooseproperties_OutputFcn, ...
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


% --- Executes just before chooseproperties is made visible.
function chooseproperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chooseproperties (see VARARGIN)

% ARGUMENTS:
% PROPERTIES: cell of strings
% [UNITS]: struct
% [PARAMS]: parameters output in previous running

% Choose default command line output for chooseproperties
handles.output = hObject;

% initialize properties
arg = 2;
if length(varargin) >= 2 && isstruct(varargin{2}),
  fnunits = fieldnames(varargin{2});
  if ~isempty(fnunits) && ...
      isfield(varargin{2}.(fnunits{1}),'den') && ...
      isfield(varargin{2}.(fnunits{1}),'num'),
    handles.units = varargin{2};
    arg = 3;
  end
end

% set constants
handles = setconstants(handles);
  
% input property names
handles.properties = setdiff(varargin{1},handles.ignore);
handles.nproperties = length(handles.properties);

if length(varargin) >= arg,
  handles = initializeproperties(handles,varargin{arg});
else
  handles = initializeproperties(handles);
end

% initialize window
handles = initializewindow(handles);

% update window
updatewindow(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chooseproperties wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = chooseproperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

params.minxfns = handles.properties(handles.minproperties(:,1));
params.maxxfns = handles.properties(handles.maxproperties(:,1));
params.minxclosefns = handles.properties(handles.minproperties(:,2));
params.maxxclosefns = handles.properties(handles.maxproperties(:,2));
params.minsumxfns = handles.properties(handles.minproperties(:,3));
params.maxsumxfns = handles.properties(handles.maxproperties(:,3));
params.minmeanxfns = handles.properties(handles.minproperties(:,4));
params.maxmeanxfns = handles.properties(handles.maxproperties(:,4));
params.minr = handles.minr;
params.maxr = handles.maxr;
params.options = {};

% find forced minimums
idx = handles.minproperties & (handles.minlower == handles.minupper);
i = 1;
if any(idx(:,i)),
  forceminx = {};
  for j = find(idx(:,i))',
    forceminx{end+1} = handles.properties{j};
    forceminx{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'forceminx';
  params.options{end+1} = forceminx;
end
i = 2;
if any(idx(:,i)),
  forceminxclose = {};
  for j = find(idx(:,i))',
    forceminxclose{end+1} = handles.properties{j};
    forceminxclose{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'forceminxclose';
  params.options{end+1} = forceminxclose;
end
i = 3;
if any(idx(:,i)),
  forceminsumx = {};
  for j = find(idx(:,i))',
    forceminsumx{end+1} = handles.properties{j};
    forceminsumx{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'forceminsumx';
  params.options{end+1} = forceminsumx;
end
i = 4;
if any(idx(:,i)),
  forceminmeanx = {};
  for j = find(idx(:,i))',
    forceminmeanx{end+1} = handles.properties{j};
    forceminmeanx{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'forceminmeanx';
  params.options{end+1} = forceminmeanx;
end
% remove
if any(idx(:)),
  handles.minlower(idx) = -inf;
  handles.minupper(idx) = inf;
end

% find forced maximums
idx = handles.maxproperties & (handles.maxlower == handles.maxupper);
i = 1;
if any(idx(:,i)),
  forcemaxx = {};
  for j = find(idx(:,i))',
    forcemaxx{end+1} = handles.properties{j};
    forcemaxx{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'forcemaxx';
  params.options{end+1} = forcemaxx;
end
i = 2;
if any(idx(:,i)),
  forcemaxxclose = {};
  for j = find(idx(:,i))',
    forcemaxxclose{end+1} = handles.properties{j};
    forcemaxxclose{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'forcemaxxclose';
  params.options{end+1} = forcemaxxclose;
end
i = 3;
if any(idx(:,i)),
  forcemaxsumx = {};
  for j = find(idx(:,i))',
    forcemaxsumx{end+1} = handles.properties{j};
    forcemaxsumx{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'forcemaxsumx';
  params.options{end+1} = forcemaxsumx;
end
i = 4;
if any(idx(:,i)),
  forcemaxmeanx = {};
  for j = find(idx(:,i))',
    forcemaxmeanx{end+1} = handles.properties{j};
    forcemaxmeanx{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'forcemaxmeanx';
  params.options{end+1} = forcemaxmeanx;
end
% remove
if any(idx(:)),
  handles.maxlower(idx) = -inf;
  handles.maxupper(idx) = inf;
end

% find lower bounds placed on minimums
idx = handles.minproperties & handles.minlower > -inf;
i = 1;
if any(idx(:,i)),
  minminx = {};
  for j = find(idx(:,i))',
    minminx{end+1} = handles.properties{j};
    minminx{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'minminx';
  params.options{end+1} = minminx;
end
i = 2;
if any(idx(:,i)),
  minminxclose = {};
  for j = find(idx(:,i))',
    minminxclose{end+1} = handles.properties{j};
    minminxclose{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'minminxclose';
  params.options{end+1} = minminxclose;
end
i = 3;
if any(idx(:,i)),
  minminsumx = {};
  for j = find(idx(:,i))',
    minminsumx{end+1} = handles.properties{j};
    minminsumx{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'minminsumx';
  params.options{end+1} = minminsumx;
end
i = 4;
if any(idx(:,i)),
  minminmeanx = {};
  for j = find(idx(:,i))',
    minminmeanx{end+1} = handles.properties{j};
    minminmeanx{end+1} = handles.minlower(j,i);
  end
  params.options{end+1} = 'minminmeanx';
  params.options{end+1} = minminmeanx;
end

% find upper bounds placed on minimums
idx = handles.minproperties & handles.minupper < inf;
i = 1;
if any(idx(:,i)),
  maxminx = {};
  for j = find(idx(:,i))',
    maxminx{end+1} = handles.properties{j};
    maxminx{end+1} = handles.minupper(j,i);
  end
  params.options{end+1} = 'maxminx';
  params.options{end+1} = maxminx;
end
i = 2;
if any(idx(:,i)),
  maxminxclose = {};
  for j = find(idx(:,i))',
    maxminxclose{end+1} = handles.properties{j};
    maxminxclose{end+1} = handles.minupper(j,i);
  end
  params.options{end+1} = 'maxminxclose';
  params.options{end+1} = maxminxclose;
end
i = 3;
if any(idx(:,i)),
  maxminsumx = {};
  for j = find(idx(:,i))',
    maxminsumx{end+1} = handles.properties{j};
    maxminsumx{end+1} = handles.minupper(j,i);
  end
  params.options{end+1} = 'maxminsumx';
  params.options{end+1} = maxminsumx;
end
i = 4;
if any(idx(:,i)),
  maxminmeanx = {};
  for j = find(idx(:,i))',
    maxminmeanx{end+1} = handles.properties{j};
    maxminmeanx{end+1} = handles.minupper(j,i);
  end
  params.options{end+1} = 'maxminmeanx';
  params.options{end+1} = maxminmeanx;
end

% find lower bounds placed on maximums
idx = handles.maxproperties & handles.maxlower > -inf;
i = 1;
if any(idx(:,i)),
  minmaxx = {};
  for j = find(idx(:,i))',
    minmaxx{end+1} = handles.properties{j};
    minmaxx{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'minmaxx';
  params.options{end+1} = minmaxx;
end
i = 2;
if any(idx(:,i)),
  minmaxxclose = {};
  for j = find(idx(:,i))',
    minmaxxclose{end+1} = handles.properties{j};
    minmaxxclose{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'minmaxxclose';
  params.options{end+1} = minmaxxclose;
end
i = 3;
if any(idx(:,i)),
  minmaxsumx = {};
  for j = find(idx(:,i))',
    minmaxsumx{end+1} = handles.properties{j};
    minmaxsumx{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'minmaxsumx';
  params.options{end+1} = minmaxsumx;
end
i = 4;
if any(idx(:,i)),
  minmaxmeanx = {};
  for j = find(idx(:,i))',
    minmaxmeanx{end+1} = handles.properties{j};
    minmaxmeanx{end+1} = handles.maxlower(j,i);
  end
  params.options{end+1} = 'minmaxmeanx';
  params.options{end+1} = minmaxmeanx;
end

% find upper bounds placed on maximums
idx = handles.maxproperties & handles.maxupper < inf;
i = 1;
if any(idx(:,i)),
  maxmaxx = {};
  for j = find(idx(:,i))',
    maxmaxx{end+1} = handles.properties{j};
    maxmaxx{end+1} = handles.maxupper(j,i);
  end
  params.options{end+1} = 'maxmaxx';
  params.options{end+1} = maxmaxx;
end
i = 2;
if any(idx(:,i)),
  maxmaxxclose = {};
  for j = find(idx(:,i))',
    maxmaxxclose{end+1} = handles.properties{j};
    maxmaxxclose{end+1} = handles.maxupper(j,i);
  end
  params.options{end+1} = 'maxmaxxclose';
  params.options{end+1} = maxmaxxclose;
end
i = 3;
if any(idx(:,i)),
  maxmaxsumx = {};
  for j = find(idx(:,i))',
    maxmaxsumx{end+1} = handles.properties{j};
    maxmaxsumx{end+1} = handles.maxupper(j,i);
  end
  params.options{end+1} = 'maxmaxsumx';
  params.options{end+1} = maxmaxsumx;
end
i = 4;
if any(idx(:,i)),
  maxmaxmeanx = {};
  for j = find(idx(:,i))',
    maxmaxmeanx{end+1} = handles.properties{j};
    maxmaxmeanx{end+1} = handles.maxupper(j,i);
  end
  params.options{end+1} = 'maxmaxmeanx';
  params.options{end+1} = maxmaxmeanx;
end

varargout{1} = params;
delete(handles.figure1);

function lowerbnd_Callback(hObject, eventdata, handles)
% hObject    handle to lowerbnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerbnd as text
%        str2double(get(hObject,'String')) returns contents of lowerbnd as
%        a double
v = str2double(get(hObject,'String'));
i = handles.editing;
j = handles.lastproperty(i);
ismax = handles.lastpropertytype(i);
if ismax,
  prevv = handles.maxlower(j,i);
  upperv = handles.maxupper(j,i);
else
  prevv = handles.minlower(j,i);
  upperv = handles.minupper(j,i);
end
if v == prevv,
  return;
end
if isnan(v) || v > upperv,
  set(hObject,'string',num2str(prevv));
  return;
end
if ismax,
  handles.maxlower(j,i) = v;
else
  handles.minlower(j,i) = v;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lowerbnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerbnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperbnd_Callback(hObject, eventdata, handles)
% hObject    handle to upperbnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperbnd as text
%        str2double(get(hObject,'String')) returns contents of upperbnd as a double

v = str2double(get(hObject,'String'));
i = handles.editing;
j = handles.lastproperty(i);
ismax = handles.lastpropertytype(i);
if ismax,
  prevv = handles.maxupper(j,i);
  lowerv = handles.maxlower(j,i);
else
  prevv = handles.minupper(j,i);
  lowerv = handles.minlower(j,i);
end
if v == prevv,
  return;
end
if isnan(v) || v < lowerv,
  set(hObject,'string',num2str(prevv));
  return;
end
if ismax,
  handles.maxupper(j,i) = v;
else
  handles.minupper(j,i) = v;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function upperbnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperbnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edit2.
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_Callback(hObject,eventdata,handles,2);

% --- Executes on button press in edit1.
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit_Callback(hObject,eventdata,handles,1);

% --- Executes on button press in edit3.
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit_Callback(hObject,eventdata,handles,3);

function edit_Callback(hObject,eventdata,handles,i)

newediting = handles.notediting(i);
handles.notediting(i) = handles.editing;
handles.editing = newediting;
guidata(hObject,handles);
updatenoteditable(handles,i);
updateeditable(handles);

% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);

function handles = setconstants(handles)

% condition definitions

handles.conditiontitle{1} = '1. Per-Frame Bounds';
handles.conditiontitle{2} = '2. Near-Frame Bounds';
handles.conditiontitle{3} = '3. Sequence Sum Bounds';
handles.conditiontitle{4} = '4. Sequence Mean Bounds';

handles.maintext{1} = {'In each frame t1,...,t2, properties ',...
  'of the fly are within given ranges.'};
handles.maintext{2} = {'In each frame t1,...,t2, properties ',...
  'are temporally near (within a given',...
  'number of frames) frames in which',...
  'the properties are within tighter',...
  'ranges'};
handles.maintext{3} = {'The summed properties of the fly''s ',...
  'trajectory in t1:t2 are within ',...
  'given ranges.'};
handles.maintext{4} = {'The mean properties of the fly''s ',...
  'trajectory in t1:t2 are within ',...
  'given ranges.'};

handles.propdescrs.du_ctr = ...
  {'Forward velocity of center'};
handles.propdescrs.du_cor = ...
  {'Forward velocity of center',...
  'of rotation'};
handles.propdescrs.corfrac = ...
  {'Location of center of ',...
  'rotation, in half body lengths'};
handles.propdescrs.velmag_ctr = ...
  {'Speed of center'};
handles.propdescrs.velmag = ...
  {'Speed of center of rotation'};
handles.propdescrs.accmag = ...
  {'Acceleration magnitude of ',...
  'center'};
handles.propdescrs.absdv_cor = ...
  {'Sideways speed of center,',...
  'of rotation'};
handles.propdescrs.flipdv_cor = ...
  {'Sideways speed of center,',...
  'of rotation',...
  'flipped if dtheta < 0'};
handles.propdescrs.absdtheta = ...
  {'Absol change in orientation'};
handles.propdescrs.absd2theta = ...
  {'Absol angular acceleration'};
handles.propdescrs.abssmoothdtheta = ...
  {'Smoothed angular speed'};
handles.propdescrs.abssmoothd2theta = ...
  {'Smoothed absol angular acc'};
handles.propdescrs.absyaw = ...
  {'Absol difference btn heading',...
  'and orientation'};
handles.propdescrs.du_tail = ...
  {'Forward vel. of tail'};
handles.propdescrs.absdu_tail = ...
  {'Forward speed of tail'};
handles.propdescrs.absdv_tail = ...
  {'Sideways speed of tail'};
handles.propdescrs.absdtheta_tail = ...
  {'Absolute rotation of nose ',...
  'around tail'};
handles.propdescrs.absphisideways = ...
  {'Abs angle between sideways ',...
  'and heading'};
handles.propdescrs.dcenter = ...
  {'Distance between flys'' centroids'};
handles.propdescrs.dnose2ell = ...
  {'Distance from first fly''s nose ',...
  'to second fly''s ellipse'};
handles.propdescrs.distnose2ell = ...
  {'Distance from first fly''s nose ',...
  'to second fly''s ellipse'};
handles.propdescrs.dell2nose = ...
  {'Distance from second fly''s nose ',...
  'to first fly''s ellipse'};
handles.propdescrs.magveldiff = ...
  {'Magnitude of difference in ',...
  'velocity of two flies'};
handles.propdescrs.veltoward = ...
  {'Velocity of first fly in ',...
  'direction of second fly'};
handles.propdescrs.absthetadiff = ...
  {'Absolute difference in ',...
  'orientation of two flies'};
handles.propdescrs.minvelmag = ...
  {'Minimum speed of either of',...
  'two flies'};
handles.propdescrs.maxvelmag = ...
  {'Maximum speed of either of',...
  'two flies'};
handles.propdescrs.anglesub = ...
  {'Fly 1''s view',...
  'covered by fly 2'};
handles.propdescrs.absphidiff = ...
  {'Absolute difference in ',...
  'heading of two flies'};
handles.propdescrs.absanglefrom1to2 = ...
  {'Absolute direction to fly 2',...
  'from fly 1'};
if isfield(handles,'units'),
  fns = fieldnames(handles.propdescrs);
  for i = 1:length(fns),
    if isfield(handles.units,fns{i}),
      s = unitsstring(handles,fns{i});
      handles.propdescrs.(fns{i}){end+1} = sprintf('[ %s ]',s);
    end
  end
end
  
handles.ignore = {'x','y','theta','a','b','id','moviename','firstframe','arena','f2i',...
  'nframes','endframe','dx','dy','v','matname','annname','dtheta','dv_ctr','x_cor_curr',...
  'y_cor_curr','x_cor_next','y_cor_next','dv_cor','signdtheta','smoothdtheta','smoothd2theta',...
  'phi','yaw','dv_tail','dtheta_tail','phisideways','d2theta','smooththeta',...
  'thetadiff','phidiff','anglefrom1to2','fps','pxpermm','units'};

handles.nconditions = 4;

function handles = initializeproperties(handles,params)

% we are currently editing the first condition
handles.editing = 1;
handles.notediting = 2:handles.nconditions;

% currently, we have no properties and no bounds
handles.minproperties = false(handles.nproperties,handles.nconditions);
handles.maxproperties = false(handles.nproperties,handles.nconditions);

% and no bounds
handles.minlower = -inf(handles.nproperties,handles.nconditions);
handles.maxlower = -inf(handles.nproperties,handles.nconditions);
handles.minupper = inf(handles.nproperties,handles.nconditions);
handles.maxupper = inf(handles.nproperties,handles.nconditions);
handles.minr = 1;
handles.maxr = inf;

% last property edited
handles.lastproperty = ones(1,handles.nconditions);
handles.lastpropertytype = false(1,handles.nconditions);

if exist('params','var'),
  if isfield(params,'minr'),
    handles.minr = params.minr;
  end
  if isfield(params,'maxr'),
    handles.maxr = params.maxr;
  end
  
  i = 1;
  for j = 1:length(params.minxfns),
    k = find(strcmpi(handles.properties,params.minxfns{j}),1);
    handles.minproperties(k,i) = true;
  end
  i = 2;
  for j = 1:length(params.minxclosefns),
    k = find(strcmpi(handles.properties,params.minxclosefns{j}),1);
    handles.minproperties(k,i) = true;
  end
  i = 3;
  for j = 1:length(params.minsumxfns),
    k = find(strcmpi(handles.properties,params.minsumxfns{j}),1);
    handles.minproperties(k,i) = true;
  end
  i = 4;
  for j = 1:length(params.minmeanxfns),
    k = find(strcmpi(handles.properties,params.minmeanxfns{j}),1);
    handles.minproperties(k,i) = true;
  end
  i = 1;
  for j = 1:length(params.maxxfns),
    k = find(strcmpi(handles.properties,params.maxxfns{j}),1);
    handles.maxproperties(k,i) = true;
  end
  i = 2;
  for j = 1:length(params.maxxclosefns),
    k = find(strcmpi(handles.properties,params.maxxclosefns{j}),1);
    handles.maxproperties(k,i) = true;
  end
  i = 3;
  for j = 1:length(params.maxsumxfns),
    k = find(strcmpi(handles.properties,params.maxsumxfns{j}),1);
    handles.maxproperties(k,i) = true;
  end
  i = 4;
  for j = 1:length(params.maxmeanxfns),
    k = find(strcmpi(handles.properties,params.maxmeanxfns{j}),1);
    handles.maxproperties(k,i) = true;
  end

  for j = 1:2:length(params.options),
    n = params.options{j};
    v = params.options{j+1};
    
    try

      % first will say what it is doing to the bound
      if strcmpi(n(1:3),'min'),
        islower = true;
        isupper = false;
        n = n(4:end);
      elseif strcmpi(n(1:3),'max')
        islower = false;
        isupper = true;
        n = n(4:end);
      elseif strcmpi(n(1:5),'force')
        islower = true;
        isupper = true;
        n = n(6:end);
      else
        continue;
      end
      
      % then will say whether min/max bound
      if strcmpi(n(1:3),'min'),
        ismax = false;
      elseif strcmpi(n(1:3),'max'),
        ismax = true;
      else
        continue;
      end
      n = n(4:end);
      
      % then will say which condition
      if strcmpi(n,'x'),
        i = 1;
      elseif strcmpi(n,'xclose'),
        i = 2;
      elseif strcmpi(n,'sumx'),
        i = 3;
      elseif strcmpi(n,'meanx'),
        i = 4;
      else
        continue;
      end
      
      % find which properties are being bounded
      for k = 1:2:length(v),
        try
          l = strcmpi(v{k},handles.properties);
          if ~any(l) || ~isnumeric(v{k+1}),
            continue;
          end
          if islower,
            if ismax,
              handles.maxlower(l,i) = v{k+1};
            else
              handles.minlower(l,i) = v{k+1};
            end
          end
          if isupper,
            if ismax,
              handles.maxupper(l,i) = v{k+1};
            else
              handles.minupper(l,i) = v{k+1};
            end
          end
        catch
          continue;
        end
      end

    catch
      continue;
    end
  end
end

function handles = initializewindow(handles)

panelpos = get(handles.minproppanel,'position');
panelht = panelpos(4);
checkboxpos = get(handles.mincheckbox1,'position');
checkboxht = checkboxpos(4);
handles.ncheckshow = min(floor(panelht/checkboxht),handles.nproperties);
un = get(handles.mincheckbox1,'units');

% create checkboxes
delete(handles.mincheckbox1);
handles.mincheckbox = zeros(1,handles.ncheckshow);
handles.maxcheckbox = zeros(1,handles.ncheckshow);
handles.checkboxprops = 1:handles.ncheckshow;
for i = 1:handles.ncheckshow,
  st = handles.properties{i};
  y = panelht - checkboxht*i;
  cbp = checkboxpos;
  cbp(2) = y;
  handles.mincheckbox(i) = uicontrol('style','checkbox',...
    'parent',handles.minproppanel,...
    'string',st,'value',0,'units',un,...
    'Callback',@(hObject,eventdata)chooseproperties('mincheckbox_Callback',hObject,eventdata,guidata(hObject)));
  set(handles.mincheckbox(i),'position',cbp);
  handles.maxcheckbox(i) = uicontrol('style','checkbox',...
    'parent',handles.maxproppanel,...
    'string',st,'value',0,'units',un,...
    'Callback',@(hObject,eventdata)chooseproperties('maxcheckbox_Callback',hObject,eventdata,guidata(hObject)));
  set(handles.maxcheckbox(i),'position',cbp);
end

% initialize scrollbar
if handles.ncheckshow >= handles.nproperties,
  delete(handles.minscrollbar);
  delete(handles.maxscrollbar);
else
  nsteps = handles.nproperties-handles.ncheckshow;
  set(handles.minslider,'min',0,...
    'max',nsteps,'sliderstep',[1/nsteps,5/nsteps],'value',nsteps);
  set(handles.maxslider,'min',0,...
    'max',nsteps,'sliderstep',[1/nsteps,5/nsteps],'value',nsteps);
end

% initialize bounds on r
set(handles.minrtext,'string',num2str(handles.minr));
set(handles.maxrtext,'string',num2str(handles.maxr));

function updatewindow(handles)

updateeditable(handles);
for i = 1:handles.nconditions-1,
  updatenoteditable(handles,i);
end

function updateeditable(handles)

i = handles.editing;

% title of description
set(handles.editablepanel,'title',handles.conditiontitle{i});

% description of condition
set(handles.editablemaintext,'String',handles.maintext{i});

% selected properties
minormaxslider_Callback(handles.editablepanel,[],handles);
handles = guidata(handles.editablepanel);

updatepropertybox(handles);

function updatepropertybox(handles)

i = handles.editing;
% last selected property
j = handles.lastproperty(i);
ismax = handles.lastpropertytype(i);
pr = handles.properties{j};

% title of panel
set(handles.propertydescrpanel,'title',pr);

% description
if isfield(handles.propdescrs,pr),
  de = handles.propdescrs.(pr);
elseif isfield(handles,'units') && isfield(handles.units,pr),
  de = ['[[ ',unitsstring(handles,pr),' ]]'];
else
  de = '[[No desription available.]]';
end
set(handles.propertydescrtext,'String',de);

% lower bound
if ismax,
  v = handles.maxlower(j,i);
else
  v = handles.minlower(j,i);
end
set(handles.lowerbnd,'string',num2str(v));

% property name again
if ismax,
  v = {'max',pr};
else
  v = {'min',pr};
end
set(handles.boundstext,'string',v);

% upper bound
if ismax,
  v = handles.maxupper(j,i);
else
  v = handles.minupper(j,i);
end
set(handles.upperbnd,'string',num2str(v));

function s = unitsstring(handles,pr)

s = '';
if ~isfield(handles,'units') || ~isfield(handles.units,pr),
  return;
end
if isempty(handles.units.(pr).den) && isempty(handles.units.(pr).num),
  s = 'units';
elseif isempty(handles.units.(pr).den),
  s = sprintf('%s ',handles.units.(pr).num{:});
elseif isempty(handles.units.(pr).num),
  s = ['1 / (',sprintf('%s ',handles.units.(pr).den{:}),')'];
else
  s = ['(',sprintf('%s ',handles.units.(pr).num{:}),') / ( ',...
    sprintf('%s ',handles.units.(pr).den{:}),')'];
end

s = strrep(s,'rad','deg');

function updatenoteditable(handles,k)

hpanel = handles.(sprintf('panel%d',k));
hmaintext = handles.(sprintf('maintext%d',k));
hproperties = handles.(sprintf('properties%d',k));

i = handles.notediting(k);

set(hpanel,'title',handles.conditiontitle{i});
set(hmaintext,'string',handles.maintext{i});

pr = {};
for j = find((handles.minproperties(:,i)|handles.maxproperties(:,i))'),
  prc = handles.properties{j};
  if handles.minproperties(j,i) && handles.maxproperties(j,i),
    prc = [prc,' [max and min]'];
  elseif handles.minproperties(j,i),
    prc = [prc,' [min only]'];
  else
    prc = [prc,' [max only]'];
  end
  pr{end+1} = prc;
end
set(hproperties,'string',pr);


% --- Executes on slider movement.
function minslider_Callback(hObject, eventdata, handles)
% hObject    handle to minslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
minormaxslider_Callback(hObject,eventdata,handles);

function minormaxslider_Callback(hObject,eventdata,handles)

if hObject == handles.minslider || hObject == handles.maxslider,
  pos = round(get(hObject,'Value'));
else
  pos = round(get(handles.minslider,'value'));
end

nsteps = get(handles.minslider,'Max');
nfromtop = nsteps - pos;
handles.checkboxprops = nfromtop+(1:handles.ncheckshow);
i = handles.editing;
for j = 1:handles.ncheckshow,
  st = handles.properties{j+nfromtop};
  set(handles.mincheckbox(j),'string',st,...
    'value',handles.minproperties(j+nfromtop,i));
  set(handles.maxcheckbox(j),'string',st,...
    'value',handles.maxproperties(j+nfromtop,i));
end

if hObject == handles.maxslider,
  set(handles.minslider,'value',pos);
end
if hObject == handles.minslider,
  set(handles.maxslider,'value',pos);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function maxslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in mincheckbox2.
function mincheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to mincheckbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mincheckbox2
i = handles.editing;
j = find(handles.mincheckbox == hObject);
k = handles.checkboxprops(j);
handles.minproperties(k,i) = get(hObject,'Value');
handles.lastproperty(i) = k;
handles.lastpropertytype(i) = false;
updatepropertybox(handles);
guidata(hObject,handles);

function maxcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to mincheckbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mincheckbox2
i = handles.editing;
j = find(handles.maxcheckbox == hObject);
k = handles.checkboxprops(j);
handles.maxproperties(k,i) = get(hObject,'Value');
handles.lastproperty(i) = k;
handles.lastpropertytype(i) = true;
updatepropertybox(handles);
guidata(hObject,handles);

% --- Executes on slider movement.
function maxslider_Callback(hObject, eventdata, handles)
% hObject    handle to maxslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minormaxslider_Callback(hObject,eventdata,handles);


% --- Executes on button press in mincheckbox1.
function mincheckbox1_Callback(hObject, eventdata, handles)
% hObject    handle to mincheckbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mincheckbox1


% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

keyboard;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
%delete(hObject);



function maxrtext_Callback(hObject, eventdata, handles)
% hObject    handle to maxrtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxrtext as text
%        str2double(get(hObject,'String')) returns contents of maxrtext as a double

v = str2double(get(hObject,'String'));
if isnan(v),
  set(hObject,'string',num2str(handles.maxr));
  return;
end
v = max(1,round(v));
if v < handles.minr,
  set(hObject,'string',num2str(handles.maxr));
  return;
end
handles.maxr = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function maxrtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxrtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minrtext_Callback(hObject, eventdata, handles)
% hObject    handle to minrtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minrtext as text
%        str2double(get(hObject,'String')) returns contents of minrtext as a double

v = str2double(get(hObject,'String'));
if isnan(v),
  set(hObject,'string',num2str(handles.minr));
  return;
end
v = max(1,round(v));
if v > handles.maxr,
  set(hObject,'string',num2str(handles.minr));
  return;
end
handles.minr = v;
guidata(hObject,handles);
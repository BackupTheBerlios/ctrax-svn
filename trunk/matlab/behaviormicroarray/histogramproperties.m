function varargout = histogramproperties(varargin)
% HISTOGRAMPROPERTIES M-file for histogramproperties.fig
%      HISTOGRAMPROPERTIES, by itself, creates a new HISTOGRAMPROPERTIES or raises the existing
%      singleton*.
%
%      H = HISTOGRAMPROPERTIES returns the handle to a new HISTOGRAMPROPERTIES or the handle to
%      the existing singleton*.
%
%      HISTOGRAMPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HISTOGRAMPROPERTIES.M with the given input arguments.
%
%      HISTOGRAMPROPERTIES('Property','Value',...) creates a new HISTOGRAMPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before histogramproperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to histogramproperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help histogramproperties

% Last Modified by GUIDE v2.5 19-Dec-2008 08:51:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @histogramproperties_OpeningFcn, ...
                   'gui_OutputFcn',  @histogramproperties_OutputFcn, ...
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


% --- Executes just before histogramproperties is made visible.
function histogramproperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to histogramproperties (see VARARGIN)

% inputs:
% trx
% [params]

% Choose default command line output for histogramproperties
handles.output = hObject;
handles.trx = varargin{1};
if length(varargin) > 1,
  params0 = varargin{2};
else
  params0 = struct;
end

% initialize
handles = initializeparams(handles,params0);
handles = initializegui(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes histogramproperties wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = initializeparams(handles,params)

handles.allpropnames = fieldnames(handles.trx);
ignorefns = {'x','y','a','b','id','moviename','matname','annname','firstframe',...
  'arena','f2i','nframes','endframe','dx','dy','v','fps','pxpermm'};
handles.allpropnames = setdiff(handles.allpropnames,ignorefns);
handles.allpropnames{end+1} = 'duration';

handles.paramsvars = {'nprops','isbehavior','prop1idx',...
  'prop2idx','nbins1','nbins2','lb1','ub1','lb2','ub2','segfile',...
  'averagingidx','doinvert'};

handles.defaultsfile = which('histogramproperties');
handles.defaultsfile = strrep(handles.defaultsfile,'histogramproperties.m','.histogrampropertiesrc.mat');

% set hard-coded defaults
handles.nprops = 1;
handles.isbehavior = false;
handles.prop1idx = 1;
handles.prop2idx = 1;
handles.nbins1 = 50;
handles.nbins2 = 50;
handles.lb1 = struct;
handles.ub1 = struct;
handles.lb2 = struct;
handles.ub2 = struct;
handles.segfile = '';
handles.averagingidx = 1;
handles.doinvert = false;

% read last-used values
if exist(handles.defaultsfile,'file'),
  tmp = load(handles.defaultsfile);
  fns = fieldnames(tmp);
  for i = 1:length(fns),
    fn = fns{i};
    handles.(fn) = tmp.(fn);
  end
end

% set from input
for i = 1:length(handles.paramsvars),
  fn = handles.paramsvars{i};
  if isfield(params,fn),
    handles.(fn) = params.(fn);
  end
end

% set isbehavior to false if we don't have a file yet
if ~exist(handles.segfile,'file'),
  fprintf('Default seg file %s does not exist\n',handles.segfile);
  handles.isbehavior = false;
end

% load in the behavior file
if handles.isbehavior,
   segcurr = load(handles.segfile);
   % check seg variable
  if ~isfield(segcurr,'seg') || length(handles.trx) ~= length(segcurr.seg),
    fprintf('Could not parse default seg file %s\n',handles.segfile);
    handles.isbehavior = false;
  else
    % store seg & compute duration
    for i = 1:length(handles.trx),
      handles.trx(i).seg = segcurr.seg(i);
      handles.trx(i).duration = (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
    end
  end
end

function handles = initializegui(handles)

% number of properties
set(handles.onepropbutton,'value',double(handles.nprops==1));
set(handles.twopropbutton,'value',double(handles.nprops==2));

% if only one property, disable second property panel
if handles.nprops==1,
  set(handles.prop2panel,'visible','off');
else
  set(handles.prop2panel,'visible','on');
end

% properties possible
set(handles.prop1menu,'string',handles.allpropnames);
set(handles.prop2menu,'string',handles.allpropnames);

% which property
set(handles.prop1menu,'value',handles.prop1idx);
set(handles.prop2menu,'value',handles.prop2idx);

% number of bins
set(handles.nbins1edit,'string',num2str(handles.nbins1));
set(handles.nbins2edit,'string',num2str(handles.nbins2));

% initialize range
handles = setrange(handles);

% whether we only want to count during a behavior
set(handles.behavior1box,'value',handles.isbehavior);
if handles.isbehavior,
  set(handles.behavior1panel,'visible','on');
else
  set(handles.behavior1panel,'visible','off');
end

% seg name
set(handles.segfile1edit,'string',handles.segfile);

% type of averaging
set(handles.averaging1menu,'value',handles.averagingidx);

% invert segmentation
set(handles.invertbehavior1,'value',handles.doinvert);

% data output
handles.histstuff = struct;

function handles = setrange(handles)

fn = handles.allpropnames{handles.prop1idx};
if ~isfield(handles.lb1,fn),
  tmp = [handles.trx.(fn)];
  handles.lb1.(fn) = prctile(double(tmp(:)),1);
  handles.ub1.(fn) = prctile(double(tmp(:)),99);
end

if isanglename(fn),
  lb = handles.lb1.(fn)*180/pi;
  ub = handles.ub1.(fn)*180/pi;
else
  lb = handles.lb1.(fn);
  ub = handles.ub1.(fn);
end
set(handles.lb1edit,'string',num2str(lb));
set(handles.ub1edit,'string',num2str(ub));

fn = handles.allpropnames{handles.prop2idx};

if ~isfield(handles.lb2,fn),
  tmp = [handles.trx.(fn)];
  handles.lb2.(fn) = prctile(tmp(:),1);
  handles.ub2.(fn) = prctile(tmp(:),99);
end

if isanglename(fn),
  lb = handles.lb2.(fn)*180/pi;
  ub = handles.ub2.(fn)*180/pi;
else
  lb = handles.lb2.(fn);
  ub = handles.ub2.(fn);
end
set(handles.lb2edit,'string',num2str(lb));
set(handles.ub2edit,'string',num2str(ub));

function savedefaults(handles)

save(handles.defaultsfile,'-struct','handles',handles.paramsvars{:});

% --- Outputs from this function are returned to the command line.
function varargout = histogramproperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.histstuff;


% --- Executes on selection change in prop1menu.
function prop1menu_Callback(hObject, eventdata, handles)
% hObject    handle to prop1menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prop1menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        prop1menu
tmp = get(hObject,'Value');
if strcmpi(handles.allpropnames{tmp},'duration') && ~handles.isbehavior,
  set(hObject,'value',handles.prop1idx);
  msgbox('Duration can only be selected if "During Behavior" is selected');
else
  handles.prop1idx = tmp;
end
handles = setrange(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function prop1menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prop1menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbins1edit_Callback(hObject, eventdata, handles)
% hObject    handle to nbins1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbins1edit as text
%        str2double(get(hObject,'String')) returns contents of nbins1edit as a double

tmp = str2double(get(hObject,'String'));
if isnan(tmp) || round(tmp)~= tmp || tmp <= 0,
  set(hObject,'string',num2str(handles.nbins1));
else
  handles.nbins1 = tmp;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nbins1edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbins1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lb1edit_Callback(hObject, eventdata, handles)
% hObject    handle to lb1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lb1edit as text
%        str2double(get(hObject,'String')) returns contents of lb1edit as a double
tmp = str2double(get(hObject,'String'));
fn = handles.allpropnames{handles.prop1idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.lb1.(fn)));
else
  % convert back to radians
  if isanglename(fn),
    tmp = tmp *pi/180;
  end
  if tmp > handles.ub1.(fn),
    tmp = handles.ub1.(fn);
    set(hObject,'string',num2str(tmp));
  end
  handles.lb1.(fn) = tmp;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lb1edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ub1edit_Callback(hObject, eventdata, handles)
% hObject    handle to ub1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ub1edit as text
%        str2double(get(hObject,'String')) returns contents of ub1edit as a double

tmp = str2double(get(hObject,'String'));
fn = handles.allpropnames{handles.prop1idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.ub1.(fn)));
else
  % convert back to radians
  if isanglename(fn),
    tmp = tmp *pi/180;
  end
  if tmp < handles.lb1.(fn),
    tmp = handles.lb1.(fn);
    set(hObject,'string',num2str(tmp));
  end
  handles.ub1.(fn) = tmp;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ub1edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ub1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in behavior1box.
function behavior1box_Callback(hObject, eventdata, handles)
% hObject    handle to behavior1box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of behavior1box

v = get(hObject,'Value');
if v,
  % get seg file
  if ~exist(handles.segfile,'file'),
    [handles,v] = getsegfile(handles);
  end
  if v,
    set(handles.behavior1panel,'visible','on');
    handles.isbehavior = true;
  else
    set(hObject,'value',0);
  end
else
  set(handles.behavior1panel,'visible','off');
  handles.isbehavior = false;
end
guidata(hObject,handles);

function [handles,issegfile] = getsegfile(handles)

issegfile = false;
while true,
  [matname,matpath] = uigetfile('.mat','Choose segmentation file',handles.segfile);
  if isnumeric(matname) && matname == 0,
    return;
  end
  matname = [matpath,matname];
  if ~exist(matname,'file'),
    msgbox(sprintf('File %s does not exist',matname));
    continue;
  end
  segcurr = load(matname);
  if ~isfield(segcurr,'seg'),
    msgbox(sprintf('File %s does not contain variable seg',matname));
    continue;
  end
  if length(handles.trx) ~= length(segcurr.seg)
    msgbox(sprintf('Number of flies in trx = %d, number of flies in seg = %d',...
      length(handles.trx),length(segcurr.seg)));
    continue;
  end
  handles.segfile = matname;
  break;
end
issegfile = true;
for i = 1:length(handles.trx),
  handles.trx(i).seg = segcurr.seg(i);
  handles.trx(i).duration = (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
end
set(handles.segfile1edit,'string',handles.segfile);  

function segfile1edit_Callback(hObject, eventdata, handles)
% hObject    handle to segfile1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segfile1edit as text
%        str2double(get(hObject,'String')) returns contents of segfile1edit as a double
[handles,issegfile] = getsegfile(handles);
set(handles.segfile1edit,'string',handles.segfile);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function segfile1edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segfile1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in averaging1menu.
function averaging1menu_Callback(hObject, eventdata, handles)
% hObject    handle to averaging1menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns averaging1menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from averaging1menu
handles.averagingidx = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function averaging1menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to averaging1menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in prop2menu.
function prop2menu_Callback(hObject, eventdata, handles)
% hObject    handle to prop2menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prop2menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prop2menu

tmp = get(hObject,'Value');
if strcmpi(handles.allpropnames{tmp},'duration') && ~handles.isbehavior2,
  set(hObject,'value',handles.prop2idx);
  msgbox('Interval duration can only be selected if "During Behavior" is selected');
else
  handles.prop2idx = tmp;
end
handles = setrange(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function prop2menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prop2menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbins2edit_Callback(hObject, eventdata, handles)
% hObject    handle to nbins2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbins2edit as text
%        str2double(get(hObject,'String')) returns contents of nbins2edit as a double
tmp = str2double(get(hObject,'String'));
if isnan(tmp) || round(tmp)~= tmp || tmp <= 0,
  set(hObject,'string',num2str(handles.nbins2));
else
  handles.nbins2 = tmp;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nbins2edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbins2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lb2edit_Callback(hObject, eventdata, handles)
% hObject    handle to lb2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lb2edit as text
%        str2double(get(hObject,'String')) returns contents of lb2edit as a double

tmp = str2double(get(hObject,'String'));
fn = handles.allpropnames{handles.prop2idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.lb2.(fn)));
else
  % convert back to radians
  if isanglename(fn),
    tmp = tmp *pi/180;
  end
  if tmp > handles.ub2.(fn),
    tmp = handles.ub2.(fn);
    set(hObject,'string',num2str(tmp));
  end
  handles.lb2.(fn) = tmp;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lb2edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ub2edit_Callback(hObject, eventdata, handles)
% hObject    handle to ub2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ub2edit as text
%        str2double(get(hObject,'String')) returns contents of ub2edit as a double

tmp = str2double(get(hObject,'String'));
fn = handles.allpropnames{handles.prop2idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.ub2.(fn)));
else
  % convert back to radians
  if isanglename(fn),
    tmp = tmp *pi/180;
  end
  if tmp < handles.lb2.(fn),
    tmp = handles.lb2.(fn);
    set(hObject,'string',num2str(tmp));
  end
  handles.ub2.(fn) = tmp;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ub2edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ub2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in updatebutton.
function updatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% compute behavior interval averages if necessary

nflies = length(handles.trx);
handles.x = cell(1,nflies);

% properties we're histogramming
props = cell(1,handles.nprops);
props{1} = handles.allpropnames{handles.prop1idx};
if handles.nprops == 2,
  props{2} = handles.allpropnames{handles.prop2idx};
end
% number of bins
nbins = zeros(1,handles.nprops);
nbins(1) = handles.nbins1;
if handles.nprops == 2,
  nbins(2) = handles.nbins2;
end
% lower and upper bounds
lb = zeros(1,handles.nprops);
ub = zeros(1,handles.nprops);
lb(1) = handles.lb1.(handles.allpropnames{handles.prop1idx});
ub(1) = handles.ub1.(handles.allpropnames{handles.prop1idx});
if handles.nprops == 2,
  lb(2) = handles.lb2.(handles.allpropnames{handles.prop2idx});
  ub(2) = handles.ub2.(handles.allpropnames{handles.prop2idx});
end

avprops = get(handles.averaging1menu,'string');
if handles.isbehavior,
  averaging = avprops{handles.averagingidx};
else
  averaging = '';
end
[handles.data,handles.histstuff] = extractdata(handles.trx,props,nbins,lb,ub,handles.isbehavior,averaging,handles.doinvert);

if ~isfield(handles,'plotstuff'),
  handles.plotstuff.fig = figure;
end
% get plot mode
contents = get(handles.plotstatisticmenu,'string');
plotmode = contents{get(handles.plotstatisticmenu,'value')};
plotindivs = get(handles.plotindivsbox,'value') == 1;
contents = get(handles.plotstdmenu,'string');
plotstd = contents{get(handles.plotstdmenu,'value')};

handles.plotstuff = plothistogram(handles.histstuff,'fighandles',handles.plotstuff,...
  'plotmode',plotmode,'plotindivs',plotindivs,'plotstd',plotstd);

guidata(hObject,handles);

% --- Executes on button press in exportbutton.
function exportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in onepropbutton.
function onepropbutton_Callback(hObject, eventdata, handles)
% hObject    handle to onepropbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onepropbutton
v = get(hObject,'Value');
if v == 1,
  handles.nprops = 1;
  set(handles.prop2panel,'visible','off');
  set(handles.twopropbutton,'value',0)
else
  handles.nprops = 2;
  set(handles.prop2panel,'visible','on');
  set(handles.twopropbutton,'value',1)
end
guidata(hObject,handles);


% --- Executes on button press in twopropbutton.
function twopropbutton_Callback(hObject, eventdata, handles)
% hObject    handle to twopropbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twopropbutton

% Hint: get(hObject,'Value') returns toggle state of onepropbutton
v = get(hObject,'Value');
if v == 0,
  handles.nprops = 1;
  set(handles.prop2panel,'visible','off');
  set(handles.onepropbutton,'value',1)
else
  handles.nprops = 2;
  set(handles.prop2panel,'visible','on');
  set(handles.onepropbutton,'value',0)
end
guidata(hObject,handles);


% --- Executes on button press in invertbehavior1.
function invertbehavior1_Callback(hObject, eventdata, handles)
% hObject    handle to invertbehavior1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertbehavior1

handles.doinvert = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on selection change in plotstatisticmenu.
function plotstatisticmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotstatisticmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotstatisticmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotstatisticmenu


% --- Executes during object creation, after setting all properties.
function plotstatisticmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotstatisticmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotstdmenu.
function plotstdmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotstdmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotstdmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotstdmenu


% --- Executes during object creation, after setting all properties.
function plotstdmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotstdmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotindivsbox.
function plotindivsbox_Callback(hObject, eventdata, handles)
% hObject    handle to plotindivsbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotindivsbox


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

% save defaults to file
savedefaults(handles);

delete(hObject);
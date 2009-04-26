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

% Last Modified by GUIDE v2.5 23-Apr-2009 08:30:45

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
setuppath;
handles.output = hObject;
if length(varargin) < 1,
  fprintf('Select mat file containing trajectories with per-frame properties\n');
  [handles.trx,handles.matname,loadsucceeded] = load_tracks();
  if ~loadsucceeded,
    error('Could not load tracks');
  end
else
  handles.trx = varargin{1};
  if isfield(handles.trx,'matname'),
    handles.matname = handles.trx(1).matname;
  end
end

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

function handles = setpropnames(handles)

% store previous values
prop1idx = get(handles.prop1menu,'Value');
s = get(handles.prop1menu,'string');
prop1name = s{prop1idx};
prop2idx = get(handles.prop2menu,'Value');
prop2name = s{prop2idx};

handles.allpropnames = fieldnames(handles.trx);
ignoreidx = [];
for i = 1:length(handles.allpropnames),
  for j = 1:length(handles.trx),
    expectedlength = handles.trx(j).nframes-1;
    if numel(handles.trx(j).(handles.allpropnames{i})) < expectedlength,
      ignoreidx(end+1) = i;
      break;
    end
  end
end
handles.allpropnames(ignoreidx) = [];
ignorefns = {'x','y','a','b'};
handles.allpropnames = setdiff(handles.allpropnames,ignorefns);
handles.allpropnames{end+1} = 'duration';

% set idx to previous values
set(handles.prop1menu,'string',handles.allpropnames);
set(handles.prop2menu,'string',handles.allpropnames);
handles.prop1idx = find(strcmp(handles.allpropnames,prop1name),1);
if isempty(handles.prop1idx),
  handles.prop1idx = 1;
end
handles.prop2idx = find(strcmp(handles.allpropnames,prop2name),1);
if isempty(handles.prop2idx),
  handles.prop2idx = 1;
end
set(handles.prop1menu,'value',handles.prop1idx);
set(handles.prop2menu,'value',handles.prop2idx);

function handles = initializeparams(handles,params)

handles.allpropnames = fieldnames(handles.trx);
ignoreidx = [];
for i = 1:length(handles.allpropnames),
  for j = 1:length(handles.trx),
    expectedlength = handles.trx(j).nframes-1;
    if numel(handles.trx(j).(handles.allpropnames{i})) < expectedlength,
      ignoreidx(end+1) = i;
      break;
    end
  end
end
handles.allpropnames(ignoreidx) = [];
ignorefns = {'x','y','a','b'};
handles.allpropnames = setdiff(handles.allpropnames,ignorefns);
handles.allpropnames{end+1} = 'duration';

handles.paramsvars = {'nprops','isbehavior','nbins1','nbins2','lb1','ub1',...
                    'lb2','ub2','segfile','averagingidx','doinvert','prop1name','prop2name',...
                    'transform1','transform2','rangeunits1','rangeunits2','flytype',...
                    'dataname','savename'};

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
handles.segfile = {''};
handles.averagingidx = 1;
handles.doinvert = false;
handles.transform1 = 'None (Identity)';
handles.transform2 = 'None (Identity)';
handles.rangeunits1 = struct;
handles.rangeunits2 = struct;
handles.savename = '';

% get names of fly types
if isfield(handles.trx,'type') && isfield(handles.trx,'sex'),
  for i = 1:length(handles.trx),
    handles.trx(i).type = sprintf('%s, %s',handles.trx(i).type,handles.trx(i).sex);
  end
elseif isfield(handles.trx,'sex') && ~isfield(handles.trx,'type'),
  for i = 1:length(handles.trx),
    handles.trx(i).type = handles.trx(i).sex;
  end
end
if isfield(handles.trx,'type'),
  handles.alltypes = unique({handles.trx.type});
else
  handles.alltypes = {};
end
handles.alltypes = [{'All'},handles.alltypes];

% initialize fly type to be all
handles.flytype = {{'All'}};
handles.dataname = {'Data 1'};
handles.datashow = 1;

% read last-used values
if exist(handles.defaultsfile,'file'),
  tmp = load(handles.defaultsfile);
  fns = fieldnames(tmp);
  nosecondprop = false;
  for i = 1:length(fns),
    fn = fns{i};
    if strcmpi(fn,'prop1name'),
      handles.prop1idx = find(strcmpi(tmp.prop1name,handles.allpropnames),1);
      if isempty(handles.prop1idx), 
        handles.prop1idx = 1;
      end
    elseif strcmpi(fn,'prop2name'),
      handles.prop2idx = find(strcmpi(tmp.prop2name,handles.allpropnames),1);
      if isempty(handles.prop2idx), 
        handles.prop2idx = 1;
        nosecondprop = true;
      end
    elseif strcmpi(fn,'prop1idx'),
      continue;
    elseif strcmpi(fn,'prop2idx'),
      continue;
    else
      handles.(fn) = tmp.(fn);
    end
  end
  
  % reset nprops to 1 if couldn't figure out the second property
  if nosecondprop,
    handles.nprops = 1;
  end
end

% set from input
for i = 1:length(handles.paramsvars),
  fn = handles.paramsvars{i};
  if isfield(params,fn),
    handles.(fn) = params.(fn);
  end
end

% check fly type to make sure it is allowed
for i = 1:length(handles.flytype),
  badidx = [];
  badidx = ~ismember(handles.flytype{i},handles.alltypes);
  handles.flytype{i}(badidx) = [];
  if isempty(handles.flytype{i}) || ismember('All',handles.flytype{i}),
    handles.flytype{i} = {'All'};
  end
end

% set isbehavior to false if we don't have a file yet
for i = 1:length(handles.segfile),
  if ~exist(handles.segfile{i},'file'),
    if ~isempty(handles.segfile{i}),
      fprintf('Default seg file >%s< does not exist\n',handles.segfile{i});
    end
    handles.isbehavior(i) = false;
  end
end

% check for and remove duplicate data types
remove = [];
for i = 1:length(handles.segfile),
  for j = i+1:length(handles.segfile),
    if (handles.isbehavior(i) == handles.isbehavior(j)) && ...
       (~handles.isbehavior(i) || strcmpi(handles.segfile{i},handles.segfile{j})) && ...
       (handles.averagingidx(i) == handles.averagingidx(j)) && ...
       (handles.doinvert(i) == handles.doinvert(j)),
      % check flytype
      tmp1 = handles.flytype{i};
      tmp2 = handles.flytype{j};
      if ismember('All',tmp1),
        tmp1 = union(tmp1,setdiff(handles.alltypes,{'All'}));
      end
      if ismember('All',tmp2),
        tmp2 = union(tmp2,setdiff(handles.alltypes,{'All'}));
      end
      issame = isempty(setxor(tmp1,tmp2));
      if issame, 
        remove(end+1) = j;
        fprintf('Fly types for %s = ',handles.dataname{i});
        fprintf('%s ',handles.flytype{i}{:});
        fprintf('\nFly types for %s = ',handles.dataname{j});
        fprintf('%s ',handles.flytype{j}{:});
        fprintf('Duplicate, removing %s\n',handles.dataname{j});
      end
    end
  end
end

handles.segfile(remove) = [];
handles.isbehavior(remove) = [];
handles.averagingidx(remove) = [];
handles.doinvert(remove) = [];
handles.flytype(remove) = [];
handles.dataname(remove) = [];

handles.ndata = length(handles.dataname);

% load in the behavior file
for d = 1:handles.ndata,

  if handles.isbehavior(d),
    segcurr = load(handles.segfile{d});
    % check seg variable
    if ~isfield(segcurr,'seg') || length(handles.trx) ~= length(segcurr.seg),
      fprintf('Could not parse default seg file %s for %s\n',handles.segfile{d},handles.dataname{d});
      handles.isbehavior(d) = false;
    else
      % store seg & compute duration
      for i = 1:length(handles.trx),
        handles.trx(i).seg{d} = segcurr.seg(i);
        handles.trx(i).duration{d} = (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
      end
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

% transform to data
contents = get(handles.prop1transformmenu,'string');
i = find(strcmpi(handles.transform1,contents),1);
if isempty(i),
  i = 1;
end
set(handles.prop1transformmenu,'value',i);
contents = get(handles.prop2transformmenu,'string');
i = find(strcmpi(handles.transform2,contents),1);
if isempty(i),
  i = 1;
end
set(handles.prop2transformmenu,'value',i);

% initialize range
handles = setrange(handles);

% data box
set(handles.flytypebox,'string',handles.alltypes);
set(handles.datamenu,'string',[handles.dataname,{'Add new data'}]);
handles = initializedatabox(handles);

% data output
handles.histstuff = struct;

function handles = initializedatabox(handles)

% whether we only want to count during a behavior
set(handles.behavior1box,'value',handles.isbehavior(handles.datashow));
if handles.isbehavior(handles.datashow),
  set(handles.behavior1panel,'visible','on');
else
  set(handles.behavior1panel,'visible','off');
end

% seg name
[paths,files] = split_path_and_filename(handles.segfile{handles.datashow});
set(handles.segfile1edit,'string',files,'tooltipstring',handles.segfile{handles.datashow});  

% type of averaging
set(handles.averaging1menu,'value',handles.averagingidx(handles.datashow));

% invert segmentation
set(handles.invertbehavior1,'value',handles.doinvert(handles.datashow));

% fly type
v = [];
for i = 1:length(handles.flytype{handles.datashow}),
  j = find(strcmpi(handles.flytype{handles.datashow}{i},handles.alltypes));
  v = [v,j];
end
if isempty(v),
  v = 1;
end
set(handles.flytypebox,'value',v);

% data names
set(handles.datamenu,'value',handles.datashow);
set(handles.datanameedit,'string',handles.dataname{handles.datashow});

% can't delete only data line
if handles.ndata == 1,
  set(handles.deletebutton,'enable','off');
else
  set(handles.deletebutton,'enable','on');
end

function handles = setrange(handles)

fn = handles.allpropnames{handles.prop1idx};
fn1 = getboundname1(handles);

if ~isfield(handles.lb1,fn1),
  handles.lb1.(fn1) = 1;
  handles.ub1.(fn1) = 99;
  handles.rangeunits1.(fn1) = 'Percent';
end

if any(strcmpi(handles.trx(1).units.(fn).num,'rad')) && ...
      strcmpi(handles.rangeunits1.(fn1),'units'),
  if strcmpi(handles.transform1,'log absolute value'),
    lb = handles.lb1.(fn1) + log(180/pi);
    ub = handles.ub1.(fn1) + log(180/pi);
  else
    lb = handles.lb1.(fn1)*180/pi;
    ub = handles.ub1.(fn1)*180/pi;
  end
else
  lb = handles.lb1.(fn1);
  ub = handles.ub1.(fn1);
end
set(handles.lb1edit,'string',num2str(lb));
set(handles.ub1edit,'string',num2str(ub));
contents = get(handles.rangeunits1menu,'string');
i = find(strcmpi(contents,handles.rangeunits1.(fn1)),1);
if isempty(i), i = 1; end
set(handles.rangeunits1menu,'value',i);

fn = handles.allpropnames{handles.prop2idx};
fn1 = getboundname2(handles);

if ~isfield(handles.lb2,fn1),  
  handles.lb2.(fn1) = 1;
  handles.ub2.(fn1) = 99;
  handles.rangeunits2.(fn1) = 'Percent';
end

if any(strcmpi(handles.trx(1).units.(fn).num,'rad')) && ... 
      strcmpi(handles.rangeunits2.(fn1),'units'),
   if strcmpi(handles.transform2,'log absolute value'),
    lb = handles.lb2.(fn1) + log(180/pi);
    ub = handles.ub2.(fn1) + log(180/pi);
   else
    lb = handles.lb2.(fn1)*180/pi;
    ub = handles.ub2.(fn1)*180/pi;
   end
else
  lb = handles.lb2.(fn1);
  ub = handles.ub2.(fn1);
end
set(handles.lb2edit,'string',num2str(lb));
set(handles.ub2edit,'string',num2str(ub));
contents = get(handles.rangeunits2menu,'string');
i = find(strcmpi(contents,handles.rangeunits2.(fn1)),1);
if isempty(i), i = 1; end
set(handles.rangeunits2menu,'value',i);

function savedefaults(handles)

handles.prop2name = handles.allpropnames{handles.prop2idx};
handles.prop1name = handles.allpropnames{handles.prop1idx};
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
if strcmpi(handles.allpropnames{tmp},'duration') && ~all(handles.isbehavior),
  set(hObject,'value',handles.prop1idx);
  msgbox('Duration can only be selected if "During Behavior" is selected for all data');
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
fn = getboundname1(handles);
fntrx = handles.allpropnames{handles.prop1idx};

if isnan(tmp),
  set(hObject,'string',num2str(handles.lb1.(fn)));
else
  % convert back to radians
  if any(strcmpi(handles.trx(1).units.(fntrx).num,'rad')) && ...
        strcmpi(handles.rangeunits1.(fn),'units'),
    if strcmpi(handles.transform1,'log absolute value'),
      tmp = tmp + log(pi/180);
    else
      tmp = tmp * pi/180;
    end
  end
  if tmp > handles.ub1.(fn),
    tmp = handles.ub1.(fn);
    set(hObject,'string',num2str(tmp));
  end
  if strcmpi(handles.rangeunits1.(fn),'percent'),
    tmp = max(min(tmp,100),0);
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
fn = getboundname1(handles);
fntrx = handles.allpropnames{handles.prop1idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.ub1.(fn)));
else
  % convert back to radians
  
  if any(strcmpi(handles.trx(1).units.(fntrx).num,'rad')) && ... 
        strcmpi(handles.rangeunits1.(fn),'units'),
    if strcmpi(handles.transform1,'log absolute value'),
      tmp = tmp + log(pi/180);
    else
      tmp = tmp * pi/180;
    end
  end
  if tmp < handles.lb1.(fn),
    tmp = handles.lb1.(fn);
    set(hObject,'string',num2str(tmp));
  end
  if strcmpi(handles.rangeunits1.(fn),'percent'),
    tmp = max(min(tmp,100),0);
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
  [handles,v] = getsegfile(handles);
  if v,
    set(handles.behavior1panel,'visible','on');
    handles.isbehavior(handles.datashow) = true;
  else
    set(hObject,'value',0);
  end
else
  set(handles.behavior1panel,'visible','off');
  handles.isbehavior(handles.datashow) = false;
end
guidata(hObject,handles);

function [handles,issegfile] = getsegfile(handles)

issegfile = false;
while true,
  if exist(handles.segfile{handles.datashow},'file'),
    matname = handles.segfile{handles.datashow};
  else
    [matname,matpath] = uigetfile('.mat','Choose segmentation file',handles.segfile{handles.datashow});
    if isnumeric(matname) && matname == 0,
      return;
    end
    matname = [matpath,matname];
  end
  if ~exist(matname,'file'),
    handles.segfile{handles.datashow} = '';
    uiwait(msgbox(sprintf('File %s does not exist',matname)));
    continue;
  end
  segcurr = load(matname);
  if ~isfield(segcurr,'seg'),
    handles.segfile{handles.datashow} = '';
    uiwait(msgbox(sprintf('File %s does not contain variable seg',matname)));
    continue;
  end
  if length(handles.trx) ~= length(segcurr.seg)
    handles.segfile{handles.datashow} = '';
    uiwait(msgbox(sprintf('Number of flies in trx = %d, number of flies in seg = %d',...
      length(handles.trx),length(segcurr.seg))));
    continue;
  end
  handles.segfile{handles.datashow} = matname;
  break;
end
issegfile = true;
for i = 1:length(handles.trx),
  handles.trx(i).seg{handles.datashow} = segcurr.seg(i);
  handles.trx(i).duration{handles.datashow} = ...
      (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
end
[paths,files] = split_path_and_filename(handles.segfile{handles.datashow});
set(handles.segfile1edit,'string',files,'tooltipstring',handles.segfile{handles.datashow});  

function segfile1edit_Callback(hObject, eventdata, handles)
% hObject    handle to segfile1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segfile1edit as text
%        str2double(get(hObject,'String')) returns contents of segfile1edit as a double
[handles,issegfile] = getsegfile(handles);
[paths,files] = split_path_and_filename(handles.segfile{handles.datashow});
set(handles.segfile1edit,'string',files,'tooltipstring',handles.segfile{handles.datashow});  
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
handles.averagingidx(handles.datashow) = get(hObject,'Value');
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
fn = getboundname2(handles);
fntrx = handles.allpropnames{handles.prop1idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.lb2.(fn)));
else
  % convert back to radians
  if any(strcmpi(handles.trx(1).units.(fntrx).num,'rad')) && ...
        strcmpi(handles.rangeunits2.(fn),'units'),
    if strcmpi(handles.transform2,'log absolute value'),
      tmp = tmp + log(pi/180);
    else
      tmp = tmp * pi/180;
    end
  end
  if tmp > handles.ub2.(fn),
    tmp = handles.ub2.(fn);
    set(hObject,'string',num2str(tmp));
  end
  if strcmpi(handles.rangeunits2.(fn),'percent'),
    tmp = max(min(tmp,100),0);
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
fn = getboundname2(handles);
fntrx = handles.allpropnames{handles.prop1idx};
if isnan(tmp),
  set(hObject,'string',num2str(handles.ub2.(fn)));
else
  % convert back to radians
  if any(strcmpi(handles.trx(1).units.(fntrx).num,'rad')) && ...
        strcmpi(handles.rangeunits2.(fn),'units'),
    if strcmpi(handles.transform2,'log absolute value'),
      tmp = tmp + log(pi/180);
    else
      tmp = tmp * pi/180;
    end
  end
  if tmp < handles.lb2.(fn),
    tmp = handles.lb2.(fn);
    set(hObject,'string',num2str(tmp));
  end
  if strcmpi(handles.rangeunits2.(fn),'percent'),
    tmp = max(min(tmp,100),0);
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

% properties we're histogramming
props = cell(1,handles.nprops);
props{1} = handles.allpropnames{handles.prop1idx};
if handles.nprops == 2,
  props{2} = handles.allpropnames{handles.prop2idx};
end
handles.datastuff.props = props;
% number of bins
nbins = zeros(1,handles.nprops);
nbins(1) = handles.nbins1;
if handles.nprops == 2,
  nbins(2) = handles.nbins2;
end
handles.datastuff.nbins = nbins;
% lower and upper bounds
lb = zeros(1,handles.nprops);
ub = zeros(1,handles.nprops);
fn = getboundname1(handles);
fntrx = handles.allpropnames{handles.prop1idx};
if strcmpi(handles.rangeunits1.(fn),'percent'),
  tmp = [handles.trx.(fntrx)];
  switch lower(handles.transform1),
   case 'absolute value',
    tmp = abs(tmp);
   case 'log absolute value',
    tmp = log(abs(tmp));
  end
  lb(1) = myprctile(tmp(:),handles.lb1.(fn));
  ub(1) = myprctile(tmp(:),handles.ub1.(fn));
else
  lb(1) = handles.lb1.(fn);
  ub(1) = handles.ub1.(fn);
end
if handles.nprops == 2,
  fn = getboundname2(handles);
  fntrx = handles.allpropnames{handles.prop2idx};
  if strcmpi(handles.rangeunits2.(fn),'percent'),
    tmp = [handles.trx.(fntrx)];
    switch lower(handles.transform2),
     case 'absolute value',
      tmp = abs(tmp);
     case 'log absolute value',
      tmp = log(abs(tmp));
    end
    lb(2) = myprctile(tmp(:),handles.lb2.(fn));
    ub(2) = myprctile(tmp(:),handles.ub2.(fn));
  else
    lb(2) = handles.lb2.(fn);
    ub(2) = handles.ub2.(fn);
  end
end
transforms = cell(1,handles.nprops);
transforms{1} = handles.transform1;
if handles.nprops == 2,
  transforms{2} = handles.transform2;
end
handles.datastuff.lb = lb;
handles.datastuff.ub = ub;
handles.datastuff.transforms = transforms;
avprops = get(handles.averaging1menu,'string');
averaging = cell(1,handles.ndata);
for i = 1:handles.ndata,
  if handles.isbehavior(i),
    averaging{i} = avprops{handles.averagingidx(i)};
  else
    averaging{i} = '';
  end
end
handles.datastuff.isbehavior = handles.isbehavior;
handles.datastuff.segfile = handles.segfile;
handles.datastuff.averaging = averaging;
handles.datastuff.doinvert = handles.doinvert;
flytype = handles.flytype;
for i = 1:handles.ndata,
  if ismember('All',handles.flytype{i}),
    flytype{i} = setdiff(handles.alltypes,{'All'});
  end
end
handles.datastuff.flytype = flytype;
handles.datastuff.dataname = handles.dataname;
[handles.data,handles.histstuff] = extractdata(handles.trx,props,nbins,lb,ub,transforms,handles.isbehavior,averaging,handles.doinvert,flytype,handles.dataname);

if ~isfield(handles,'plotstuff'),
  handles.plotstuff = struct;
end
% get plot mode
contents = get(handles.plotstatisticmenu,'string');
plotmode = contents{get(handles.plotstatisticmenu,'value')};
plotindivs = get(handles.plotindivsbox,'value') == 1;
contents = get(handles.plotstdmenu,'string');
plotstd = contents{get(handles.plotstdmenu,'value')};
logplot = get(handles.doplotlogbox,'value') == 1;

handles.plotstuff = plothistogram(handles.histstuff,'fighandles',handles.plotstuff,...
  'plotmode',plotmode,'plotindivs',plotindivs,'plotstd',plotstd,'logplot',logplot);

guidata(hObject,handles);

% --- Executes on button press in exportbutton.
function exportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('Enter name of mat file to save histogrammed data to.\n');
[savename,savepath] = uiputfile('*.mat','Save histogrammed results',handles.savename);
if ~ischar(savename),
  return;
end
handles.savename = [savepath,savename];
if ~isfield(handles,'histstuff') || isempty(fieldnames(handles.histstuff)),
  guidata(hObject,handles);
  updatebutton_Callback(hObject, eventdata, handles);
  handles = guidata(hObject);
end
tmp = handles.histstuff;
save(handles.savename,'-struct','tmp');
tmp = handles.datastuff;
save('-append',handles.savename,'-struct','tmp');
%handles.prop2name = handles.allpropnames{handles.prop2idx};
%handles.prop1name = handles.allpropnames{handles.prop1idx};
%save(handles.savename,'-append','-struct','handles',handles.paramsvars{:});

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

handles.doinvert(handles.datashow) = get(hObject,'Value');
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


% --- Executes on button press in doplotlogbox.
function doplotlogbox_Callback(hObject, eventdata, handles)
% hObject    handle to doplotlogbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doplotlogbox


% --- Executes on selection change in prop2transformmenu.
function prop2transformmenu_Callback(hObject, eventdata, handles)
% hObject    handle to prop2transformmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prop2transformmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prop2transformmenu

contents = get(hObject,'string');
handles.transform2 = contents{get(hObject,'value')};
handles = setrange(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function prop2transformmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prop2transformmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in prop1transformmenu.
function prop1transformmenu_Callback(hObject, eventdata, handles)
% hObject    handle to prop1transformmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prop1transformmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prop1transformmenu
contents = get(hObject,'string');
handles.transform1 = contents{get(hObject,'value')};
handles = setrange(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function prop1transformmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prop1transformmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fn = getboundname1(handles)

fn = handles.allpropnames{handles.prop1idx};
if strcmpi(handles.transform1,'absolute value')
  fn = ['ABSVAL_',fn];
elseif strcmpi(handles.transform1,'log absolute value'),
  fn = ['LOGABSVAL',fn];
end

function fn = getboundname2(handles)

fn = handles.allpropnames{handles.prop2idx};
if strcmpi(handles.transform2,'absolute value')
  fn = ['ABSVAL_',fn];
elseif strcmpi(handles.transform2,'log absolute value'),
  fn = ['LOGABSVAL',fn];
end

function y = myprctile(x,p)
% only for 1D data
% fixes errors when we are averaging -infs, infs

y = prctile(x,p);

if ~isnan(y),
  return;
end

n = length(x);
q = min(max(1,round(n*p)),n);
y = x(q);


% --- Executes on selection change in rangeunits1menu.
function rangeunits1menu_Callback(hObject, eventdata, handles)
% hObject    handle to rangeunits1menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns rangeunits1menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rangeunits1menu

fn = handles.allpropnames{handles.prop1idx};
fn1 = getboundname1(handles);

contents = get(hObject,'String');
newunits = contents{get(hObject,'Value')};
if strcmpi(newunits,handles.rangeunits1.(fn1)),
  return;
end

% get current property data 
tmp = [handles.trx.(fn)];
tmp = tmp(:);
switch lower(handles.transform1),
 case 'absolute value',
  tmp = abs(tmp);
 case 'log absolute value',
  tmp = log(abs(tmp));
end
if strcmpi(newunits,'percent'),
  lb = (nnz(tmp<handles.lb1.(fn1))+nnz(tmp==handles.lb1.(fn1))/2)/length(tmp)*100;
  ub = 100-(nnz(tmp>handles.ub1.(fn1))+nnz(tmp==handles.ub1.(fn1))/2)/length(tmp)*100;
else
  % units
  lb = myprctile(tmp,handles.lb1.(fn1));
  ub = myprctile(tmp,handles.ub1.(fn1));
end
handles.lb1.(fn1) = lb;
handles.ub1.(fn1) = ub;
handles.rangeunits1.(fn1) = newunits;
if any(strcmpi(handles.trx(1).units.(fn).num,'rad')) && ...
      strcmpi(handles.rangeunits1.(fn1),'units'),
  if strcmpi(handles.transform2,'log absolute value'),
    lb = lb + log(180/pi);
    ub = ub + log(180/pi);
  else
    lb = lb*180/pi;
    ub = ub*180/pi;
  end
else
  lb = lb;
end
set(handles.lb1edit,'string',num2str(lb));
set(handles.ub1edit,'string',num2str(ub));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rangeunits1menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangeunits1menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in rangeunits2menu.
function rangeunits2menu_Callback(hObject, eventdata, handles)
% hObject    handle to rangeunits2menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns rangeunits2menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rangeunits2menu

fn = handles.allpropnames{handles.prop2idx};
fn1 = getboundname2(handles);

contents = get(hObject,'String');
newunits = contents{get(hObject,'Value')};
if strcmpi(newunits,handles.rangeunits2.(fn1)),
  return;
end

% get current property data 
tmp = [handles.trx.(fn)];
tmp = tmp(:);
switch lower(handles.transform2),
 case 'absolute value',
  tmp = abs(tmp);
 case 'log absolute value',
  tmp = log(abs(tmp));
end
if strcmpi(newunits,'percent'),
  lb = (nnz(tmp<handles.lb2.(fn1))+nnz(tmp==handles.lb2.(fn1))/2)/length(tmp)*100;
  ub = 100-(nnz(tmp>handles.ub2.(fn1))+nnz(tmp==handles.ub2.(fn1))/2)/length(tmp)*100;
else
  % units
  lb = myprctile(tmp,handles.lb2.(fn1));
  ub = myprctile(tmp,handles.ub2.(fn1));
end
handles.lb2.(fn1) = lb;
handles.ub2.(fn1) = ub;
handles.rangeunits2.(fn1) = newunits;

% for display, need to convert to degrees
if any(strcmpi(handles.trx(1).units.(fn).num,'rad')) && ...
      strcmpi(handles.rangeunits2.(fn1),'units'),
  if strcmpi(handles.transform2,'log absolute value'),
    lb = lb + log(180/pi);
    ub = ub + log(180/pi);
  else
    lb = lb*180/pi;
    ub = ub*180/pi;
  end
else
  lb = lb;
end
set(handles.lb2edit,'string',num2str(lb));
set(handles.ub2edit,'string',num2str(ub));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rangeunits2menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangeunits2menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in datamenu.
function datamenu_Callback(hObject, eventdata, handles)
% hObject    handle to datamenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns datamenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from datamenu

olddatashow = handles.datashow;
handles.datashow = get(hObject,'Value');
if handles.datashow > handles.ndata,
  newdata = sprintf('Data %d',handles.datashow);
  handles.dataname = [handles.dataname,{newdata}];
  set(handles.datamenu,'string',[handles.dataname,{'Add new data'}]);
  handles.ndata = handles.ndata + 1;
  handles.isbehavior(handles.ndata) = false;
  handles.doinvert(handles.ndata) = false;
  handles.averagingidx(handles.ndata) = handles.averagingidx(olddatashow);
  handles.segfile{handles.ndata} = handles.segfile{olddatashow};
  handles.flytype{handles.ndata} = {'All'};
end

handles = initializedatabox(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function datamenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datamenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in flytypebox.
function flytypebox_Callback(hObject, eventdata, handles)
% hObject    handle to flytypebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns flytypebox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from flytypebox
v = get(hObject,'value');
newtypes = cell(1,length(v));
for i = 1:length(v),
  newtypes{i} = handles.alltypes{v(i)};
end
handles.flytype{handles.datashow} = newtypes;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function flytypebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flytypebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deletebutton.
function deletebutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.ndata == 1,
  return;
end

handles.dataname(handles.datashow) = [];  
set(handles.datamenu,'string',[handles.dataname,{'Add new data'}]);
handles.isbehavior(handles.datashow) = [];
handles.doinvert(handles.datashow) = [];
handles.averagingidx(handles.datashow) = [];
handles.segfile(handles.datashow) = [];
handles.flytype(handles.ndata) = [];
handles.ndata = handles.ndata - 1;
if handles.datashow > 1,
  handles.datashow = handles.datashow - 1;
end
handles = initializedatabox(handles);
guidata(hObject,handles);


function datanameedit_Callback(hObject, eventdata, handles)
% hObject    handle to datanameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datanameedit as text
%        str2double(get(hObject,'String')) returns contents of datanameedit as a double
newname = get(hObject,'String');
if strcmpi(newname,handles.dataname{handles.datashow}),
  return;
end
if ismember(newname,handles.dataname),
  msgbox(sprintf('Data name %s already chosen, please choose a unique name',newname));
  set(hObject,'string',handles.dataname{handles.datashow});
  return;
end
handles.dataname{handles.datashow} = newname;
set(handles.datamenu,'string',[handles.dataname,{'Add new data'}]);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function datanameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datanameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over segfile1edit.
function segfile1edit_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to segfile1edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[matname,matpath] = uigetfile('*.mat','Choose behavior segmentation mat file',handles.segfile{handles.datashow});
if ~ischar(matname),
  return;
end
matname = [matpath,matname];
if ~exist(matname,'file'),
  msgbox(sprintf('File %s does not exist',matname));
  return;
end

segcurr = load(matname);
if ~isfield(segcurr,'seg'),
  msgbox(sprintf('File %s does not contain variable seg',matname));
  return;
end
if length(handles.trx) ~= length(segcurr.seg)
  msgbox(sprintf('Number of flies in trx = %d, number of flies in seg = %d',...
    length(handles.trx),length(segcurr.seg)));
  return;
end
handles.segfile{handles.datashow} = matname;

for i = 1:length(handles.trx),
  handles.trx(i).seg{handles.datashow} = segcurr.seg(i);
  handles.trx(i).duration{handles.datashow} = ...
    (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
end
[paths,files] = split_path_and_filename(handles.segfile{handles.datashow});
set(handles.segfile1edit,'string',files,'tooltipstring',handles.segfile{handles.datashow});  

guidata(hObject,handles);


% --- Executes on button press in changesegfilebutton.
function changesegfilebutton_Callback(hObject, eventdata, handles)
% hObject    handle to changesegfilebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[matname,matpath] = uigetfile('*.mat','Choose behavior segmentation mat file',handles.segfile{handles.datashow});
if ~ischar(matname),
  return;
end
matname = [matpath,matname];
if ~exist(matname,'file'),
  msgbox(sprintf('File %s does not exist',matname));
  return;
end

segcurr = load(matname);
if ~isfield(segcurr,'seg'),
  msgbox(sprintf('File %s does not contain variable seg',matname));
  return;
end
if length(handles.trx) ~= length(segcurr.seg)
  msgbox(sprintf('Number of flies in trx = %d, number of flies in seg = %d',...
    length(handles.trx),length(segcurr.seg)));
  return;
end
handles.segfile{handles.datashow} = matname;

for i = 1:length(handles.trx),
  handles.trx(i).seg{handles.datashow} = segcurr.seg(i);
  handles.trx(i).duration{handles.datashow} = ...
    (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
end
[paths,files] = split_path_and_filename(handles.segfile{handles.datashow});
set(handles.segfile1edit,'string',files,'tooltipstring',handles.segfile{handles.datashow});  

guidata(hObject,handles);


% --- Executes on button press in quitbutton.
function quitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

savedefaults(handles);

delete(handles.figure1);

% --- Executes on button press in computeperframebutton.
function computeperframebutton_Callback(hObject, eventdata, handles)
% hObject    handle to computeperframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'matname'),
  [matpath,matname] = split_path_and_filename(handles.matname);
  opts = {'matname',matname,'matpath',matpath};
else
  opts = {};
end
[loadsucceeded,savename,trx] = compute_perframe_stats_f('trx',handles.trx,opts{:});
if ~loadsucceeded,
  return;
end
handles.matname = savename;
handles.trx = trx;

handles = setpropnames(handles);

guidata(hObject,handles);


% --- Executes on button press in detectbehaviorbutton.
function detectbehaviorbutton_Callback(hObject, eventdata, handles)
% hObject    handle to detectbehaviorbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[succeeded,savenames] = detect_behaviors_f('trxnames',{handles.matname});

if ~succeeded,
  return;
end

savename = savenames{1};

olddatashow = handles.datashow;
handles.datashow = handles.ndata + 1;

newdata = sprintf('Data %d',handles.datashow);
handles.dataname = [handles.dataname,{newdata}];
set(handles.datamenu,'string',[handles.dataname,{'Add new data'}]);
handles.ndata = handles.ndata + 1;
handles.isbehavior(handles.ndata) = true;
handles.doinvert(handles.ndata) = false;
handles.averagingidx(handles.ndata) = handles.averagingidx(olddatashow);
handles.segfile{handles.ndata} = savename;
handles.flytype{handles.ndata} = {'All'};

segcurr = load(savename);
if ~isfield(segcurr,'seg'),
  handles.segfile{handles.datashow} = '';
  uiwait(msgbox(sprintf('File %s does not contain variable seg',savename)));
  return;
end
if length(handles.trx) ~= length(segcurr.seg)
  handles.segfile{handles.datashow} = '';
  uiwait(msgbox(sprintf('Number of flies in trx = %d, number of flies in seg = %d',...
    length(handles.trx),length(segcurr.seg))));
  return;
end
for i = 1:length(handles.trx),
  handles.trx(i).seg{handles.datashow} = segcurr.seg(i);
  handles.trx(i).duration{handles.datashow} = ...
    (segcurr.seg(i).t2 - segcurr.seg(i).t1 + 1)/handles.trx(i).fps;
end

%[paths,files] = split_path_and_filename(handles.segfile{handles.datashow});
%set(handles.segfile1edit,'string',files,'tooltipstring',handles.segfile{handles.datashow});  

handles = initializedatabox(handles);

guidata(hObject,handles);

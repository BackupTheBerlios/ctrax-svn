function varargout = plotbehaviormicroarray(varargin)
% PLOTBEHAVIORMICROARRAY M-file for plotbehaviormicroarray.fig
%      PLOTBEHAVIORMICROARRAY, by itself, creates a new PLOTBEHAVIORMICROARRAY or raises the existing
%      singleton*.
%
%      H = PLOTBEHAVIORMICROARRAY returns the handle to a new PLOTBEHAVIORMICROARRAY or the handle to
%      the existing singleton*.
%
%      PLOTBEHAVIORMICROARRAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTBEHAVIORMICROARRAY.M with the given input arguments.
%
%      PLOTBEHAVIORMICROARRAY('Property','Value',...) creates a new PLOTBEHAVIORMICROARRAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotbehaviormicroarray_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotbehaviormicroarray_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotbehaviormicroarray

% Last Modified by GUIDE v2.5 08-Jan-2009 16:57:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotbehaviormicroarray_OpeningFcn, ...
                   'gui_OutputFcn',  @plotbehaviormicroarray_OutputFcn, ...
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


% --- Executes just before plotbehaviormicroarray is made visible.
function plotbehaviormicroarray_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotbehaviormicroarray (see VARARGIN)

% Choose default command line output for plotbehaviormicroarray
handles.output = hObject;
setuppath;

if length(varargin) >= 1,
  params0 = varargin{1};
else
  params0 = struct;
end

handles = initializeparams(handles,params0);
handles = initializegui(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plotbehaviormicroarray wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = hardcoded_initializeparams(handles)

handles.lasttrxname = '';
handles.lastsegname = '';
handles.trxnames = {};
handles.propnamespertrx = {};
handles.allpropnames = {};
handles.trxselected = 1;
handles.behaviors = struct('name',cell(1,0),'segnames',cell(1,0));
handles.behaviorselected = 1;
%handles.segselected = 1;
handles.props.name = 'Property 1';
handles.props.behavior = 1;
handles.props.proptype = 'perframe';
handles.props.perframe = 'velmag';
handles.props.during = 'allframes';
handles.allaveraging = {'None (per-frame)','interval mean','interval median',...
                    'interval start','interval end'};
handles.props.averaging = 'None (per-frame)';
handles.propselected = 1;
handles.autochooseallownames = ...
    {'duration','fractime','freq','perframe','noaveraging','mean',...
     'median','start','end','duringbehavior',...
     'invertbehavior','allframes'};
for i = 1:length(handles.autochooseallownames),
  handles.autochoose.(handles.autochooseallownames{i}) = 'allow';
end
handles.autochoose.ignoreprops = {};
handles.autochoose.nshow = 10;
handles.autochoose.isindep = false;
handles.types.name = 'All';
handles.types.fields = {'All'};
handles.typeselected = 1;
handles.plotsettings.errorbars = 'none';
handles.plotsettings.plotindivs = false;
handles.plotsettings.normalize = true;
handles.defaultsfile = which('plotbehaviormicroarray');
handles.defaultsfile = strrep(handles.defaultsfile,'plotbehaviormicroarray.m','.plotbehaviormicroarrayrc.mat');

function handles = initializeparams(handles,params)

handles.paramsvars = {'autochoose','plotsettings','lasttrxname','lastsegname'};

% intialize all parameters if there is nothing to read in
handles = hardcoded_initializeparams(handles);
if exist(handles.defaultsfile,'file'),
  tmp = load(handles.defaultsfile);
  fns = fieldnames(tmp);
  for i = 1:length(fns),
    handles.(fns{i}) = tmp.(fns{i});
  end
end

% set from input
for i = 1:length(handles.paramsvars),
  fn = handles.paramsvars{i};
  if isfield(params,fn),
    handles.(fn) = params.(fn);
  end
end

% set trxnames, behaviors, microarray
handles = initializetrxnames(handles);

function handles = initializetrxnames(handles)

% get at least one trx file
if isempty(handles.trxnames),
  handles = initializetrxnames_empty(handles);
else
  handles = checktrxandsegnames(handles);
  if isempty(handles.trxnames),
    handles = initializetrxnames_empty(handles);
  end
end

function handles = checktrxandsegnames(handles)

ntrx = length(handles.trxnames);
nbehaviors = length(handles.behaviors);

% check that all trx files exist
doesnotexist = [];
for i = 1:ntrx,
  if ~exist(handles.trxnames{i},'file'),
    warning('File %s does not exist, ignoring\n',handles.trxnames{i});
    doesnotexist(end+1) = i;
  end
end
% delete trx files that don't exist
handles.trxnames(doesnotexist) = [];
for i = 1:nbehaviors,
  removei = doesnotexist(doesnotexist <= length(handles.behaviors(i).segnames));
  handles.behaviors(i).segnames(removei) = [];
end

% there has to be one segfile for each trx file for each behavior
deleteidx = [];
for i = 1:nbehaviors,
  nseg = length(handles.behaviors(i).segnames);
  % delete extra seg files
  if nseg > ntrx,
    handles.behaviors(i).segnames = handles.behaviors(i).segnames(1:ntrx);
  elseif nseg < ntrx,
    % delete this behavior if not enough seg files
    warning('Not enough seg files for behavior %s, ignoring\n',handles.behaviors(i).name);
    deleteidx(end+1) = i;
  end
end
handles.behaviors(deleteidx) = [];

% check that all seg files exist
doesnotexist = [];
for i = 1:nbehaviors,
  for j = 1:length(handles.behaviors(i).segnames),
    if ~exist(handles.behaviors(i).segnames{j},'file'),
      warning('File %s does not exist, ignoring behavior %s\n',...
              handles.behaviors(i).segnames{j},handles.behaviors(i).name);
      doesnotexist(end+1) = i;
      break;
    end
  end
end
handles.behaviors(doesnotexist) = [];
nbehaviors = length(handles.behaviors);

% add each trx file
removetrxfile = [];
for i = 1:ntrx,
  trxname = handles.trxnames{i};
  segnames = cell(1,nbehaviors);
  for j = 1:nbehaviors,
    segnames{j} = handles.behaviors(j).segname{i};
  end
  [handles,succeeded] = addnewtrxfile(handles,trxname,segnames);
  if ~succeeded,
    warning('Could not add trx file %s, ignoring\n',trxname);
    removetrxfile(end+1) = i;
  end
end
handles.trxnames(removetrxfile) = [];
for i = 1:nbehaviors,
  handles.behaviors(i).segnames(removetrxfile) = [];
end

function handles = initializetrxnames_empty(handles)
% force no behaviors
handles.behaviors = struct('name',cell(1,0),'segnames',cell(1,0));
handles.behaviorselected = 1;

helpmsg = 'Choose mat file(s) containing trajectories with per-frame properties';
[matnames,matpath] = uigetfilehelp('*.mat','Choose per-frame properties file',handles.lasttrxname,'MultiSelect','on','helpmsg',helpmsg);

if isnumeric(matnames) && matnames == 0,
  % cancel
  if isempty(handles.trxnames),
    error('At least one per-frame properties file must be selected');
  else
    return;
  end
elseif ischar(matnames),
  matnames = {matnames};
end

for i = 1:length(matnames),
  matname = [matpath,matnames{i}];
  [handles,succeeded] = addnewtrxfile(handles,matname,{});
  if ~succeeded,
    continue;
  end
end

if isempty(handles.trxnames),
  error('No trajectories could be loaded from selected file(s).');
end

function [handles,succeeded] = addnewtrxfile(handles,matname,segnames)

succeeded = false;
if ~exist(matname,'file'),
  msgbox(sprintf('File %s does not exist\n',matname));
  return;
end

fprintf('Checking %s...\n',matname);
[trx,matname,loadsucceeded] = load_tracks(matname);

if ~loadsucceeded,
  fprintf('Could not load trajectories from file %s\n',matname);
  return;
end

if ~isperframe(trx),
  b = questdlg('Per-frame properties have not yet been computed for %s. Compute now, or skip this file?',...
    'Compute Per-Frame Properties?','Compute','Skip','Compute');
  if strcmpi(b,'skip'),
    return;
  end
  [matpath,matname0] = split_path_and_filename(matname);
  [computesucceeded,matname,trx] = compute_perframe_stats_f('matname',matname0,'matpath',matpath);
  if ~computesucceeded,
    return;
  end
end

% initialize microarray
[newmicroarray,isok] = makemicroarraypertrx(trx,segnames,getbasename(matname));
if ~isok,
  msgbox(sprintf('Error creating microarray for file %s\n',matname));
  return;
end
if isempty(handles.trxnames),
  handles.microarrays = newmicroarray;
  handles.nflies = length(newmicroarray);
else
  handles.nflies(end+1) = length(newmicroarray);
  
  % make the fields the same
  newfields = setdiff(fieldnames(newmicroarray),fieldnames(handles.microarrays));
  for i = 1:length(newfields),
    for j = 1:length(handles.microarrays),
      handles.microarrays(j).(newfields{i}) = nan(size(newmicroarray(1).(newfields{i})));
    end
  end
  newfields = setdiff(fieldnames(handles.microarrays),fieldnames(newmicroarray));
  for i = 1:length(newfields),
    for j = 1:length(newmicroarray),
      newmicroarray(j).(newfields{i}) = nan(size(handles.microarrays(1).(newfields{i})));
    end
  end   
  handles.microarrays = orderfields(handles.microarrays);
  newmicroarray = orderfields(newmicroarray);
  handles.microarrays(end+1:end+handles.nflies(end)) = newmicroarray;
end

handles.propnamespertrx{end+1} = getperframepropnames(trx);
if isempty(handles.trxnames),
  handles.allpropnames = handles.propnamespertrx{end};
else
  handles.allpropnames = intersect(handles.allpropnames,handles.propnamespertrx{end});
end

handles.trxnames{end+1} = matname;
handles.lasttrxname = matname;
succeeded = true;

if exist(handles.defaultsfile,'file'),
  save(handles.defaultsfile,'-append','-struct','handles','lasttrxname');
else
  save(handles.defaultsfile,'-struct','handles','lasttrxname');
end

function [handles,succeeded] = addnewbehavior(handles,segnames,behaviorname)

succeeded = false;
ntrx = length(handles.trxnames);
if length(segnames) ~= ntrx,
  msgbox(sprintf('Number of segfiles = %d does not match ntrx = %d\n',length(segnames),ntrx));
  return;
end

for i = 1:ntrx,
  if ~exist(segnames{i},'file')
    msgbox(sprintf('File %s does not exist\n',segnames{i}));
    return;
  end
end

for i = 1:ntrx,
  fprintf('Checking %s...\n',segnames{i});
  datacurr = load(segnames{i});
  isok = isfield(datacurr,'seg') && length(datacurr.seg) == handles.nflies(i);
  if ~isok,
    msgbox(sprintf('File %s does not contain seg with %d entries\n',segnames{i},handles.nflies(i)));
    return;
  end
end

% update microarray
nbehaviors = length(handles.behaviors);
off = 0;
for i = 1:ntrx,
  segnamescurr = cell(1,nbehaviors+1);
  for j = 1:nbehaviors,
    segnamescurr{j} = handles.behaviors(j).segnames{i};
  end
  segnamescurr{end} = segnames{i};
  idx = off + (1:handles.nflies(i));
  fprintf('Adding behavior "%s" for trx file %s\n',behaviorname,handles.trxnames{i});
  [newmicroarray,isok] = makemicroarraypertrx(handles.trxnames{i},...
                  segnamescurr,getbasename(handles.trxnames{i}),handles.microarrays(idx));
  if ~isok,
    msgbox(sprintf('Error creating microarray for file %s\n',handles.trxnames{i}));
    return;
  end
  % make the fields the same
  newfields = setdiff(fieldnames(newmicroarray),fieldnames(handles.microarrays));
  for k = 1:length(newfields),
    for j = 1:length(handles.microarrays),
      handles.microarrays(j).(newfields{k}) = nan(size(newmicroarray(1).(newfields{k})));
    end
  end
  newfields = setdiff(fieldnames(handles.microarrays),fieldnames(newmicroarray));
  for k = 1:length(newfields),
    for j = 1:length(newmicroarray),
      newmicroarray(j).(newfields{k}) = nan(size(handles.microarrays(1).(newfields{k})));
    end
  end
  handles.microarrays = orderfields(handles.microarrays);
  newmicroarray = orderfields(newmicroarray);
  handles.microarrays(idx) = newmicroarray;
  
  off = off + handles.nflies(i);
end

handles.behaviors(end+1).name = behaviorname;
handles.behaviors(end).segnames = segnames;
succeeded = true;

function handles = initializetrajectories(handles)

if isempty(handles.trxnames),
  handles = initializetrxnames_empty(handles);
end
set(handles.trxfilelist,'string',handles.trxnames);
set(handles.trxfilelist,'value',handles.trxselected);

function handles = initializebehaviors(handles)

if isempty(handles.behaviors),
  set(handles.behaviormenu,'string',{'Empty'},'value',1,'enable','off');
  set(handles.behaviornameedit,'string','','enable','off');
  set(handles.segfilelist,'string',{'Empty'},'value',1,'enable','off');
  set(handles.changesegfilebutton,'enable','off');
  set(handles.addbehaviorbutton,'enable','on');  
  set(handles.deletebehaviorbutton,'enable','off');
else
  behaviornames = {handles.behaviors.name};
  set(handles.behaviormenu,'string',behaviornames,'value',handles.behaviorselected,...
      'enable','on');
  set(handles.behaviornameedit,'string',handles.behaviors(handles.behaviorselected).name,...
                    'enable','on');
  set(handles.segfilelist,'string',handles.behaviors(handles.behaviorselected).segnames,...
                    'value',handles.trxselected,'enable','on');
  set(handles.changesegfilebutton,'enable','on');
  set(handles.addbehaviorbutton,'enable','on');  
  set(handles.deletebehaviorbutton,'enable','on');
end

function handles = initializeprops(handles)

if isempty(handles.props),
  set(handles.propertymenu,'string',{'Empty'},'value',1,'enable','off');
  set(handles.propnameedit,'string','','enable','off');
  set(handles.propbehaviormenu,'string',{handles.behaviors.name},'value',1,'enable','off');
  set(handles.durationbutton,'value',0,'enable','off');
  set(handles.fractimebutton,'value',0,'enable','off');
  set(handles.freqbutton,'value',0,'enable','off');
  set(handles.perframebutton,'value',1,'enable','off');
  set(handles.perframepanel,'visible','off');
  set(handles.addpropbutton,'enable','on');
  set(handles.deletepropbutton,'enable','off');
else
  prop = handles.props(handles.propselected);
  nobehaviors = isempty(handles.behaviors);
  set(handles.propertymenu,'string',{handles.props.name},'enable','on');
  set(handles.propertymenu,'value',handles.propselected);
  set(handles.propnameedit,'string',prop.name,'enable','on');
  if nobehaviors,
    set(handles.propbehaviormenu,'string',{'Empty'});
    prop.proptype = 'perframe';
    prop.during = 'allframes';
  else
    set(handles.propbehaviormenu,'string',{handles.behaviors.name},'enable','on');
  end
  % behavior is irrelevant
  if strcmpi(prop.proptype,'perframe') && ...
        strcmpi(prop.during,'allframes'),
    set(handles.propbehaviormenu,'enable','off');
  else
    set(handles.propbehaviormenu,'enable','on','value',prop.behavior);
  end
  
  if nobehaviors,
    set(handles.durationbutton,'value',0,'enable','off');
    set(handles.fractimebutton,'value',0,'enable','off');
    set(handles.freqbutton,'value',0,'enable','off');
    set(handles.perframebutton,'value',1,'enable','on');
  else
    set(handles.durationbutton,'enable','on','value',strcmpi(prop.proptype,'duration'));
    set(handles.fractimebutton,'enable','on','value',strcmpi(prop.proptype,'fractime'));
    set(handles.freqbutton,'enable','on','value',strcmpi(prop.proptype,'frequency'));
    set(handles.perframebutton,'enable','on','value',strcmpi(prop.proptype,'perframe'));
  end
  set(handles.perframemenu,'string',handles.allpropnames);
  i = find(strcmpi(handles.allpropnames,prop.perframe));
  if isempty(i),
    warning('Could not find property %s on allpropnames\n',prop.perframe);
    i = 1;
    prop.perframe = handles.allpropnames{1};
  end
  set(handles.perframemenu,'value',i);
  set(handles.isbehaviorbutton,'value',strcmpi(prop.during,'during'));
  set(handles.isnotbehaviorbutton,'value',strcmpi(prop.during,'invert'));
  set(handles.allframesbutton,'value',strcmpi(prop.during,'allframes'));
  if strcmpi(prop.during,'allframes'),
    set(handles.averagingmenu,'enable','off');
    prop.averaging = 'none (per-frame)';
  else
    set(handles.averagingmenu,'enable','on');
  end
  i = find(strcmpi(handles.allaveraging,prop.averaging));
  if isempty(i),
    warning('Could not find averaging %s on allaveraging\n',prop.averaging);
    i = 1;
    prop.averaging = handles.allaveraging{i};
  end
  set(handles.averagingmenu,'string',handles.allaveraging,'value',i);
  if strcmpi(prop.proptype,'perframe'),
    set(handles.perframepanel,'visible','on');
  else
    set(handles.perframepanel,'visible','off');
  end  
  
  if nobehaviors
    set(handles.isbehaviorbutton,'enable','off');
    set(handles.isnotbehaviorbutton,'enable','off');
  else
    set(handles.isbehaviorbutton,'enable','on');
    set(handles.isnotbehaviorbutton,'enable','on');
  end
  
  % save changes
  handles.props(handles.propselected) = prop;
end

function nautochoose = computenautochoose(handles)

% compute number of properties to choose from
% nbehaviors * [isduration + isfractime + isfrequency + 
%   isperframe * (isduring+isinvert) * (nprops * (isnoaveraging+ismean+ismedian+isstart+isend))]
% + isperframe * isallframes * nprops
nbehaviors = length(handles.behaviors);
nprops = length(handles.allpropnames);
isduration = double(strcmpi(handles.autochoose.duration,'allow'));
isfractime = double(strcmpi(handles.autochoose.fractime,'allow'));
isfrequency = double(strcmpi(handles.autochoose.freq,'allow'));
isperframe = double(strcmpi(handles.autochoose.perframe,'allow'));
isduring = double(strcmpi(handles.autochoose.duringbehavior,'allow'));
isinvert = double(strcmpi(handles.autochoose.invertbehavior,'allow'));
isnoaveraging = double(strcmpi(handles.autochoose.noaveraging,'allow'));
ismean = double(strcmpi(handles.autochoose.mean,'allow'));
ismedian = double(strcmpi(handles.autochoose.median,'allow'));
isstart = double(strcmpi(handles.autochoose.start,'allow'));
isend = double(strcmpi(handles.autochoose.end,'allow'));
isallframes = double(strcmpi(handles.autochoose.allframes,'allow'));
nignoreprops = nnz(ismember(handles.autochoose.ignoreprops,handles.allpropnames));
nautochoose = ...
  nbehaviors*(isduration + isfractime + isfrequency + ...
              isperframe*(isduring+isinvert)*...
              ((nprops-nignoreprops)*(isnoaveraging+ismean+ismedian+isstart+isend))) + ...
    isperframe*isallframes*(nprops-nignoreprops);


function handles = initializeautochoose(handles);

handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));

for i = 1:length(handles.autochooseallownames),
  s = handles.autochooseallownames{i};
  if strcmpi(handles.autochoose.(s),'allow'),
    v = 1;
  else
    v = 2;
  end
  set(handles.(sprintf('allow%smenu',s)),'value',v);
end
ignorelist = [{'None'},handles.allpropnames];
set(handles.ignorepropertieslist,'string',ignorelist,...
                  'min',0,'max',length(ignorelist));
v = [];
for i = 1:length(handles.autochoose.ignoreprops),
  j = find(strcmpi(handles.autochoose.ignoreprops{i},handles.allpropnames));
  v = [v,j];
end
set(handles.ignorepropertieslist,'value',v+1);
set(handles.npropsautoshowedit,'string',num2str(handles.autochoose.nshow));
set(handles.chooseindepbox,'value',handles.autochoose.isindep);

function handles = initializegui(handles)

if isempty(handles.trxnames),
  return;
end

% trajectories
handles = initializetrajectories(handles);

% behaviors
handles = initializebehaviors(handles);

% properties
handles = initializeprops(handles);

% auto choose
handles = initializeautochoose(handles);

% types
handles = initializetypes(handles);

% plot
handles = initializeplotpanel(handles);

function handles = initializeplotpanel(handles)

contents = get(handles.errorbarsmenu,'string');
i = find(strcmpi(contents,handles.plotsettings.errorbars));
if isempty(i),
  warning('could not find errorbars = %s on errorbarsmenu',plotsettings.errorbars);
  i = 1;
  handles.plotsettings.errorbars = contents{i};
end
set(handles.errorbarsmenu,'value',i);
set(handles.plotindivsbox,'value',handles.plotsettings.plotindivs);
set(handles.normalizebox,'value',handles.plotsettings.normalize);

function handles = initializetypes(handles)

if isempty(handles.types),
  set(handles.typemenu,'string',{'Empty'},'enable','off','value',1);
  set(handles.typenameedit,'string','','enable','off');
  set(handles.typeslist,'string',handles.alltypes,'enable','off');
  %set(handles.edittypesbutton,'enable','off');
  set(handles.deletetypebutton,'enable','off');
  set(handles.addnewtypebutton,'enable','on');
else
  set(handles.typemenu,'enable','on','value',1);
  set(handles.typenameedit,'enable','on');
  set(handles.typeslist,'enable','on');
  %set(handles.edittypesbutton,'enable','on');
  set(handles.deletetypebutton,'enable','on');
  set(handles.addnewtypebutton,'enable','on');
  typecurr = handles.types(handles.typeselected);
  handles.alltypes = unique({handles.microarrays.type});
  handles.alltypes = [{'All'},handles.alltypes];
  set(handles.typemenu,'string',{handles.types.name});
  set(handles.typemenu,'value',handles.typeselected);
  set(handles.typenameedit,'string',typecurr.name);
  set(handles.typeslist,'string',handles.alltypes);
  v = [];
  for i = 1:length(typecurr.fields),
    j = find(strcmpi(typecurr.fields(i),handles.alltypes));
    v = [v,j];
  end
  if isempty(v),
    warning('Could not find any types in alltype');
    v = 1;
    typecurr.fields = handles.alltypes(v);
    handles.types(handles.typeselected) = typecurr;
  end
  set(handles.typeslist,'value',v);
end

% --- Outputs from this function are returned to the command line.
function varargout = plotbehaviormicroarray_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in perframemenu.
function perframemenu_Callback(hObject, eventdata, handles)
% hObject    handle to perframemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns perframemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from perframemenu
contents = get(hObject,'String');
handles.props(handles.propselected).perframe = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function perframemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perframemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in propbehaviormenu.
function propbehaviormenu_Callback(hObject, eventdata, handles)
% hObject    handle to propbehaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns propbehaviormenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from propbehaviormenu
handles.props(handles.propselected).behavior = get(hObject,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function propbehaviormenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propbehaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trxfilelist.
function trxfilelist_Callback(hObject, eventdata, handles)
% hObject    handle to trxfilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns trxfilelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trxfilelist
handles.trxselected = get(hObject,'value');
if ~isempty(handles.behaviors),
  set(handles.segfilelist,'value',handles.trxselected);
end

% --- Executes during object creation, after setting all properties.
function trxfilelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trxfilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addtrxbutton.
function addtrxbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addtrxbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get name of trx file
helpmsg = {};
helpmsg{1} = 'Choose mat file(s) containing trajectories with per-frame properties. Mat files loaded in so far:';
helpmsg(2:1+length(handles.trxnames)) = handles.trxnames;
while true,
  [matname,matpath] = uigetfilehelp('*.mat','Choose per-frame properties file',handles.lasttrxname,'helpmsg',helpmsg);
  if isnumeric(matname) && matname == 0,
    return;
  end
  trxname = [matpath,matname];
  if ~exist(trxname,'file'),
    msgbox(sprintf('Trx file %s does not exist',trxname));
    continue;
  end
  break;
end
handles.lasttrxname = trxname;

% get names of seg files corresponding to this trx file
nbehaviors = length(handles.behaviors);
segnames = cell(1,nbehaviors);
for i = 1:nbehaviors,
  s = sprintf('Choose seg file for behavior %s and trx file %s',...
              handles.behaviors(i).name,getbasename(trxname));
  helpmsg = sprintf('Choose seg file for behavior %s and trx file %s',...
    handles.behaviors(i).name,trxname);
  while true
    [matname,matpath] = uigetfilehelp('*.mat',s,handles.lastsegname,'helpmsg',helpmsg);
    if isnumeric(matname) && matname == 0,
      return;
    end
    segnames{i} = [matpath,matname];
    if ~exist(segnames{i},'file'),
      msgbox(sprintf('Seg file %s does not exist',segnames{i}));
      continue;
    end
    break;
  end
  handles.lastsegname = segnames{i};
end
  
% compute the microarray  
[handles,succeeded] = addnewtrxfile(handles,trxname,segnames);
if ~succeeded,
  return;
end

% save segnames to behaviors
for i = 1:nbehaviors,
  handles.behaviors(i).segnames{end+1} = segnames{i};
end

% property names may have changed
handles = updateallpropnames(handles);

% update gui
handles.trxselected = length(handles.trxnames);
handles = initializetrajectories(handles);
handles = initializebehaviors(handles);
handles = initializetypes(handles);

guidata(hObject,handles);

% --- Executes on button press in deletetrxbutton.
function deletetrxbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletetrxbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% can't delete last trxfile
if length(handles.trxnames) == 1,
  msgbox({'Cannot delete last trx file.',...
          'If you want to delete this file, add another file then delete it'});
  return;
end

idx = handles.trxselected;

% remove from trxnames
handles.trxnames(idx) = [];

% remove from microarray
off0 = sum(handles.nflies(1:idx-1));
handles.microarrays(off0+(1:handles.nflies(idx))) = [];

% remove from nflies
handles.nflies(idx) = [];

% remove from behaviors.segnames
for i = 1:length(handles.behaviors),
  handles.behaviors(i).segnames(idx) = [];
end

% remove from propnamespertrx
handles.propnamespertrx(idx) = [];

% may need to change trxselected
if handles.trxselected > length(handles.trxnames),
  handles.trxselected = length(handles.trxnames);
end

% may need to reset allpropnames
newallpropnames = handles.propnamespertrx{1};
for i = 2:length(handles.propnamespertrx),
  newallpropnames = intersect(newallpropnames,handles.propnamespertrx{i});
end
ischange = length(intersect(newallpropnames,handles.allpropnames)) < length(newallpropnames);
if ischange,
  handles.allpropnames = newallpropnames;
  handles = updateallpropnames(handles);
end

% update trajectories
handles = initializetrajectories(handles);
% update behaviors
handles = initializebehaviors(handles);
% update types
handles = initializetypes(handles);

% save
guidata(hObject,handles);

function handles = updateallpropnames(handles)

handles = update_perframemenu(handles);
handles = update_ignoreproperties(handles);

function handles = update_ignoreproperties(handles)

handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
ignorelist = [{'None'},handles.allpropnames];
set(handles.ignorepropertieslist,'string',ignorelist,...
                  'min',0,'max',length(ignorelist));
set(handles.ignorepropertieslist,'string',ignorelist,...
                  'min',0,'max',length(handles.allpropnames));
v = [];
for i = 1:length(handles.autochoose.ignoreprops),
  j = find(strcmpi(handles.autochoose.ignoreprops{i},handles.allpropnames));
  v = [v,j];
end
set(handles.ignorepropertieslist,'value',v+1);

function handles = update_perframemenu(handles)

if isempty(handles.props),
  return;
end
prop = handles.props(handles.propselected);
set(handles.perframemenu,'string',handles.allpropnames);
i = find(strcmpi(handles.allpropnames,prop.perframe));
if isempty(i),
  warning('Could not find property %s on allpropnames\n',prop.perframe);
  i = 1;
  prop.perframe = handles.allpropnames{1};
end
set(handles.perframemenu,'value',i);
handles.props(handles.propselected) = prop;

% --- Executes on button press in addpropbutton.
function addpropbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addpropbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propnames = {handles.props.name};
nprops = length(handles.props);
for i = 1:inf,
  prop.name = sprintf('Property %d',i);
  if ~ismember(prop.name,propnames),
    break;
  end
end
prop.behavior = 1;
prop.proptype = 'perframe';
prop.perframe = 'velmag';
prop.during = 'allframes';
prop.averaging = 'None (per-frame)';
handles.props(nprops+1) = prop;
handles.propselected = nprops+1;

handles = initializeprops(handles);

guidata(hObject,handles);

% --- Executes on button press in deletepropbutton.
function deletepropbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletepropbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(handles.props) == 1,
  msgbox('Cannot delete only property');
  return;
end
idx = handles.propselected;
handles.props(idx) = [];
if idx > length(handles.props),
  handles.propselected = handles.propselected - 1;
end
handles = initializeprops(handles);

guidata(hObject,handles);

% --- Executes on selection change in propertymenu.
function propertymenu_Callback(hObject, eventdata, handles)
% hObject    handle to propertymenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns propertymenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from propertymenu

handles.propselected = get(hObject,'value');
handles = initializeprops(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function propertymenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propertymenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in behaviormenu.
function behaviormenu_Callback(hObject, eventdata, handles)
% hObject    handle to behaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns behaviormenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from behaviormenu

handles.behaviorselected = get(hObject,'value');
handles = initializebehaviors(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function behaviormenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to behaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segfilelist.
function segfilelist_Callback(hObject, eventdata, handles)
% hObject    handle to segfilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns segfilelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segfilelist
segi = get(hObject,'value');
handles.trxselected = segi;

% trxfilelist and segfilelist are linked
set(handles.trxfilelist,'value',segi);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function segfilelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segfilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in changesegfilebutton.
function changesegfilebutton_Callback(hObject, eventdata, handles)
% hObject    handle to changesegfilebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change seg file for trxselected

ntrx = length(handles.trxnames);
nbehaviors = length(handles.behaviors);
segi = min(length(handles.behaviors(1).segnames),handles.trxselected);
behi = min(length(handles.behaviors),handles.behaviorselected);

% get names of seg files corresponding to this behavior file
segnames = cell(1,nbehaviors);
for i = 1:nbehaviors,
  segnames{i} = handles.behaviors(i).segnames{segi};
end

s = sprintf('Choose seg file for trx file %s and behavior "%s"',...
            getbasename(handles.trxnames{segi}),...
            handles.behaviors(behi).name);
helpmsg = sprintf('Choose seg file for trx file %s and behavior "%s"',...
  handles.trxnames{segi},...
  handles.behaviors(behi).name);
while true
  fprintf([s,'\n']);
  [matname,matpath] = uigetfilehelp('*.mat',s,segnames{segi},'helpmsg',helpmsg);
  if isnumeric(matname) && matname == 0,
    return;
  end
  segname = [matpath,matname];
  if ~exist(segname,'file'),
    msgbox(sprintf('Seg file %s does not exist',segname));
    continue;
  end
  break;
end
segnames{segi} = segname;
handles.lastsegname = segnames{segi};
  
if exist(handles.defaultsfile,'file'),
  save(handles.defaultsfile,'-append','-struct','handles','lastsegname');
else
  save(handles.defaultsfile,'-struct','handles','lastsegname');
end

% compute the microarray  
off = sum(handles.nflies(1:segi-1));
idx = off + (1:handles.nflies(segi));
fprintf('Editing behavior "%s" for trx file %s\n',handles.behaviors(behi).name,...
        handles.trxnames{segi});
[handles.microarrays(idx),succeeded] = ...
    makemicroarraypertrx(handles.trxnames{segi},...
                         segnames,...
                         getbasename(handles.trxnames{segi}),...
                         handles.microarrays(idx),...
                         behi);
if ~succeeded,
  return;
end

handles.behaviors(behi).segnames{segi} = segnames{behi};
% update the gui
handles = initializebehaviors(handles);

guidata(hObject,handles);

% --- Executes on button press in addbehaviorbutton.
function addbehaviorbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addbehaviorbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ntrx = length(handles.trxnames);
nbehaviors = length(handles.behaviors);

% get an unused name for this behavior
allnames = {handles.behaviors.name};
for i = 1:inf,
  name = sprintf('Behavior %d',i);
  if ~ismember(name,allnames),
    break;
  end
end

% get names of seg files corresponding to this behavior file
segnames = cell(1,ntrx);
for i = 1:ntrx,
  s = sprintf('Choose seg file for trx file %s',getbasename(handles.trxnames{i}));
  helpmsg = sprintf('Choose seg file for trx file %s',handles.trxnames{i});
  while true
    fprintf([s,'\n']);
    [matname,matpath] = uigetfilehelp('*.mat',s,handles.lastsegname,'helpmsg',helpmsg);
    if isnumeric(matname) && matname == 0,
      return;
    end
    segnames{i} = [matpath,matname];
    if ~exist(segnames{i},'file'),
      msgbox(sprintf('Seg file %s does not exist',segnames{i}));
      continue;
    end
    break;
  end
  handles.lastsegname = segnames{i};
end
  
if exist(handles.defaultsfile,'file'),
  save(handles.defaultsfile,'-append','-struct','handles','lastsegname');
else
  save(handles.defaultsfile,'-struct','handles','lastsegname');
end

% compute the microarray  
[handles,succeeded] = addnewbehavior(handles,segnames,name);
if ~succeeded,
  return;
end

% update the gui
nbehaviors = nbehaviors + 1;
handles.behaviorselected = nbehaviors;
handles = initializebehaviors(handles);
handles = initializeprops(handles);
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));

guidata(hObject,handles);


% --- Executes on button press in deletebehaviorbutton.
function deletebehaviorbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletebehaviorbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

behi = handles.behaviorselected;
name = handles.behaviors(behi).name;

% remove from behaviors
handles.behaviors(behi) = [];

% delete properties requiring this behavior
dodelete = [];
for i = 1:length(handles.props),
  prop = handles.props(i);
  % doesn't depend on behavior or depends on a different behavior
  if (strcmpi(prop.proptype,'perframe') && strcmpi(prop.during,'allframes')) || ...
        prop.behavior ~= behi,
    % decrement behavior so that things line up
    if prop.behavior >= behi,
      prop.behavior = max(1,prop.behavior - 1);
      handles.props(i) = prop;
    end
  else
    % delete this behavior
    warning('Deleting property %s as it relies on behavior %s\n',prop.name,name);
    dodelete(end+1) = i;  
  end
end
handles.props(dodelete) = [];

% remove from microarrays
durings = {'during','invert'};
averagings = {'none','mean','median','start','end'};
nflies = length(handles.microarrays);
off = 0;
trxi = 0;
for fly = 1:nflies,
  if fly > off,
    trxi = trxi + 1;
    off = off + handles.nflies(trxi);
    % need names of all propetries for this trx file
    props = handles.propnamespertrx{trxi};
    nprops = length(props);
  end
  handles.microarrays(fly).duration(behi) = [];
  handles.microarrays(fly).fractime(behi) = [];
  handles.microarrays(fly).frequency(behi) = [];
  
  for propi = 1:nprops,
    prop = props{propi};
    for duringi = 1:length(durings),
      during = durings{duringi};
      for avei = 1:length(averagings),
        averaging = averagings{avei};
        s = sprintf('%s_%s_%s',prop,during,averaging);
        handles.microarrays(fly).(s)(behi) = [];
      end
    end
  end
end

% update gui

% behavior box
if handles.behaviorselected > length(handles.behaviors) && ...
      handles.behaviorselected > 1,
  handles.behaviorselected = handles.behaviorselected - 1;
end
handles = initializebehaviors(handles);

% properties
handles = initializeprops(handles);

% update nautochoose
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));

guidata(hObject,handles);

function npropsautoshowedit_Callback(hObject, eventdata, handles)
% hObject    handle to npropsautoshowedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of npropsautoshowedit as text
%        str2double(get(hObject,'String')) returns contents of npropsautoshowedit as a double
v = str2double(get(hObject,'String'));
if isnan(v) || v <= 0,
  set(hObject,'string',num2str(handles.autochoose.nshow));
  return;
end
if v > handles.maxnautochoose,
  set(hObject,'string',num2str(handles.maxnautochoose));
  v = handles.maxnautochoose;
end
handles.autochoose.nshow = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function npropsautoshowedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to npropsautoshowedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autochoosebutton.
function autochoosebutton_Callback(hObject, eventdata, handles)
% hObject    handle to autochoosebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nallowed = handles.maxnautochoose;
nchoose = min(nallowed,handles.autochoose.nshow);
nflies = length(handles.microarrays);
nbehaviors = length(handles.behaviors);
nprops = length(handles.allpropnames);

isduration = nbehaviors > 0 && strcmpi(handles.autochoose.duration,'allow');
isfractime = nbehaviors > 0 && strcmpi(handles.autochoose.fractime,'allow');
isfrequency = nbehaviors > 0 && strcmpi(handles.autochoose.freq,'allow');
isperframe = strcmpi(handles.autochoose.perframe,'allow');
isduring = nbehaviors > 0 && strcmpi(handles.autochoose.duringbehavior,'allow');
isinvert = nbehaviors > 0 && strcmpi(handles.autochoose.invertbehavior,'allow');
isnoaveraging = nbehaviors > 0 && strcmpi(handles.autochoose.noaveraging,'allow');
ismean = nbehaviors > 0 && strcmpi(handles.autochoose.mean,'allow');
ismedian = nbehaviors > 0 && strcmpi(handles.autochoose.median,'allow');
isstart = nbehaviors > 0 && strcmpi(handles.autochoose.start,'allow');
isend = nbehaviors > 0 && strcmpi(handles.autochoose.end,'allow');
isallframes = strcmpi(handles.autochoose.allframes,'allow');

% grab out allowed fields
X = zeros(nflies,nallowed);
allowedprops.name = 'Property 1';
allowedprops.behavior = 1;
allowedprops.proptype = 'perframe';
allowedprops.perframe = 'velmag';
allowedprops.during = 'allframes';
allowedprops.averaging = 'None (per-frame)';
allowedprops = repmat(allowedprops,[1,nallowed]);

off = 0;
tmpbehaviors = mat2cell(1:nbehaviors,1,ones(1,nbehaviors));
behaviornames = {handles.behaviors.name};
if isduration,
  idx = off+(1:nbehaviors);
  X(:,idx) = reshape([handles.microarrays.duration],[nbehaviors,nflies])';
  [allowedprops(idx).proptype] = 'duration';
  [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
  for i = 1:nbehaviors,
    allowedprops(idx(i)).name = sprintf('%s duration',behaviornames{i});
  end
  off = off + nbehaviors;
end
if isfractime,
  idx = off+(1:nbehaviors);
  X(:,idx) = reshape([handles.microarrays.fractime],[nbehaviors,nflies])';
  [allowedprops(idx).proptype] = 'fractime';
  [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
  for i = 1:nbehaviors,
    allowedprops(idx(i)).name = sprintf('%s fractime',behaviornames{i});
  end
  off = off + nbehaviors;
end
if isfrequency,
  idx = off+(1:nbehaviors);
  X(:,idx) = reshape([handles.microarrays.frequency],[nbehaviors,nflies])';
  [allowedprops(idx).proptype] = 'frequency';
  [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
  for i = 1:nbehaviors,
    allowedprops(idx(i)).name = sprintf('%s frequency',behaviornames{i});
  end
  off = off + nbehaviors;
end
if isperframe,
  propsallowed = setdiff(handles.allpropnames,handles.autochoose.ignoreprops);
  if isduring,
    if isnoaveraging,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_during_none',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('during');
        [allowedprops(idx).averaging] = deal('None (per-frame)');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        off = off + nbehaviors;
      end      
    end
    if ismean,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_during_mean',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('during');
        [allowedprops(idx).averaging] = deal('interval mean');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        off = off + nbehaviors;
      end      
    end
    if ismedian,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_during_median',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('during');
        [allowedprops(idx).averaging] = deal('interval median');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end
    if isstart,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_during_start',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('during');
        [allowedprops(idx).averaging] = deal('interval start');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end
    if isend,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_during_end',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('during');
        [allowedprops(idx).averaging] = deal('interval end');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end

  end % end isduring
  
  if isinvert,
    if isnoaveraging,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_invert_none',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('invert');
        [allowedprops(idx).averaging] = deal('None (per-frame)');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end
    if ismean,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_invert_mean',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('invert');
        [allowedprops(idx).averaging] = deal('interval mean');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end
    if ismedian,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_invert_median',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('invert');
        [allowedprops(idx).averaging] = deal('interval median');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end
    if isstart,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_invert_start',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('invert');
        [allowedprops(idx).averaging] = deal('interval start');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end
    if isend,
      for propi = 1:length(propsallowed),
        prop = propsallowed{propi};
        s = sprintf('%s_invert_end',prop);
        idx = off+(1:nbehaviors);
        X(:,idx) = reshape([handles.microarrays.(s)],[nbehaviors,nflies])';
        [allowedprops(idx).proptype] = deal('perframe');
        [allowedprops(idx).behavior] = deal(tmpbehaviors{:});
        [allowedprops(idx).perframe] = deal(prop);
        [allowedprops(idx).during] = deal('invert');
        [allowedprops(idx).averaging] = deal('interval end');
        for i = 1:nbehaviors,
          allowedprops(idx(i)).name = sprintf('%s %s',behaviornames{i},s);
        end
        %[fieldsallowed{idx}] = deal(s);
        %fieldsallowedbehi(idx) = 1:nbehaviors;
        off = off + nbehaviors;
      end      
    end

  end % end isinvert
  
  if isallframes,
    
    for propi = 1:length(propsallowed),
      prop = propsallowed{propi};
      s = sprintf('%s_allframes',prop);
      idx = off+1;
      X(:,idx) = [handles.microarrays.(s)]';
      allowedprops(idx).proptype = 'perframe';
      allowedprops(idx).perframe = prop;
      allowedprops(idx).during = 'allframes';
      allowedprops(idx).name = s;
      %[fieldsallowed{idx}] = deal(s);
      %fieldsallowedbehi(idx) = nan;
      off = off + 1;
    end
    
  end % end isallframes
  
end % end isperframe

% get types 
ntypes = length(handles.types);
isall = false(1,ntypes);
types = false(nflies,ntypes);
for i = 1:ntypes,
  % if all is one of the types, then include all flies
  if any(strcmpi('All',handles.types(i).fields)),
    types(:,i) = true;
    isall(i) = true;
  else
    for fly = 1:nflies,
      if any(strcmpi(handles.microarrays(fly).type,handles.types(i).fields)),
        types(fly,i) = true;
      end
    end
  end
end
y = double(types);
y = y ./ repmat(sum(y,2),[1,ntypes]);

% remove flies with no type
ignoreflies = any(isnan(X),2) | any(isnan(y),2);
  
chosenprops = choosediscriminativeprops(X(~ignoreflies,:),y(~ignoreflies,:),nchoose,...
                                        'isindependent',handles.autochoose.isindep);
nchosen = length(chosenprops);

% create props structure
handles.props = allowedprops(chosenprops);
handles = initializeprops(handles);
guidata(hObject,handles);

% --- Executes on selection change in allowdurationmenu.
function allowdurationmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowdurationmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowdurationmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowdurationmenu

contents = get(hObject,'String');
handles.autochoose.duration = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function allowdurationmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowdurationmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowfractimemenu.
function allowfractimemenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowfractimemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowfractimemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowfractimemenu
contents = get(hObject,'String');
handles.autochoose.fractime = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function allowfractimemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowfractimemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowfreqmenu.
function allowfreqmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowfreqmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowfreqmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowfreqmenu
contents = get(hObject,'String');
handles.autochoose.freq = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowfreqmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowfreqmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowperframemenu.
function allowperframemenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowperframemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowperframemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowperframemenu
contents = get(hObject,'String');
handles.autochoose.perframe = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowperframemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowperframemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ignorepropertieslist.
function ignorepropertieslist_Callback(hObject, eventdata, handles)
% hObject    handle to ignorepropertieslist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ignorepropertieslist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ignorepropertieslist
contents = get(hObject,'String');
handles.autochoose.ignoreprops = contents(get(hObject,'Value'));
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ignorepropertieslist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ignorepropertieslist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allownoaveragingmenu.
function allownoaveragingmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allownoaveragingmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allownoaveragingmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allownoaveragingmenu

contents = get(hObject,'String');
handles.autochoose.noaveraging = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function allownoaveragingmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allownoaveragingmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowmeanmenu.
function allowmeanmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowmeanmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowmeanmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowmeanmenu
contents = get(hObject,'String');
handles.autochoose.mean = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowmeanmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowmeanmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowstartmenu.
function allowstartmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowstartmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowstartmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowstartmenu
contents = get(hObject,'String');
handles.autochoose.start = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowstartmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowstartmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowendmenu.
function allowendmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowendmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowendmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowendmenu
contents = get(hObject,'String');
handles.autochoose.end = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowendmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowendmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowduringbehaviormenu.
function allowduringbehaviormenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowduringbehaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowduringbehaviormenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowduringbehaviormenu
contents = get(hObject,'String');
handles.autochoose.duringbehavior = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowduringbehaviormenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowduringbehaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowinvertbehaviormenu.
function allowinvertbehaviormenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowinvertbehaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowinvertbehaviormenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowinvertbehaviormenu
contents = get(hObject,'String');
handles.autochoose.invertbehavior = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowinvertbehaviormenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowinvertbehaviormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowallframesmenu.
function allowallframesmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowallframesmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowallframesmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowallframesmenu
contents = get(hObject,'String');
handles.autochoose.allframes = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowallframesmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowallframesmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allowmedianmenu.
function allowmedianmenu_Callback(hObject, eventdata, handles)
% hObject    handle to allowmedianmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns allowmedianmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allowmedianmenu
contents = get(hObject,'String');
handles.autochoose.median = contents{get(hObject,'Value')};
handles.maxnautochoose = computenautochoose(handles);
set(handles.npropertiestext,'string',num2str(handles.maxnautochoose));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function allowmedianmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowmedianmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in errorbarsmenu.
function errorbarsmenu_Callback(hObject, eventdata, handles)
% hObject    handle to errorbarsmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns errorbarsmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from errorbarsmenu
contents = get(hObject,'String');
handles.plotsettings.errorbars = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function errorbarsmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to errorbarsmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updatebutton.
function updatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% plot the microarray
nprops = length(handles.props);
ntypes = length(handles.types);
nflies = length(handles.microarrays);
propsperfly = zeros(nflies,nprops);
types = false(nflies,ntypes);

% grab properties from microarrays
for i = 1:nprops,
  prop = handles.props(i);
  if ~strcmpi(prop.proptype,'perframe'),
    % not perframe
    s = prop.proptype;
    for fly = 1:nflies,
      propsperfly(fly,i) = handles.microarrays(fly).(s)(prop.behavior);
    end
  else
    % per-frame property
    if strcmpi(prop.during,'allframes'),
      % per-frame and allframes
      s = sprintf('%s_allframes',prop.perframe);
      for fly = 1:nflies,
        propsperfly(fly,i) = handles.microarrays(fly).(s);
      end
    else % per-frame, but not allframes
      if strcmpi(prop.averaging,'None (per-frame)'),
        prop.averaging = 'none';
      elseif strcmpi(prop.averaging,'interval mean'),
        prop.averaging = 'mean';
      elseif strcmpi(prop.averaging,'interval median'),
        prop.averaging = 'median';
      elseif strcmpi(prop.averaging,'interval start'),
        prop.averaging = 'start';
      elseif strcmpi(prop.averaging,'interval end'),
        prop.averaging = 'end';
      end
      s = sprintf('%s_%s_%s',prop.perframe,prop.during,prop.averaging);
      for fly = 1:nflies,
        propsperfly(fly,i) = handles.microarrays(fly).(s)(prop.behavior);
      end
    end % end per-frame, but not allframes
  end % end per-frame
end

% types
for i = 1:ntypes,
  % if all is one of the types, then include all flies
  if any(strcmpi('All',handles.types(i).fields)),
    types(:,i) = true;
  else
    for fly = 1:nflies,
      if any(strcmpi(handles.microarrays(fly).type,handles.types(i).fields)),
        types(fly,i) = true;
      end
    end
  end
end

propnames = {handles.props.name};
typenames = {handles.types.name};
if isfield(handles,'plotstuff'),
  hfig = handles.plotstuff.hfig;
  hfig2 = handles.plotstuff.hfig2;
else
  hfig = nan;
  hfig2 = nan;
end

handles.plotstuff = plotbehaviormicroarray_fcn(propsperfly,types,...
    'propnames',propnames,'typenames',typenames,...
    'plotindivs',handles.plotsettings.plotindivs,...
    'normalize',handles.plotsettings.normalize,...
    'errorbars',handles.plotsettings.errorbars,...
    'hfig',hfig,'hfig2',hfig2);

guidata(hObject,handles);

% --- Executes on button press in exportbutton.
function exportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'savename'),
  helpmsg = 'Choose file to save microarray statistics to.';
  [paths,names] = split_path_and_filename(handles.trxnames{end});
  [savename,savepath] = uiputfilehelp('*.mat','Save microarray as',paths,'helpmsg',helpmsg);
  if ~ischar(savename),
    return;
  end
  handles.savename = [savepath,savename];
  guidata(hObject,handles);
end
fprintf('Exporting current plot properties to file %s\n',handles.savename);
save(handles.savename,'-struct','handles','microarrays','trxnames','propnamespertrx',...
  'allpropnames','behaviors','props','autochoose','types','plotsettings','nflies',...
  'maxnautochoose','alltypes');
plotstuff = handles.plotstuff;
save('-append',handles.savename,'-struct','plotstuff','mu_all','sig_all','nflies_all','stderr_all','stderr');


% --- Executes on button press in plotindivsbox.
function plotindivsbox_Callback(hObject, eventdata, handles)
% hObject    handle to plotindivsbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotindivsbox
handles.plotsettings.plotindivs = get(hObject,'value');
guidata(hObject,handles);

% --- Executes on button press in normalizebox.
function normalizebox_Callback(hObject, eventdata, handles)
% hObject    handle to normalizebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalizebox
handles.plotsettings.normalize = get(hObject,'value');
guidata(hObject,handles);

function propnameedit_Callback(hObject, eventdata, handles)
% hObject    handle to propnameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of propnameedit as text
%        get(hObject,'String')) returns contents of propnameedit as a double
handles.props(handles.propselected).name = get(hObject,'String');
set(handles.propertymenu,'string',{handles.props.name});
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function propnameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propnameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in typemenu.
function typemenu_Callback(hObject, eventdata, handles)
% hObject    handle to typemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns typemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from typemenu
handles.typeselected = get(hObject,'value');
typecurr = handles.types(handles.typeselected);
v = [];
for i = 1:length(typecurr.fields),
  j = find(strcmpi(typecurr.fields(i),handles.alltypes));
  v = [v,j];
end
set(handles.typeslist,'value',v);
set(handles.typenameedit,'string',typecurr.name);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function typemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function typenameedit_Callback(hObject, eventdata, handles)
% hObject    handle to typenameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of typenameedit as text
%        str2double(get(hObject,'String')) returns contents of typenameedit as a double

handles.types(handles.typeselected).name = get(hObject,'string');
set(handles.typemenu,'string',{handles.types.name});
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function typenameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typenameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in typeslist.
function typeslist_Callback(hObject, eventdata, handles)
% hObject    handle to typeslist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns typeslist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from typeslist
contents = get(hObject,'String');
handles.types(handles.typeselected).fields = contents(get(hObject,'Value'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function typeslist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typeslist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edittypesbutton.
function edittypesbutton_Callback(hObject, eventdata, handles)
% hObject    handle to edittypesbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('TO DO\n');
% to do

% --- Executes on button press in addnewtypebutton.
function addnewtypebutton_Callback(hObject, eventdata, handles)
% hObject    handle to addnewtypebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.typeselected = length(handles.types)+1;
names = {handles.types.name};
for i = 1:inf,
  handles.types(handles.typeselected).name = sprintf('Population %d',i);
  if ~ismember(handles.types(handles.typeselected).name,names),
    break;
  end
end
handles.types(handles.typeselected).fields = {'All'};
handles = initializetypes(handles);
guidata(hObject,handles);

% --- Executes on button press in deletetypebutton.
function deletetypebutton_Callback(hObject, eventdata, handles)
% hObject    handle to deletetypebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(handles.types) == 1,
  msgbox('Cannot delete only type');
  return;
end

idx = handles.typeselected;
handles.types(idx) = [];
if idx > length(handles.types),
  handles.typeselected = handles.typeselected - 1;
end
handles = initializetypes(handles);
guidata(hObject,handles);

function behaviornameedit_Callback(hObject, eventdata, handles)
% hObject    handle to behaviornameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of behaviornameedit as text
%        str2double(get(hObject,'String')) returns contents of behaviornameedit as a double
handles.behaviors(handles.behaviorselected).name = get(hObject,'string');
behaviornames = {handles.behaviors.name};
set(handles.behaviormenu,'string',behaviornames);
set(handles.propbehaviormenu,'string',behaviornames);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function behaviornameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to behaviornameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
delete(hObject);

function basename = getbasename(path)

i = find(path == '/',1,'last');
if isempty(i),
  basename = path;
else
  basename = path(i+1:end);
end


% --- Executes on button press in perframebutton.
function perframebutton_Callback(hObject, eventdata, handles)
% hObject    handle to perframebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of perframebutton

v = get(hObject,'value');
if v,
  if strcmpi(props(handles.propselected).during,'allframes'),
    set(handles.propbehaviormenu,'enable','off');
  else
    set(handles.propbehaviormenu,'enable','on');
  end
  set(handles.perframepanel,'visible','on');
  handles.props(handles.propselected).proptype = 'perframe';
end
guidata(hObject,handles);

% --- Executes on button press in durationbutton.
function durationbutton_Callback(hObject, eventdata, handles)
% hObject    handle to durationbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of durationbutton
v = get(hObject,'value');
if v,
  set(handles.propbehaviormenu,'enable','on');
  set(handles.perframepanel,'visible','off');
  handles.props(handles.propselected).proptype = 'duration';
end
guidata(hObject,handles);

% --- Executes on button press in fractimebutton.
function fractimebutton_Callback(hObject, eventdata, handles)
% hObject    handle to fractimebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fractimebutton
v = get(hObject,'value');
if v,
  set(handles.propbehaviormenu,'enable','on');
  handles.props(handles.propselected).proptype = 'fractime';
  set(handles.perframepanel,'visible','off');
end
guidata(hObject,handles);

% --- Executes on button press in freqbutton.
function freqbutton_Callback(hObject, eventdata, handles)
% hObject    handle to freqbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of freqbutton
v = get(hObject,'value');
if v,
  set(handles.propbehaviormenu,'enable','on');
  handles.props(handles.propselected).proptype = 'frequency';
  set(handles.perframepanel,'visible','off');
end
guidata(hObject,handles);

% --- Executes on button press in isbehaviorbutton.
function isbehaviorbutton_Callback(hObject, eventdata, handles)
% hObject    handle to isbehaviorbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isbehaviorbutton
v = get(hObject,'value');
if v,
  set(handles.propbehaviormenu,'enable','on');
  set(handles.averagingmenu,'enable','on');
  handles.props(handles.propselected).during = 'during';
end
guidata(hObject,handles);


% --- Executes on button press in isnotbehaviorbutton.
function isnotbehaviorbutton_Callback(hObject, eventdata, handles)
% hObject    handle to isnotbehaviorbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isnotbehaviorbutton
v = get(hObject,'value');
if v,
  set(handles.propbehaviormenu,'enable','on');
  set(handles.averagingmenu,'enable','on');
  handles.props(handles.propselected).during = 'invert';
end
guidata(hObject,handles);

% --- Executes on button press in allframesbutton.
function allframesbutton_Callback(hObject, eventdata, handles)
% hObject    handle to allframesbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allframesbutton
v = get(hObject,'value');
if v,
  if strcmpi(handles.props(handles.propselected).proptype,'perframe'),
    set(handles.propbehaviormenu,'enable','off');
  else
    set(handles.propbehaviormenu,'enable','on');
  end
  set(handles.averagingmenu,'enable','off');
  handles.props(handles.propselected).during = 'allframes';
end
guidata(hObject,handles);


% --- Executes on button press in chooseindepbox.
function chooseindepbox_Callback(hObject, eventdata, handles)
% hObject    handle to chooseindepbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chooseindepbox
handles.autochoose.isindep = get(hObject,'value');
guidata(hObject,handles);


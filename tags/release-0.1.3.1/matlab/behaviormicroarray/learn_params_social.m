% learn_params_social
setuppath;

%% set all defaults

px2mm = 4;
labelmatname = '';
labelmatpath = '';
paramsmatpath = '';
minseqlengthorder = 5;

%% load settings

pathtolearnparams = which('learn_params');
savedsettingsfile = strrep(pathtolearnparams,'learn_params.m','.learnparamsrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% get the labeled data

labelmode = questdlg('Label new data or load labeled data?','Data source',...
  'Label','Load','Cancel','Label');

if strcmpi(labelmode,'cancel'),
  return;
end

if strcmpi(labelmode,'load'),
  fprintf('Choose a mat file with the labeled data\n');
  labelmatname = [labelmatpath,labelmatname];
  [labelmatname,labelmatpath,filterindex] = uigetfile('.mat','Choose labeled data file',labelmatname);
  if isnumeric(labelmatname) && labelmatname == 0,
    return;
  end
  load(labelmatname);
else  
  label_data_social;
end

%% crop out labeled data

clear labeledtrx;
nmovies = length(movienames);
starts = {};
ends = {};
load(labelmatname);
for i = 1:nmovies,
  load(matnames{i});
  % assume this has already been done
  %trx = process_data(trx,matnames{i},movienames{i});
  %trx = process_data_crabwalks(trx);
  nflies = length(trx);
  for j = 1:length(fliestolabel{i}),
    fly = fliestolabel{i}(j);
    if length(labeledbehavior{i}) < fly || ~isfield(labeledbehavior{i}(fly),'starts'),
      fprintf('No labels for fly %d, movie %d, skipping\n',j,i);
      continue;
    end
    t0 = t0tolabel{i}(j) + trx(fly).firstframe - 1;
    t1 = t1tolabel{i}(j) + trx(fly).firstframe - 1;
    [labeledbehavior{i}(fly).starts,order] =  sort(labeledbehavior{i}(fly).starts);
    labeledbehavior{i}(fly).ends =  labeledbehavior{i}(fly).ends(order);
    labeledbehavior{i}(fly).otherfly = labeledbehavior{i}(fly).otherfly(order);
    if ~isempty(labeledbehavior{i}(fly).starts) && labeledbehavior{i}(fly).starts(1) <= t0,
      fprintf('The first frame for fly %d, movie %d is labeled.\n',j,i);
      fprintf('Assuming that the sequence has been artificially cropped\n');
      t0 = labeledbehavior{i}(fly).ends(1)+1;
      fprintf('Shortening to frames %d to %d\n',t0,t1);
      labeledbehavior{i}(fly).starts(1) = [];
      labeledbehavior{i}(fly).ends(1) = [];
      labeledbehavior{i}(fly).otherfly(1) = [];
    end
    if ~isempty(labeledbehavior{i}(fly).ends) && labeledbehavior{i}(fly).ends(end) >= t1,
      fprintf('The last frame for fly %d, movie %d is labeled.\n',j,i);
      fprintf('Assuming that the sequence has been artificially cropped\n');
      t1 = labeledbehavior{i}(fly).starts(end)-1;
      fprintf('Shortening to frames %d to %d\n',t0,t1);      
      labeledbehavior{i}(fly).starts(end) = [];
      labeledbehavior{i}(fly).ends(end) = [];
      labeledbehavior{i}(fly).otherfly(end) = [];
    end
    fprintf('Processing data for movie %d, fly %d\n',i,j);
    trk = process_data_social(trx,fly,t0,t1);
    
    % split into one track per pair of flies
    newtrx = struct([]);
    fns = fieldnames(trk);
    for ii = 1:length(fns),
      fn = fns{ii};
      kk = 1;
      for jj = 1:nflies,
        if jj == fly, continue; end
        t0curr = max(t0,trx(jj).firstframe);
        t1curr = min(t1,trx(jj).endframe);
        if t0curr >= t1curr,
          continue;
        end
        if size(trk.(fn),1) == nflies && size(trk.(fn),2) >= trk.nframes - 10,
          newtrx(kk).(fn) = trk.(fn)(jj,:);
        else
          newtrx(kk).(fn) = trk.(fn);
        end
        kk = kk + 1;
      end
    end
    if kk == 1, continue; end
    for jj = 1:length(newtrx),
      t0curr = max(t0,trx(jj).firstframe);
      t1curr = min(t1,trx(jj).endframe);
      newtrx(jj) = GetPartOfTrack(newtrx(jj),t0curr,t1curr);
    end
    if ~exist('labeledtrx','var'),
      labeledtrx = newtrx;
    else
      fieldsremove = setdiff(fieldnames(labeledtrx),fieldnames(newtrx));
      labeledtrx = rmfield(labeledtrx,fieldsremove);
      fieldsremove = setdiff(fieldnames(newtrx),fieldnames(labeledtrx));
      newtrx = rmfield(newtrx,fieldsremove);
      newtrx = orderfields(newtrx,labeledtrx);
      labeledtrx = [labeledtrx,newtrx];
    end
    
    % create labels for each other fly
    newstarts = cell(1,nflies-1);
    newends = cell(1,nflies-1);
    kk = 1;
    for jj = 1:nflies,
      if jj == fly,
        continue;
      end
      idx = labeledbehavior{i}(fly).otherfly == jj;
      newstarts{kk} = labeledbehavior{i}(fly).starts(idx);
      newends{kk} = labeledbehavior{i}(fly).ends(idx);
      kk = kk + 1;
    end
    starts = [starts,newstarts];
    ends = [ends,newends];
  end
end

fprintf('Done processing data\n');

% put in format used by systematic_learn_params
datalearn = labeledtrx;
labels = struct('starts',starts,'ends',ends,'notes',cell(size(starts)));
for pair = 1:length(datalearn),
  t0 = datalearn(pair).firstframe;
  datalearn(pair).firstframe = 1;
  datalearn(pair).endframe = datalearn(pair).nframes;
  datalearn(pair).f2i = @(f) f;
  labels(pair).starts = labels(pair).starts - t0 + 1;
  labels(pair).ends = labels(pair).ends - t0 + 1;
end

%% choose file to save params to

if ~exist('ds','var')
  ds = datestr(now,30);
end
paramsmatname = sprintf('learnedparams_%s.mat',ds);
paramsmatname = [labelmatpath,paramsmatname];
[paramsmatname,paramsmatpath] = uiputfile('.mat','Choose file to save parameters to',paramsmatname);

%% choose parameters

try
  load(paramsmatname,'params');
catch
end
fns = fieldnames(labeledtrx);
if exist('params','var'),
  
  % convert radians to degrees
  for i = 2:2:length(params.options),
    for j = 1:2:length(params.options{i}),
      if isfield(trx(1).units,params.options{i}{j}), 
        n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).num));
        if n > 1,
          params.options{i}{j+1} = params.options{i}{j+1}*(180/pi)^n;
        end
        n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).den));
        if n > 1,
          params.options{i}{j+1} = params.options{i}{j+1}/(180/pi)^n;
        end
      end
    end
  end
  params = chooseproperties(fns,params);
else
  params = chooseproperties(fns);
end

% convert degrees to radians
for i = 2:2:length(params.options),
  for j = 1:2:length(params.options{i}),
    if isfield(trx(1).units,params.options{i}{j}), 
      n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).num));
      if n > 1,
        params.options{i}{j+1} = params.options{i}{j+1}*(180/pi)^n;
      end
      n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).den));
      if n > 1,
        params.options{i}{j+1} = params.options{i}{j+1}/(180/pi)^n;
      end
    end
  end
end

save(paramsmatname,'labelmatname','params');

%% are there any parameters that are free?

% remove the forced parameters
tmpparams = params;
if ~isempty(params.options),
  paramnames = params.options(1:2:end);
  paramvalues = params.options(2:2:end);
  for j = 1:length(paramnames),
    if ~isempty(strfind(paramnames{j},'force')),
      tmp = strrep(paramnames{j},'force','');
      tmp = [tmp,'fns'];
      tmp3 = tmpparams.(tmp);
      tmp2 = paramvalues{j}(1:2:end);
      tmpparams.(tmp) = setdiff(tmp3,tmp2);
    end
  end
end

% check if there is anything left
isunknown = false;
fns = fieldnames(tmpparams);
for i = 1:length(fns),
  if ismember(fns{i},{'options','minr','maxr'}),
    continue;
  end
  if ~isempty(tmpparams.(fns{i})),
    isunknown = true;
    break;
  end
end
isunknown = isunknown || (params.minr < params.maxr);

if ~isunknown,
  behaviorparams = struct('minx',struct,'maxx',struct,'minxclose',struct,...
    'maxxclose',struct,'minsumx',struct,'maxsumx',struct,'minmeanx',struct,...
    'maxmeanx',struct','r',params.minr,'minseqlength',nan,'maxseqlength',inf);
  for i = 1:2:length(params.options),
    n = params.options{i};
    if isempty(strfind(n,'force')),
      continue;
    end
    n = n(6:end);
    v = params.options{i+1};
    for j = 1:2:length(v),
      behaviorparams.(n).(v{i}) = v{i+1};
    end
  end
  behaviorparams.r = params.minr;
  lengths = getstructarrayfield(labels,'ends')-getstructarrayfield(labels,'starts');
  behaviorparams.minseqlength = prctile(lengths,minseqlengthorder);
  return;
end

%% set up to learn optimal parameters

% some parameters that are not yet changeable by the user
params.options{end+1} = 'maxr';
params.options{end+1} = 10;

costparams.prepend = [0,.2];
costparams.append = [0,.2];
costparams.connect = [1,.2];
costparams.spurious = [1,.2];
costparams.preshorten = [0,.2];
costparams.postshorten = [0,.2];
costparams.lost = [1,.2];
costparams.split = [1,.2];
costparams.maxextend = 5;

params.options{end+1} = 'costparams';
params.options{end+1} = costparams;

params.options{end+1} = 'nsamples';
params.options{end+1} = 100;
params.options{end+1} = 'minseqlengthorder';
params.options{end+1} = minseqlengthorder;

%% learn parameters

behaviorparams = systematic_learn_params2(datalearn,labels,params.minxfns,...
  params.maxxfns,params.minxclosefns,params.maxxclosefns,...
  params.minsumxfns,params.maxsumxfns,params.minmeanxfns,...
  params.maxmeanxfns,params.options{:});

save('-append',paramsmatname,'behaviorparams');
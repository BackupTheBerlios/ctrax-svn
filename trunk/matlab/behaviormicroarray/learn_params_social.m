% learn_params
setuppath;

%% set all defaults

labelmatname = '';
labelmatpath = '';
paramsmatpath = '';
minseqlengthorder = 5;
succeeded = false;

%% load settings

pathtolearnparams = which('learn_params');
savedsettingsfile = strrep(pathtolearnparams,'learn_params.m','.learnparamsrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% get the labeled data

fprintf('Do you want to label new data, or load pre-labeled data from a mat file?\n');
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
  fprintf('Labeled data loaded\n');
else  
  [labeldata_succeeded,labelmatname] = label_data_social_f();
  if ~labeldata_succeeded,
    return;
  end
  load(labelmatname);

end

%% compute the pairwise per-frame parameters for each fly labeled

clear labeledtrx;
nmovies = length(movienames);
starts = {};
ends = {};
for i = 1:nmovies,
  fprintf('Getting tracks for movie %s\n',movienames{i});
  load(matnames{i});
  for j = 1:length(fliestolabel{i}),
    fly = fliestolabel{i}(j);
    if length(labeledbehavior{i}) < fly || ~isfield(labeledbehavior{i}(fly),'starts'),
      fprintf('No labels for fly %d, movie %d, skipping\n',j,i);
      continue;
    end
    t0 = t0tolabel{i}(j) + trx(fly).firstframe - 1;
    t1 = t1tolabel{i}(j) + trx(fly).firstframe - 1;
    labeledbehavior{i}(fly).starts =  sort(labeledbehavior{i}(fly).starts);
    labeledbehavior{i}(fly).ends =  sort(labeledbehavior{i}(fly).ends);
    if ~isempty(labeledbehavior{i}(fly).starts) && labeledbehavior{i}(fly).starts(1) <= t0,
      fprintf('The first frame for fly %d, movie %d is labeled.\n',j,i);
      fprintf('Assuming that the sequence has been artificially cropped\n');
      t0 = labeledbehavior{i}(fly).ends(1)+1;
      fprintf('Shortening to frames %d to %d\n',t0,t1);
      labeledbehavior{i}(fly).starts(1) = [];
      labeledbehavior{i}(fly).ends(1) = [];
    end
    if ~isempty(labeledbehavior{i}(fly).ends) && labeledbehavior{i}(fly).ends(end) >= t1,
      fprintf('The last frame for fly %d, movie %d is labeled.\n',j,i);
      fprintf('Assuming that the sequence has been artificially cropped\n');
      t1 = labeledbehavior{i}(fly).starts(end)-1;
      fprintf('Shortening to frames %d to %d\n',t0,t1);      
      labeledbehavior{i}(fly).starts(end) = [];
      labeledbehavior{i}(fly).ends(end) = [];
    end
    pairwisematnames{i}{j} = strrep(matnames{i},'.mat',...
      sprintf('_pairwiseperframefly%02d%t%05dto%05d.mat',fly,t0,t1));
    fprintf('Computing pairwise per-frame properties for movie %s, fly %d, frames %d to %d\n',...
      matnames{i},fly,t0,t1);
    pairwise_succeeded = compute_pairwise_perframe_stats('matname',matnames{i},...
      'fly',fly,'t0',t0,'t1',t1,'savename',pairwisematnames{i}{j});
    if ~pairwise_succeeded,
      return;
    end

    %trk = process_data_crabwalks(trk);
    trk = GetPartOfTrack(trk,t0,t1);
    if ~exist('labeledtrx','var'),
      labeledtrx = trk;
    else
      fieldsremove = setdiff(fieldnames(labeledtrx),fieldnames(trk));
      labeledtrx = rmfield(labeledtrx,fieldsremove);
      fieldsremove = setdiff(fieldnames(trk),fieldnames(labeledtrx));
      trk = rmfield(trk,fieldsremove);
      labeledtrx(end+1) = labeledtrx(1);
      fns = fieldnames(labeledtrx);
      for tmp = 1:length(fns),
        labeledtrx(end).(fns{tmp}) = trk.(fns{tmp});
      end
    end
    starts{end+1} = labeledbehavior{i}(fly).starts;
    ends{end+1} = labeledbehavior{i}(fly).ends;
  end
end

% put in format used by systematic_learn_params
datalearn = labeledtrx;
labels = struct('starts',{},'ends',{},'notes',{});
for i = 1:nmovies,
  labels = [labels,labeledbehavior{i}(fliestolabel{i})];
end
for fly = 1:length(datalearn),
  t0 = datalearn(fly).firstframe;
  datalearn(fly).firstframe = 1;
  datalearn(fly).endframe = datalearn(fly).nframes;
  datalearn(fly).f2i = @(f) f;
  labels(fly).starts = labels(fly).starts - t0 + 1;
  labels(fly).ends = labels(fly).ends - t0 + 1;
end

%% choose file to save params to

fprintf('Choose a file to save the learned behavior parameters to.\n');
if ~exist('ds','var')
  ds = datestr(now,30);
end
paramsmatname = sprintf('learnedparams_%s.mat',ds);
paramsmatname = [labelmatpath,paramsmatname];
[paramsmatname,paramsmatpath] = uiputfile('.mat','Choose file to save parameters to',paramsmatname);
if isnumeric(paramsmatname) && paramsmatname == 0,
  return;
end

%% choose parameters

fprintf('Define the type of classifier: choose per-frame properties to bound, set lower and upper limits on the bounds\n');
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
  fprintf('All parameters set exactly, no learning necessary\n');
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
  
  save('-append',paramsmatname,'behaviorparams');

  succeeded = true;
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

fprintf('Learning parameters ... This will take a while. Hit the quit button to stop early.\n');

behaviorparams = systematic_learn_params2(datalearn,labels,params.minxfns,...
  params.maxxfns,params.minxclosefns,params.maxxclosefns,...
  params.minsumxfns,params.maxsumxfns,params.minmeanxfns,...
  params.maxmeanxfns,params.options{:});

save('-append',paramsmatname,'behaviorparams');

succeeded = true;
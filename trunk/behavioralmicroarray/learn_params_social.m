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
  helpmsg = 'Choose a mat file containing the pre-labeled data to load.';
  labelmatname = [labelmatpath,labelmatname];
  [labelmatname,labelmatpath,filterindex] = uigetfilehelp('.mat','Choose labeled data file',labelmatname,'helpmsg',helpmsg);
  if isnumeric(labelmatname) && labelmatname == 0,
    return;
  end
  labelmatname = [labelmatpath,labelmatname];
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

nmovies = length(movienames);
starts = {};
ends = {};
savenames = cell(1,nmovies);
for i = 1:nmovies,
  fprintf('Loading tracks from %s\n',matnames{i});
  [trx,matnames{i},loadsucceeded] = load_tracks(matnames{i},movienames{i});
  if ~loadsucceeded,
    msgbox('Could not load trx from file %s\n',matnames{i});
    return;
  end
  savenames{i} = cell(1,length(fliestolabel{i}));
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
    fprintf('Computing pairwise per-frame properties for movie %s, fly %d, frames %d to %d\n',...
      matnames{i},fly,t0,t1);
    [pairwise_succeeded,savenamecurr] = compute_perframe_stats_social_f('matname',matnames{i},...
      'fly1',fly,'t0',t0,'t1',t1);
    if ~pairwise_succeeded,
      return;
    end
    savenames{i}{j} = savenamecurr{1};

  end
end

%% put in a format that works with systematic_learn_params
for i = 1:nmovies,
  for j = 1:length(fliestolabel{i}),
    fly = fliestolabel{i}(j);
    load(savenames{i}{j});
    for fly2 = 1:length(pairtrx),
      pairtrx(fly2).f2i = @(f) f - pairtrx(fly2).firstframe + 1;
    end
    [labelscurr] = pairwise2unarylabels(labeledbehavior{i}(fly),pairtrx);
    % we want firstframe to be 1 for use in learn_params
    for fly2 = 1:length(pairtrx),
      if pairtrx(fly2).firstframe > pairtrx(fly2).endframe,
        pairtrx(fly2).endframe = 0;
        pairtrx(fly2).nframes = 0;
      else
        pairtrx(fly2).endframe = pairtrx(fly2).endframe - pairtrx(fly2).firstframe + 1;
      end
      pairtrx(fly2).firstframe = 1;
      pairtrx(fly2).f2i = @(f) f;
      
      % seems to be unhappy with non-speed parameters
      fns = fieldnames(pairtrx(fly2));
      for k = 1:length(fns),
        fn = fns{k};
        if length(pairtrx(fly2).(fn)) == pairtrx(fly2).nframes,
          pairtrx(fly2).(fn) = pairtrx(fly2).(fn)(1:end-1);
        end
      end
    end
    if i == 1 && j == 1,
      datalearn = pairtrx(:)';
      labels = labelscurr;
    else
      datalearn = [datalearn,pairtrx(:)'];
      labels = [labels,labelscurr];
    end
  end
end

%% choose file to save params to

fprintf('Choose a file to save the learned behavior parameters to.\n');
if ~exist('ds','var')
  ds = datestr(now,30);
end
paramsmatname = sprintf('learnedsocialparams_%s.mat',ds);
paramsmatname = [labelmatpath,paramsmatname];
helpmsg = 'Choose mat file to save the learned behavior classifier parameters to';
[paramsmatname,paramsmatpath] = uiputfilehelp('.mat','Choose file to save parameters to',paramsmatname,'helpmsg',helpmsg);
if isnumeric(paramsmatname) && paramsmatname == 0,
  return;
end
paramsmatname = [paramsmatpath,paramsmatname];

%% choose parameters

fprintf('Define the type of classifier: choose per-frame properties to bound, set lower and upper limits on the bounds\n');
try
  load(paramsmatname,'params');
catch
end
fns = fieldnames(datalearn);
if exist('params','var'),
  
  % convert radians to degrees
  for i = 2:2:length(params.options),
    for j = 1:2:length(params.options{i}),
      if isfield(trx(1).units,params.options{i}{j}), 
 
        
        n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).num));
        if n >= 1,
          params.options{i}{j+1} = params.options{i}{j+1}*(180/pi)^n;
        end
        n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).den));
        if n >= 1,
          params.options{i}{j+1} = params.options{i}{j+1}/(180/pi)^n;
        end
        
      end
      
    end
  end
  params = chooseproperties(fns,datalearn(1).units,params);
else
  params = chooseproperties(fns,datalearn(1).units);
end

% convert degrees to radians
for i = 2:2:length(params.options),
  for j = 1:2:length(params.options{i}),
    if isfield(trx(1).units,params.options{i}{j}), 

      n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).num));
      if n >= 1,
        params.options{i}{j+1} = params.options{i}{j+1}/(180/pi)^n;
      end
      n = nnz(strcmpi('rad',trx(1).units.(params.options{i}{j}).den));
      if n >= 1,
        params.options{i}{j+1} = params.options{i}{j+1}*(180/pi)^n;
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
params.fps = datalearn(1).fps;

%% learn parameters

fprintf('Learning parameters ... This will take a while.\n');
fprintf('Hit the quit button to stop early.\n');
fprintf('Typing Control-C will result in all computations being lost.\n');

behaviorparams = systematic_learn_params2(datalearn,labels,params.minxfns,...
  params.maxxfns,params.minxclosefns,params.maxxclosefns,...
  params.minsumxfns,params.maxsumxfns,params.minmeanxfns,...
  params.maxmeanxfns,params.fps,params.options{:});
behaviorparams.issocial = true;

save('-append',paramsmatname,'behaviorparams');

succeeded = true;

%% show the resulting classifications

% show at most maxn1 x maxn2 plots per figure
npairs = length(datalearn);
maxn1 = 3;
maxn2 = 4;
nperpage = maxn1*maxn2;
npages = ceil(npairs/nperpage);
nlastpage = mod(npairs,nperpage);
if nlastpage == 0, nlastpage = nperpage; end

page = 1;
plotcurr = 1;
if npages == 1,
  nperpagecurr = nlastpage;
else
  nperpagecurr = nperpage;
end

for pair = 1:npairs,

  % go to the next page
  if plotcurr == 1,
    figure('position',[360   402   886   519],'units','pixels');;
    if page == npages,
      n2 = ceil(sqrt(nlastpage));
      n1 = ceil(nlastpage / n2);
    else
      n1 = maxn1;
      n2 = maxn2;
    end
    hax = createsubplots(n1,n2,.01);
  end
  
  axes(hax(plotcurr));
  hold on;
    
  % compute and plot the other fly's position
  x = datalearn(pair).x_mm + datalearn(pair).dcenter.*cos(datalearn(pair).theta+datalearn(pair).anglefrom1to2);
  y = datalearn(pair).y_mm + datalearn(pair).dcenter.*sin(datalearn(pair).theta+datalearn(pair).anglefrom1to2);

  % get labeled and detected behaviors for this fly
  idxlabel = false(1,datalearn(pair).nframes);
  for i = 1:length(labels(pair).starts),
    idxlabel(labels(pair).starts(i):labels(pair).ends(i)) = true;
  end
  idxdetect = false(1,datalearn(pair).nframes);
  seg = systematic_detect_event4(datalearn(pair),...
    'event',true(1,datalearn(pair).nframes),...
    behaviorparams.minx,behaviorparams.maxx,...
    behaviorparams.minxclose,behaviorparams.maxxclose,...
    behaviorparams.minsumx,behaviorparams.maxsumx,...
    behaviorparams.minmeanx,behaviorparams.maxmeanx,...
    behaviorparams.r,behaviorparams.minseqlength,...
    behaviorparams.maxseqlength);
  for i = 1:length(seg.t1),
    idxdetect(seg.t1(i):seg.t2(i)) = true;
  end
  
  % plot in and around detected and labeled frames
  idxplotother = oned_binary_imdilate(idxdetect | idxlabel,ones(1,5));
  if ~any(idxplotother),
    ax = axis;  
    text(mean(ax(1:2)),mean(ax(3:4)),'No detected or labeled sequences','horizontalalignment','center');
  else
    
    % plot labels
    for i = 1:length(labels(pair).starts),
      hlabel = plot(datalearn(pair).x_mm(labels(pair).starts(i):labels(pair).ends(i)),...
        datalearn(pair).y_mm(labels(pair).starts(i):labels(pair).ends(i)),...
        'm-','linewidth',8);
    end
    
    % plot lines joining
    [idxstarts,idxends] = get_interval_ends(idxplotother);
    idxplotconnector = false(1,datalearn(pair).nframes);
    for i = 1:length(idxstarts),
      idxplotconnector([idxstarts(i):5:idxends(i)-2,idxends(i)-1]) = true;
    end
    hconnect = plot([x(idxplotconnector);datalearn(pair).x_mm(idxplotconnector)],...
      [y(idxplotconnector);datalearn(pair).y_mm(idxplotconnector)],'--','color',[1,1,0]*.7);
    
    % plot other fly's position
    hother = plot(x(idxplotother),y(idxplotother),'g.-');
    
    % plot the first fly's position
    hmain = plot(datalearn(pair).x_mm(idxplotother),datalearn(pair).y_mm(idxplotother),'k.-');
    
    % plot detected positions
    for i = 1:length(seg.t1),
      hdetect = plot(datalearn(pair).x_mm(seg.t1(i):seg.t2(i)),...
        datalearn(pair).y_mm(seg.t1(i):seg.t2(i)),'c.-');
      plot(datalearn(pair).x_mm(seg.t1(i)),...
        datalearn(pair).y_mm(seg.t1(i)),'go','markerfacecolor','g');
      plot(datalearn(pair).x_mm(seg.t2(i)),...
        datalearn(pair).y_mm(seg.t2(i)),'ro','markerfacecolor','r');
    end
    
    axisalmosttight;
  end
  set(hax(plotcurr),'xtick',[],'ytick',[]);
  
  % increment plot
  plotcurr = plotcurr + 1;
  if plotcurr > nperpagecurr && page ~= npages,
    plotcurr = 1;
    page = page + 1;
    if npages == page,
      nperpagecurr = nlastpage;
    else
      nperpagecurr = nperpage;
    end
  end
  
end

if plotcurr > 1,
  delete(hax(plotcurr:end));
end

fprintf('Illustration of detected and labeled behaviors\n');
fprintf('We plot the main and other fly''s position in and around\n');
fprintf('frames where either the behavior is labeled or detected\n');
fprintf('The main fly is plotted in black, the other fly is\n');
fprintf('plotted in green. We connect corresponding frames in yellow.\n');
fprintf('We plot labeled behaviors in magenta and detected behaviors\n');
fprintf('in cyan. At the start and end of each detected behavior, we plot\n');
fprintf('a green and red circle, respectively. There is a subplot for each\n');
fprintf('labeled pair of flies. We plot nothing in cases where no behaviors\n');
fprintf('are labeled or detected.\n');

%% print out the learned parameters

s = printparams(behaviorparams.minxclose,behaviorparams.maxxclose,datalearn(1).units);
if ~isempty(s),
  fprintf('All frame bounds:\n');
  for i = 1:length(s),
    fprintf(s{i});
  end
  fprintf('\n');
end
s = printparams(behaviorparams.minx,behaviorparams.maxx,datalearn(1).units);
if ~isempty(s),
  fprintf('Near-frame bounds:\n');
  for i = 1:length(s),
    fprintf(s{i});
  end
  fprintf('Near frame radius: %d frames\n',behaviorparams.r);
  fprintf('\n');
end
s = printparams(behaviorparams.minsumx,behaviorparams.maxsumx,datalearn(1).units);
if ~isempty(s),
  fprintf('Sequence sum bounds:\n');
  for i = 1:length(s),
    fprintf(s{i});
  end
  fprintf('\n');
end
s = printparams(behaviorparams.minmeanx,behaviorparams.maxmeanx,datalearn(1).units);
if ~isempty(s),
  fprintf('Sequence mean bounds:\n');
  for i = 1:length(s),
    fprintf(s{i});
  end
  fprintf('\n');
end
fprintf('Sequence length: %f <= seqlength <= %f [frames]\n',behaviorparams.minseqlength,behaviorparams.maxseqlength);
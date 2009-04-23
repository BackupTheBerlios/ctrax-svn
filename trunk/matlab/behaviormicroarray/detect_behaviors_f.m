function [succeeded,savenames] = detect_behaviors_f(varargin)

succeeded = false;
savenames = {};
[trxnames,behmatname] = myparse(varargin,'trxnames',{},'behmatname','');
ds = datestr(now,30);

% trajectory files
if isempty(trxnames),
  fprintf('Choose mat files containing trajectories with per-frame properties to segment.\n');
  [trxnames,trxpath] = uigetfile('*.mat','Choose trajectory mat file(s)','','MultiSelect','on');
  if isnumeric(trxnames),
    return;
  end
  if ~iscell(trxnames),
    trxnames = {trxnames};
  end
  for i = 1:length(trxnames),
    trxnames{i} = [trxpath,trxnames{i}];
  end
else
  if ~iscell(trxnames) && ischar(trxnames),
    trxnames = {trxnames};
  end
  [trxpath,tmp] = split_path_and_filename(trxnames{1});
end

% behavior classifier parameters
if isempty(behmatname),
  fprintf('Choose mat file containing parameters of behavior classifier.\n');
  [behmatname,behpath] = uigetfile('*.mat','Choose behavior parameters mat file',trxpath);
  if isnumeric(behmatname),
    return;
  end
  behmatname = [behpath,behmatname];
end

% load the classifier
if ~exist(behmatname,'file'),
  msgbox(sprintf('File %s does not exist',behmatname));
  return;
end
tmp = load(behmatname);
if ~isfield(tmp,'behaviorparams'),
  msgbox(sprintf('File %s does not contain the variable behaviorparams',behmatname));
  return;
end
behaviorparams = tmp.behaviorparams;
if ~isfield(behaviorparams,'issocial'),
  behaviorparams.issocial = false;
end

% get names of per-frame properties the classifier depends on
paramnames = {'minx','maxx','minxclose','maxxclose','minsumx','maxsumx','minmeanx','maxmeanx'};
perframeprops = {};
for i = 1:length(paramnames),
  if ~isfield(behaviorparams,paramnames{i}),
    msgbox(sprintf('behaviorparams loaded from %s does not have field %s',behmatname,paramnames{i}));
    return;
  end
  perframeprops = union(perframeprops,fieldnames(behaviorparams.(paramnames{i})));
end

% loop through all the movies to perform the classification on
nmovies = length(trxnames);
for i = 1:nmovies,
  
  % load in the trajectories
  [trx,loadname,loadsucceeded] = load_tracks(trxnames{i});
  if ~loadsucceeded
    b = questdlg(sprintf('Could not load trajectories from %s, skip it or abort?',trxnames{i}),'Load error','Skip','Abort','Skip');
    if strcmpi(b,'abort'),
      return;
    end
    continue;
  else
    trxnames{i} = loadname;
  end
  
  % check for the per-frame properties
  missingfns = setdiff(perframeprops,fieldnames(trx));
  if ~isempty(missingfns),
  
    % which process function(s) do we need to call?
    [processdata_fcns,stillmissingfns] = perframeprop2processfcn(missingfns,behaviorparams.issocial);
    
    % properties we don't know how to compute?
    if ~isempty(stillmissingfns),
      s = ['Behavior classifier relies on the following properties; unknown how to compute: ',sprintf('%s ',stillmissingfns{:})];
      b = questdlg(s,'Missing per-frame properties','Skip','Abort','Skip');
      if strcmpi(b,'abort'),
        return;
      end
      continue;
    end
    
    % pairwise properties?
    if ismember('process_data_pairs',processdata_fcns),
      b = questdlg('Social behavior classifier, but pairwise perframe parameters have not yet been computed for file. Skip or abort?',...
        'Missing pairwise per-frame parameters','Skip','Abort','Skip');
      if strcmpi(b,'skip'),
        continue;
      else
        return;
      end
    end
    
    % process
    fprintf('The following perframe properties are missing: ');
    fprintf('%s ',missingfns{:});
    fprintf('\n');
    fprintf('The following functions must be called to compute them: ');
    fprintf('%s ',processdata_fcns{:});
    fprintf('\nCompute now?\n');
    b = questdlg('Missing per-frame properties. Compute now?','Missing per-frame properties','Compute','Skip','Abort','Compute');
    if strcmpi(b,'skip'),
      continue;
    elseif strcmpi(b,'abort'),
      return;
    else
      
      % non-social per-frame parameters
      [trxpathcurr,trxnameonly] = split_path_and_filename(trxnames{i});
      docomputearena = ismember('process_data_arena',processdata_fcns);
      docomputeclosest = ismember('process_data_closestfly',processdata_fcns);
      docomputemain = ismember('process_data',processdata_fcns);
      docomputecrabwalk = ismember('process_data_crabwalks',processdata_fcns);
      [perframesucceeded,perframename,trx] = ...
        compute_perframe_stats_f('matname',trxnameonly,'matpath',trxpathcurr,...
        'docomputemain',docomputemain,'docomputecrabwalk',docomputecrabwalk,...
        'docomputearena',docomputearena,'docomputeclosest',docomputeclosest);
      if ~perframesucceeded,
        b = questdlg('Computing per-frame properties failed. Skip or abort?','Compute per-frame failed','Skip','Abort','Skip');
        if strcmpi(b,'skip'),
          continue;
        else
          return;
        end
      end
      trxnames{i} = perframename;
      
    end
    
  end
  
  % check to see if fps matches
  fpsmovie = trx(1).fps;
  fpsbeh = behaviorparams.fps;
  if fpsmovie ~= fpsbeh,
    % convert to the new fps
    fprintf('Converting from behaviorparams.fps = %f to trajectory''s fps = %f\n',fpsbeh,fpsmovie);
    behaviorparams.r = behaviorparams.r*fpsmovie/fpsbeh;
    if round(behaviorparams.r) ~= behaviorparams.r,
      fprintf('Need to round near-frame radius from %f to %d\n',behaviorparams.r,round(behaviorparams.r));
      behaviorparams.r = round(behaviorparams.r);
    end
    behaviorparams.minseqlength = behaviorparams.minseqlength*fpsmovie/fpsbeh;
    behaviorparams.maxseqlength = behaviorparams.maxseqlength*fpsmovie/fpsbeh;
    behaviorparams.fps = fpsmovie;
  end
  
  % segment
  fprintf('Segmenting...\n');
  nflies = length(trx);
  for fly = 1:nflies
    
    if behaviorparams.issocial,
      fprintf('Classifying behavior for pair %d, trajectory file %s\n',fly,trxnames{i});
    else
      fprintf('Classifying behavior for fly %d, trajectory file %s\n',fly,trxnames{i});
    end
    segcurr = systematic_detect_event4(trx(fly),...
      'event',true(1,trx(fly).nframes),...
      behaviorparams.minx,behaviorparams.maxx,...
      behaviorparams.minxclose,behaviorparams.maxxclose,...
      behaviorparams.minsumx,behaviorparams.maxsumx,...
      behaviorparams.minmeanx,behaviorparams.maxmeanx,...
      behaviorparams.r,behaviorparams.minseqlength,...
      behaviorparams.maxseqlength,behaviorparams.fps);
    if fly == 1,
      seg = segcurr;
    else
      seg(fly) = segcurr;
    end
  end
  
  % save results to file
  savenamecurr = strrep(trxnames{i},'.mat',sprintf('_seg_%s.mat',ds));
  fprintf('Choose mat file to export segmentation results for trx file %s\n',trxnames{i});
  while true,
    [savename,savepath] = uiputfile('*.mat',sprintf('Save segmentations for %s',trxnames{i}),savenamecurr);
    if ~ischar(savename),
      b = questdlg('Are you sure you don''t want to save the results? Computation will be lost.','Really cancel?','Yes','No','No');
      if strcmpi(b,'yes'),
        return;
      end
    else
      break;
    end
  end
  savename = [savepath,savename];
  trxname = trxnames{i};
  save(savename,'seg','behmatname','trxname');
  savenames{end+1} = savename;

  [plotsucceeded,plot_figs] = plot_detectbehaviors(trx,seg,behaviorparams.issocial);
  fprintf('Figures for segmentations for file %s: ',trxnames{i});
  fprintf('%d ',plot_figs);
  fprintf('\n');
  
end

succeeded = true;


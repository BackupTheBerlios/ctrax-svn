%% label_data_social

%% user preferences

lastmoviename = '';
lastmoviepath = '';
lastmatname = '';
lastmatpath = '';
labelmatpath = '';
movieexts = {'*.fmf','*.sbfmf','*.avi'}';

pathtolearnparams = which('learn_params');
savedsettingsfile = strrep(pathtolearnparams,'learn_params.m','.learnparamsrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% choose some movies to label

movienames = {};
moviepaths = {};
matnames = {};
matpaths = {};
while true,
  fprintf('Choose a movie to label\n');
  lastmoviename = [lastmoviepath,lastmoviename];
  [lastmoviename,lastmoviepath,filterindex] = uigetfile(movieexts,'Choose movie file',lastmoviename);
  if isnumeric(lastmoviename) && lastmoviename == 0,
    break;
  end
  movienames{end+1} = lastmoviename;
  moviepaths{end+1} = lastmoviepath;
  lastmatname = [lastmoviepath,strrep(lastmoviename,movieexts{filterindex}(2:end),'.mat')];
  [lastmatname,lastmatpath] = uigetfile('*.mat',sprintf('Choose per-frame stats mat file for movie %s',lastmoviename),lastmatname);
  if isnumeric(lastmatname) && lastmatname == 0,
    break;
  end
  if any(strcmp(lastmatname,matnames)&strcmp(lastmatpath,matpaths)),
    fprintf('Duplicate movie selected, ignoring.\n');
    continue;
  end
  matnames{end+1} = lastmatname;
  matpaths{end+1} = lastmatpath;
  fprintf('Movie:   %s%s\nMatfile: %s%s\n\n',moviepaths{end},movienames{end},matpaths{end},matnames{end});

  b = questdlg('Do you want to add more movies?','More?','Yes');
  if ~strcmpi('yes',b),
    break;
  end
end

nmovies = length(matnames);
movienames = movienames(1:nmovies);
moviepaths = moviepaths(1:nmovies);

lastmoviename = movienames{end};
lastmatname = matnames{end};
lastmoviepath = moviepaths{end};
lastmatpath = matpaths{end};
if exist(savedsettingsfile,'file'),
  save('-append',savedsettingsfile,'lastmoviename','lastmatname','lastmoviepath','lastmatpath');
else
  save(savedsettingsfile,'lastmoviename','lastmatname','lastmoviepath','lastmatpath');
end

for i = 1:nmovies,
  movienames{i} = [moviepaths{i},movienames{i}];
  matnames{i} = [matpaths{i},matnames{i}];
end

%% how many frames, flies / movie

% get stats for each movie
nfliespermovie = zeros(1,nmovies);
nmalepermovie = zeros(1,nmovies);
nfemalepermovie = zeros(1,nmovies);
nunknownpermovie = zeros(1,nmovies);
nframesperfly = cell(1,nmovies);
sexperfly = cell(1,nmovies);
ds = datestr(now,30);
skipped = [];
for i = 1:nmovies,
  datacurr = load(matnames{i});
  fprintf('Checking file %s\n',matnames{i});
  if ~isfield(datacurr,'trx'),
    if isfield(datacurr,'ntargets'),
      [datacurr.trx,matnames{i}] = cleanup_ctrax_data(matnames{i},movienames{i},datacurr,ds);
    else
      b = questdlg(sprintf('Could not load data from %s.',matnames{i}),'Data error','Skip','Abort','Skip');
      if strcmpi(b,'abort'),
        return;
      else
        skipped(end+1) = i;
        continue;
      end
    end
  end
  
  % convert units if necessary
  if ~isfield(datacurr.trx,'fps') || ~isfield(datacurr.trx,'pxpermm'),
    savedsettingsfile0 = savedsettingsfile;
    ISAUTOMATIC = true;
    ISMATNAME = true;
    matname = matnames{i};
    ISMOVIENAME = true;
    moviename = movienames{i};
    convert_units;
    matnames{i} = savename;
    savedsettingsfile = savedsettingsfile0;
    datacurr = load(matnames{i});
  end
  
  nflies1 = length(datacurr.trx);
  nframesperfly1 = getstructarrayfield(datacurr.trx,'nframes');
  if isfield(datacurr.trx,'sex'),
    sex1 = getstructarrayfield(datacurr.trx,'sex');
  else
    sex1 = repmat('?',[1,nflies1]);
  end
  nfliespermovie(i) = nflies1;
  nframesperfly{i} = nframesperfly1;
  nfemalepermovie(i) = nnz(sex1=='F');
  nmalepermovie(i) = nnz(sex1=='M');
  nunknownpermovie(i) = nnz(sex1=='?');
  sexperfly{i} = sex1;
  fprintf('%d flies:\n',nflies1);
  for j = 1:nflies1,
    fprintf('  Fly %d: %d frames, sex = %s\n',j,nframesperfly{i}(j),sexperfly{i}(j));
  end
end

% remove bad movies
if ~isempty(skipped),
  nfliespermovie(skipped) = [];
  nmalepermovie(skipped) = [];
  nfemalepermovie(skipped) = [];
  nunknownpermovie(skipped) = [];
  nframesperfly(skipped) = [];
  sexperfly(skipped) = [];
  movienames(skipped) = [];
  matnames(skipped) = [];
end

%% choose flies to label, number of frames to show

fprintf('Total of %d movies.\n',nmovies);
fprintf('Average flies per movie: %.1f.\n',mean(nfliespermovie));
fprintf('Average female flies per movie: %.1f.\n',mean(nfemalepermovie));
fprintf('Average male flies per movie: %.1f.\n',mean(nmalepermovie));
fprintf('Average unknown sex flies per movie: %.1f\n',mean(nunknownpermovie));
fprintf('Average frames per fly: %.1f\n',mean(cell2mat(nframesperfly)));
if ~exist('nfliestolabelpermovie','var')
  nfliestolabelpermovie = min(nfliespermovie);
end
if ~exist('nframestolabelperfly','var'),
  nframestolabelperfly = 1000;
end
while true,
  b = inputdlg({'Number of flies to label per movie:','Number of frames to label per fly',...
    'Sex of flies to label (choose F=female, M=male, ?=any)'},...
    'Frames to Label',1,{num2str(nfliestolabelpermovie),num2str(nframestolabelperfly),'?'});
  if isempty(b),
    return;
  end
  nfliestolabelpermovie = str2double(b{1});
  if isnan(nfliestolabelpermovie) || nfliestolabelpermovie <= 0 || round(nfliestolabelpermovie) ~= nfliestolabelpermovie,
    fprintf('Number of flies must be a positive whole number\n');
    continue;
  end
  nframestolabelperfly = str2double(b{2});
  if isnan(nframestolabelperfly) || nframestolabelperfly <= 0 || round(nframestolabelperfly) ~= nframestolabelperfly,
    fprintf('Number of frames to label per fly must be a positive whole number\n');
    continue;
  end
  if ~ismember(lower(b{3}),{'m','f','?'}),
    fprintf('Sex must be M, F, or ?\n');
    continue;
  end
  sextolabel = upper(b{3});
  break;
end

save('-append',savedsettingsfile,'nframestolabelperfly','nfliestolabelpermovie');

% choose flies to label at equal intervals
fliestolabel = cell(1,nmovies);
for i = 1:nmovies,
  if nfliestolabelpermovie > nfliespermovie(i),
    fliestolabel{i} = 1:nfliespermovie(i);
  else
    
    % how good is each fly to label:
    % we first perfer flies of the right sex + with right number of frames, 
    % then flies of the right sex + not enough frames,
    % then enough frames, not the right sex,
    % then not enough frames, not the right sex
    score = min(nframestolabelperfly,nframesperfly{i});
    if sextolabel ~= '?',
      score = score + nframestolabelperfly*9*double(sexperfly{i} == sextolabel);
    end

    % find the nfliestolabelpermovie^th highest score
    [tmp,order] = sort(-score);
    minscore = score(order(nfliestolabelpermovie));
    
    % choose flies with at least that score
    okflies = score >= minscore;
    tmp = find(okflies);
    fliestolabel{i} = tmp(round(linspace(1,length(tmp),nfliestolabelpermovie)));
    
  end
  
  fprintf('Stats of flies chosen to label for movie %s:\n',movienames{i});
  for j = 1:length(fliestolabel{i}),
    k = fliestolabel{i}(j);
    fprintf('Fly %d: sex=%s, nframes=%d\n',k,sexperfly{i}(k),nframesperfly{i}(k));
  end
  
end

% choose frames to label
t0tolabel = cell(1,nmovies);
t1tolabel = cell(1,nmovies);
for i = 1:nmovies,
  t0tolabel{i} = zeros(1,length(fliestolabel{i}));
  t1tolabel{i} = zeros(1,length(fliestolabel{i}));
  for j = 1:length(fliestolabel{i}),
    fly = fliestolabel{i}(j);
    if nframesperfly{i}(fly) <= nframestolabelperfly,
      t0tolabel{i}(j) = 1;
      t1tolabel{i}(j) = nframesperfly{i}(fly);
    else
      t0tolabel{i}(j) = randsample(nframesperfly{i}(fly)-nframestolabelperfly+1,1);
      t1tolabel{i}(j) = t0tolabel{i}(j) + nframestolabelperfly - 1;
    end
  end
end

%% save file name

labelmatname = sprintf('%slabeleddata%s.mat',labelmatpath,ds);
[labelmatname,labelmatpath] = uiputfile('*.mat','Save labels to...',labelmatname);
save('-append',savedsettingsfile,'labelmatname','labelmatpath');

labelmatname = [labelmatpath,labelmatname];
save(labelmatname,'movienames','matnames','fliestolabel','t0tolabel','t1tolabel');

%% label

fprintf('Labeling begun. To restart labeling later, run the script "restart_label_data_social"\n');
fprintf('and select as input %s.\n\n',labelmatname);

clear labeledbehavior;
for moviei = 1:nmovies,
  load(matnames{moviei});
  % assuming already processed now
  %trx = process_data(trx,matnames{moviei},movienames{moviei});
  trx = mindist2fly(trx);
  for flyi = 1:length(fliestolabel{moviei}),
    fly = fliestolabel{moviei}(flyi);
    t0 = t0tolabel{moviei}(flyi);
    t1 = t1tolabel{moviei}(flyi);
    t0 = t0 + trx(fly).firstframe - 1;
    t1 = t1 + trx(fly).firstframe - 1;
    for i = 1:length(trx),
      trx(i) = GetPartOfTrack(trx(i),t0,t1);
    end
    [labeledbehavior{moviei}(fly).starts,labeledbehavior{moviei}(fly).ends,...
      labeledbehavior{moviei}(fly).otherfly,labeledbehavior{moviei}(fly).notes] = ...
      labelbehaviorssocial(trx,fly,movienames{moviei},...
      [],[],[],[],'max(0,50-trk.mindist2fly)');
    save('-append',labelmatname,'labeledbehavior');
    fprintf('Done labeling fly %d of movie %s,\n(fly %d / %d for movie %d / %d)\n',fly,movienames{moviei},...
      flyi,length(fliestolabel{moviei}),moviei,nmovies);
    if ~(moviei == nmovies && flyi == length(fliestolabel{moviei})),
      b = questdlg('Keep labeling?','Keep labeling?','Next fly','Quit','Next fly');
      if strcmpi(b,'quit'),
        return;
      end
    end
  end
end

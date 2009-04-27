function [succeeded,savenames] = compute_perframe_stats_social_f(varargin)

% parse inputs
setuppath;

[matname,flies1,t0s,t1s,fov] = ...
  myparse(varargin,'matname',nan,'fly1',[],'t0',[],'t1',[],'fov',pi);

succeeded = false;
savename = '';

ISMATNAME = ischar(matname);

if ~ISMATNAME,
  matname = '';
  matpath = '';
end

%% load settings

pathtocomputeperframestats = which('compute_perframe_stats_social_f');
savedsettingsfile = strrep(pathtocomputeperframestats,'compute_perframe_stats_social_f.m','.computeperframestatssocialrc.mat');
if exist(savedsettingsfile,'file')
  defaultparams = load(savedsettingsfile);
  if ~ISMATNAME,
    if isfield(defaultparams,'matname') && isfield(defaultparams,'matpath'),
      matname = defaultparams.matname;
      matpath = defaultparams.matpath;
    end
  end
end

%% choose a mat file to analyze
if ~ISMATNAME,
  fprintf('Choose mat file to analyze.\n');
  matname = [matpath,matname];
  [matname,matpath] = uigetfile('*.mat','Choose mat file to analyze',matname);
  if isnumeric(matname) && matname == 0,
    return;
  end
else
  [matpath,matname] = split_path_and_filename(matname);
end
fprintf('Matfile: %s%s\n\n',matpath,matname);
  
if exist(savedsettingsfile,'file'),
  save('-append',savedsettingsfile,'matname','matpath');
else
  save(savedsettingsfile,'matname','matpath');
end

matnameonly = matname;
matname = [matpath,matname];

%% get conversion to mm, seconds

[convertunits_succeeded,matname] = ...
  convert_units_f('isautomatic',true,'matname',matnameonly,'matpath',matpath);

if ~convertunits_succeeded,
  return;
end

%% load in the data
fprintf('Checking file %s\n',matname);
[trx,matname,loadsucceeded] = load_tracks(matname);
if ~loadsucceeded,
  msgbox('Could not load data from %s, exiting',matname);
  return;
end

%% set flies1 if not set yet

nflies = length(trx);
if isempty(flies1),
  flies1 = 1:nflies;
end
nflies1 = length(flies1);
if nflies1 ~= length(t0s),
  if length(t0s) == 1,
    t0s = repmat(t0s,[1,nflies1]);
  else
    t0s = ones(1,nflies1);
  end
end
if nflies1 ~= length(t1s),
  if length(t1s) == 1,
    t1s = repmat(t1s,[1,nflies1]);
  else
    t1s = inf(1,nflies1);
  end
end
for i = 1:nflies1,
  t0s(i) = max(t0s(i),trx(flies1(i)).firstframe);
  t1s(i) = min(t1s(i),trx(flies1(i)).endframe);
end

%% compute per-frame stats if necessary

if ~has_perframe_props(fieldnames(trx)),
  fprintf('Computing non-social per-frame properties\n');
  % check to see that process_data has been called
  [computeperframe_succeeded,newmatname,trx] = ...
    compute_perframe_stats_f('matname',matnameonly,'matpath',matpath,'docomputeclosest',false);
  if ~computeperframe_succeeded,
    return;
  end
end

%% compute pairwise per-frame stats

savenames = cell(1,nflies1);
for i = 1:nflies1,
  fprintf('Computing pairwise per-frame stats for fly1 = %d\n',flies1(i));
  savematname = strrep(matname,'.mat',sprintf('_perframepairs_fly%02d_start%05d_end%05d.mat',flies1(i),t0s(i),t1s(i)));
  if exist(savematname,'file'),
    fprintf('File %s exists, not recreating\n',savematname);
  else
    pairtrx = process_data_pairs(trx,flies1(i),t0s(i),t1s(i));
    for j = 1:length(pairtrx),
      pairtrx(j).matname = savematname;
    end
    fprintf('Saving pair per-frame data for fly %d to %s\n',flies1(i),savematname);
    save(savematname,'pairtrx');
  end
  savenames{i} = savematname;
end

succeeded = true;
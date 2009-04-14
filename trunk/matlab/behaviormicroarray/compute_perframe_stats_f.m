function [succeeded,savename,trx] = compute_perframe_stats_f(varargin)

% parse inputs
setuppath;
trx = [];

[matname,matpath,docomputearena,docomputeclosest,fov] = ...
  myparse(varargin,'matname',nan,'matpath',nan,'docomputearena',nan,...
  'docomputeclosest',nan,'fov',nan);

succeeded = false;
savename = '';

ISMATNAME = ischar(matname);
ISMATPATH = ischar(matpath);
if ISMATNAME && ~ISMATPATH,
  matpath = '';
  ISMATNAME = true;
end
ISDOCOMPUTEARENA = ~isnan(docomputearena);
ISDOCOMPUTECLOSEST = ~isnan(docomputeclosest);
ISFOV = ~isnan(fov);

if ~ISDOCOMPUTEARENA,
  docomputearena = false;
end
if ~ISDOCOMPUTECLOSEST,
  docomputeclosest = false;
end
if ~ISFOV,
  fov = pi;
end

if ~ISMATNAME,
  matname = '';
  matpath = '';
end

%% load settings

pathtocomputeperframestats = which('compute_perframe_stats_f');
savedsettingsfile = strrep(pathtocomputeperframestats,'compute_perframe_stats_f.m','.computeperframestatsrc.mat');
if exist(savedsettingsfile,'file')
  defaultparams = load(savedsettingsfile);
  if ~ISMATNAME,
    matname = defaultparams.matname;
    matpath = defaultparams.matpath;
  end
  %fps = defaultparams.fps;
  if ~ISDOCOMPUTEARENA,
    docomputearena = defaultparams.docomputearena;
  end
  if ~ISDOCOMPUTECLOSEST,
    docomputeclosest = defaultparams.docomputeclosest;
  end
  %ppm = defaultparams.ppm;
end

%% choose a mat file to analyze
if ~ISMATNAME,
  fprintf('Choose mat file to analyze.\n');
  matname = [matpath,matname];
  [matname,matpath] = uigetfile('*.mat','Choose mat file to analyze',matname);
  if isnumeric(matname) && matname == 0,
    return;
  end
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

tmp = load(matname);
moviename = '?';
%if ~exist('ds','var'),
ds = datestr(now,30);
%end
fprintf('Checking file %s\n',matname);
if ~isfield(tmp,'trx'),
  if isfield(tmp,'ntargets'),
    [trx,matname] = cleanup_ctrax_data(matname,moviename,tmp,ds);
  else
    msgbox('Could not load data from %s, exiting',matname);
    return;
  end
else
  trx = tmp.trx;
end

%% compute per-frame stats

fprintf('Computing statistics...\n');
trx = process_data(trx,matname);
trx = process_data_crabwalks(trx);

%% compute arena-based per-frame stats?

if ~ISDOCOMPUTEARENA,

  if docomputearena,
    defaultans = 'Yes';
  else
    defaultans = 'No';
  end
  b = questdlg('Do you want to compute distance, angle to arena wall statistics?','Compute arena wall statistics?',defaultans);
  docomputearena = strcmpi(b,'yes');
end
if docomputearena;
  fprintf('Enter name of annotation file to read arena position from.\n');
  annname = [matpath,strrep(matnameonly,'.mat','.ann')];
  if ~exist(annname,'file'),
    annnamefmf = [matpath,strrep(matnameonly,'.mat','.fmf.ann')];
    if exist(annnamefmf,'file')
      annname = annnamefmf;
    else
      annnamesbfmf = [matpath,strrep(matnameonly,'.mat','.sbfmf.ann')];
      if exist(annnamesbfmf,'file')
        annname = annnamesbfmf;
      else
        annnameavi = [matpath,strrep(matnameonly,'.mat','.avi.ann')];
        if exist(annnameavi,'file')
          annname = annnameavi;
        end
      end
    end
  end
  [annname,annpath] = uigetfile('*.ann',sprintf('Annotation file corresponding to %s',matnameonly),annname);
  annname = [annpath,annname];
  if ~exist(annname,'file'),
    fprintf('Annotation file %s does not exist, not computing dist2wall\n',annname);
    docomputearena = false;
  else
    try
      [arena_center_x,arena_center_y,arena_radius] = ...
        read_ann(annname,'arena_center_x','arena_center_y','arena_radius');
    catch
      fprintf('Could not read annotation file %s, not computing dist2wall\n',annname);
      docomputearena = false;
    end
    for fly = 1:length(trx),
      trx(fly).arena.x = arena_center_x;
      trx(fly).arena.y = arena_center_y;
      trx(fly).arena.r = arena_radius;
    end
  end
  trx = process_data_arena(trx);
  
end

save('-append',savedsettingsfile,'docomputearena');

%% compute closest fly statistics

if ~ISDOCOMPUTECLOSEST,

  if docomputeclosest,
    defaultans = 'Yes';
  else
    defaultans = 'No';
  end
  b = questdlg('Do you want to compute statistics relating to closest fly? This may take some time, space.',...
    'Compute closest fly statistics?',defaultans);
  docomputeclosest = strcmpi(b,'yes');
end
if docomputeclosest,
  % get field of view of fly
  if ~ISFOV,
    b = inputdlg({'Field of view of fly in degrees'},'Field of view',1,{num2str(fov*180/pi)});
    if isempty(b),
      fprintf('Field of view not entered, using default value of %.1f\n',fov*180/pi);
    else
      fov = str2double(b{1})*pi/180;
    end
  end
  trx = process_data_closestfly(trx,fov);
end

save('-append',savedsettingsfile,'docomputeclosest');

%% Get name of file to save to 
% 
fprintf('Enter results filename\n');
savename = [matpath,'perframestats_',matnameonly];
[savename, savepath] = uiputfile('*.mat', 'Save results mat name', savename);

savename = [savepath,savename];

%% save results

fprintf('Saving results to file %s\n',savename);
save(savename,'trx');

succeeded = true;
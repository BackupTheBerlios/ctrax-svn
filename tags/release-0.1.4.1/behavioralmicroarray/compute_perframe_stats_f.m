function [succeeded,savename,trx] = compute_perframe_stats_f(varargin)

% parse inputs
setuppath;

[matname,matpath,docomputearena,docomputeclosest,fov,docomputemain,docomputecrabwalk,trx] = ...
  myparse(varargin,'matname',nan,'matpath',nan,'docomputearena',nan,...
  'docomputeclosest',nan,'fov',nan,'docomputemain',true,'docomputecrabwalk',true,'trx',[]);

succeeded = false;
savename = '';

ISTRX = ~isempty(trx);
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
  matnameonly = '';
else
  matnameonly = matname;
end

%% load settings

pathtocomputeperframestats = which('compute_perframe_stats_f');
savedsettingsfile = strrep(pathtocomputeperframestats,'compute_perframe_stats_f.m','.computeperframestatsrc.mat');
if exist(savedsettingsfile,'file')
  defaultparams = load(savedsettingsfile);
  if ~ISMATNAME,
    if isfield(defaultparams,'matname') && isfield(defaultparams,'matpath'),
      matname = defaultparams.matname;
      matpath = defaultparams.matpath;
    end
  end
  %fps = defaultparams.fps;
  if ~ISDOCOMPUTEARENA,
    if isfield(defaultparams,'docomputearena'),
      docomputearena = defaultparams.docomputearena;
    end
  end
  if ~ISDOCOMPUTECLOSEST,
    if isfield(defaultparams,'docomputeclosest'),
      docomputeclosest = defaultparams.docomputeclosest;
    end
  end
  %ppm = defaultparams.ppm;
end

%% choose a mat file to analyze
if ~ISTRX,

  if ~ISMATNAME,
    helpmsg = 'Choose file containing trajectories for which to compute per-frame statistics';
    matname = [matpath,matname];
    [matname,matpath] = uigetfilehelp('*.mat','Choose mat file to analyze',matname,'helpmsg',helpmsg);
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

end

%% get conversion to mm, seconds

if ISTRX,
  [convertunits_succeeded,matname] = ...
    convert_units_f('isautomatic',true,'trx',trx);
else
  [convertunits_succeeded,matname] = ...
    convert_units_f('isautomatic',true,'matname',matnameonly,'matpath',matpath);
end
if ~convertunits_succeeded,
  return;
end

%% load in the data
if ~ISTRX,
  fprintf('Checking file %s\n',matname);
  [trx,matname,loadsucceeded] = load_tracks(matname);
  if ~loadsucceeded,
    msgbox('Could not load data from %s, exiting',matname);
    return;
  end
end

%% compute per-frame stats

fprintf('Computing statistics...\n');
if docomputemain,
  trx = process_data(trx,matname);
end
if docomputecrabwalk,
  trx = process_data_crabwalks(trx);
end

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
  helpmsg = {'We will read in the location of the arena from the Ctrax annotation file.',...
    sprintf('Choose the annotation file corresponding to trx file %s',matname)};
  [annname,annpath] = uigetfilehelp('*.ann',sprintf('Annotation file corresponding to %s',matnameonly),annname,'helpmsg',helpmsg);
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
    if isempty(arena_center_x) || isempty(arena_center_y) || isempty(arena_radius),
      fprintf('Circular arena is not defined in annotation file %s\n',annname);
      docomputearena = false;
    else
      for fly = 1:length(trx),
        trx(fly).arena.x = arena_center_x;
        trx(fly).arena.y = arena_center_y;
        trx(fly).arena.r = arena_radius;
      end
    end
  end
  
  if docomputearena,
    trx = process_data_arena(trx);
  end
  
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
helpmsg = sprintf('Choose file to save single-fly per-frame statistics for trx loaded from %s',matname);
savename = [matpath,'perframestats_',matnameonly];
[savename, savepath] = uiputfilehelp('*.mat', 'Save results mat name', savename,'helpmsg',helpmsg);

savename = [savepath,savename];

%% save results

fprintf('Saving results to file %s\n',savename);
for i = 1:length(trx),
  trx(i).matname = savename;
end
save(savename,'trx');

succeeded = true;
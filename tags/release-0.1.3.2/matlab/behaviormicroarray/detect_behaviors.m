%% detect_behaviors

setuppath;

paramsmatname = '';
matname = '';
fov = pi;

pathtodetectbehaviors = which('detect_behaviors');
savedsettingsfile = strrep(pathtodetectbehaviors,'detect_behaviors.m','.detectbehaviorsrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% choose mat file

fprintf('Choose mat file to detect behaviors for.\n');
matname = [matpath,matname];
[matname,matpath] = uigetfile('*.mat','Choose mat file to detect behaviors for',matname);
if isnumeric(matname) && matname == 0,
  return;
end
fprintf('Matfile: %s%s\n\n',matpath,matname);

if exist(savedsettingsfile,'file'),
  save('-append',savedsettingsfile,'matname','matpath');
else
  save(savedsettingsfile,'matname','matpath');
end

matname0 = matname;
matname = [matpath,matname];

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

if ~isfield(trx,'fps') || ~isfield(trx,'pxpermm'),
  savedsettingsfile0 = savedsettingsfile;
  ISAUTOMATIC = true;
  ISMATNAME = true;
  convert_units;
  matname = savename;
  savedsettingsfile = savedsettingsfile0;
  load(matname);
end

%% choose the behavior parameters file

fprintf('Choose parameters mat file defining the behavior to detect.\n');
[paramsmatname,paramsmatpath] = uigetfile('.mat','Choose parameters mat file defining behavior',paramsmatname);
if isnumeric(paramsmatname) && paramsmatname == 0,
  return;
end
paramsmatname = [paramsmatpath,paramsmatname];
load(paramsmatname);

if ~exist('behaviorparams','var'),
  msgbox('Invalid parameters file. Aborting.');
  return;
end

if exist(savedsettingsfile,'file')
  save('-append',savedsettingsfile,'paramsmatname');
else
  save(savedsettingsfile,'paramsmatname');
end

%% process data, if necessary

fns = fieldnames(trx);
neededfns = {};
tmp = fieldnames(behaviorparams);
for i = 1:length(tmp),
  fn = tmp{i};
  if isstruct(behaviorparams.(fn)),
    tmp2 = fieldnames(behaviorparams.(fn));
    neededfns = union(neededfns,tmp2);
  end
end
missingfns = setdiff(neededfns,fns);
if ~isempty(missingfns),
  fns_processdata = {'dtheta','du_ctr','dv_ctr','corfrac','corisonfly','du_cor',...
    'dv_cor','velmag_ctr','velmag','accmag','signdtheta','absdv_cor','flipdv_cor',...
    'absdtheta','d2theta','absd2theta','smooththeta','smoothdtheta','abssmoothdtheta',...
    'smoothd2theta','abssmoothd2theta','phi','yaw','absyaw'};
  fns_processdatacrabwalks = {'absdtheta_tail','absdu_tail','absdv_tail','absphisideways',...
    'dtheta_tail','du_tail','dv_tail','phisideways'};
  fns_processdataarena = {'absdangle2wall','ddist2wall','dist2wall','theta2wall'};
  fns_processdataclosestfly = {'absanglefrom2to1_anglesub','absanglefrom2to1_center',...
    'absanglefrom2to1_ell2nose','absanglefrom2to1_nose2ell','absphidiff_anglesub',...
    'absphidiff_center','absphidiff_ell2nose','absphidiff_nose2ell','absthetadiff_anglesub',...
    'absthetadiff_center','absthetadiff_ell2nose','absthetadiff_nose2ell','anglesub',...
    'closestfly_anglesub','closestfly_center','closestfly_ell2nose','closestfly_nose2ell',...
    'danglesub','dcenter','ddcenter','ddell2nose','ddnose2ell','dell2nose','dnose2ell',...
    'magveldiff_anglesub','magveldiff_center','magveldiff_ell2nose','magveldiff_nose2ell',...
    'veltoward_anglesub','veltoward_center','veltoward_ell2nose','veltoward_nose2ell'};
  fprintf('Computing per-frame parameters:\n');
  fprintf('%s\n',missingfns{:});
  if any(ismember(missingfns,fns_processdata)),
    trx = process_data(trx);
  end
  if any(ismember(missingfns,fns_processdatacrabwalks)),
    trx = process_data_crabwalks(trx);
  end
  if any(ismember(missingfns,fns_processdataarena)),

    if ~isfield(trx,'arena'),
      fprintf('Enter name of annotation file to read arena position from.\n');
      annname = [matpath,strrep(matname0,'.mat','.ann')];
      if ~exist(annname,'file'),
        annnamefmf = [matpath,strrep(matname0,'.mat','.fmf.ann')];
        if exist(annnamefmf,'file')
          annname = annnamefmf;
        else
          annnamesbfmf = [matpath,strrep(matname0,'.mat','.sbfmf.ann')];
          if exist(annnamesbfmf,'file')
            annname = annnamesbfmf;
          else
            annnameavi = [matpath,strrep(matname0,'.mat','.avi.ann')];
            if exist(annnameavi,'file')
              annname = annnameavi;
            end
          end
        end
      end
      [annname,annpath] = uigetfile('*.ann',sprintf('Annotation file corresponding to %s',matname0),annname);
      annname = [annpath,annname];
      if ~exist(annname,'file'),
        fprintf('Annotation file %s does not exist, aborting\n',annname);
        return;
      else
        try
          [arena_center_x,arena_center_y,arena_radius] = ...
            read_ann(annname,'arena_center_x','arena_center_y','arena_radius');
        catch
          fprintf('Could not read annotation file %s, aborting\n',annname);
          return;
        end
        for fly = 1:length(trx),
          trx(fly).arena.x = arena_center_x;
          trx(fly).arena.y = arena_center_y;
          trx(fly).arena.r = arena_radius;
        end
      end
    end
    trx = process_data_arena(trx);
  end
  
  if any(ismember(missingfns,fns_processdataclosestfly)),
    b = inputdlg({'Field of view of fly in degrees'},'Field of view',1,{num2str(fov*180/pi)});
    if isempty(b),
      fprintf('Field of view not entered, using default value of %.1f\n',fov*180/pi);
    else
      fov = str2double(b{1})*pi/180;
    end
    trx = process_data_closestfly(trx,fov);
  end
  
end

%% name of file to save to 

fprintf('Enter name for behavior\n');
b = inputdlg({'Name of behavior to detect:'},'Behavior name',1,{''});
if isempty(b),
  return;
end
behaviorname = b{1};
fprintf('Enter results filename\n');
savename = sprintf('%sdetect_%s_%s',matpath,behaviorname,matname0);
[savename, savepath] = uiputfile('*.mat', 'Save results mat name', savename);
savename = [savepath,savename];

%% detect the behavior

nflies = length(trx);
seg = struct('t1',cell(1,nflies),'t2',cell(1,nflies),'score',cell(1,nflies),'params',cell(1,nflies));
params = struct2cell(behaviorparams);
for fly = 1:nflies,
  fprintf('Detecting %s for fly %d\n',behaviorname,fly);
  seg(fly) = systematic_detect_event4(trx(fly),behaviorname,true(1,trx(fly).nframes-1),params{:});
end

save(savename,'seg','matname','behaviorparams','behaviorname');
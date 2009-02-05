%% detect_behaviors

setuppath;

paramsmatname = '';
matname = '';
fov = pi;

pathtodetectbehaviors = which('detect_behaviors_social');
savedsettingsfile = strrep(pathtodetectbehaviors,'detect_behaviors_social.m','.detectbehaviorssocialrc.mat');
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
[paramsmatname,paramsmatpath] = uigetfile('.mat','Choose parameters mat file defining social behavior',paramsmatname);
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

%% process single-fly data, if necessary

fns = fieldnames(trx);
fns_processdata = {'dtheta','du_ctr','dv_ctr','corfrac','corisonfly','du_cor',...
    'dv_cor','velmag_ctr','velmag','accmag','signdtheta','absdv_cor','flipdv_cor',...
    'absdtheta','d2theta','absd2theta','smooththeta','smoothdtheta','abssmoothdtheta',...
    'smoothd2theta','abssmoothd2theta','phi','yaw','absyaw'};
fns_processdatacrabwalks = {'absdtheta_tail','absdu_tail','absdv_tail','absphisideways',...
    'dtheta_tail','du_tail','dv_tail','phisideways'};
% process_data, if necessary
if ~isempty(setdiff(fns_processdata,fns)),
  trx = process_data(trx,matname,moviename);  
end
if ~isempty(setdiff(fns_processdatacrabwalks,fns)),
  trx = process_data_crabwalks(trx);
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

%% main loop over flies

nflies = length(trx);
seg = struct('t1',cell(1,nflies),'t2',cell(1,nflies),'score',cell(1,nflies),'params',cell(1,nflies),...
             'otherfly',cell(1,nflies));
for fly = 1:nflies,

  % process data for this fly
  fprintf('Computing per-frame parameters for fly %d:\n',fly);
  trk = process_data_social(trx,fly);

  %% detect the behavior
  params = struct2cell(behaviorparams);
  off = 0;
  for fly2 = 1:nflies,
    if fly == fly2,
      continue;
    end
    fprintf('Detecting %s for fly1 = %d, fly2 = %d\n',behaviorname,fly,fly2);
    segcurr = systematic_detect_event4(trk(fly2),behaviorname,true(1,trk(fly2).nframes-1),params{:});
    ncurr = length(seg(fly1).t1);
    seg(fly1).t1 = [seg(fly1).t1,segcurr.t1];
    seg(fly1).t2 = [seg(fly1).t2,segcurr.t2];
    seg(fly1).otherfly = [seg(fly1).otherfly,zeros(1,ncurr)+fly2];
    seg(fly1).score = [seg(fly1).score,segcurr.score];
    off = off + ncurr;
  end
  save(savename,'seg','matname','behaviorparams','behaviorname');
end
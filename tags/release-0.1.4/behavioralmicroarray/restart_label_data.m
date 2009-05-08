% start up

setuppath;

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

%% get labeled data

labelmatname = [labelmatpath,labelmatname];
[labelmatname,labelmatpath] = uigetfile('*.mat','Save labels to...',labelmatname);
save('-append',savedsettingsfile,'labelmatname','labelmatpath');

labelmatname = [labelmatpath,labelmatname];
load(labelmatname);

%% label

fprintf('Restarting labeling for file %s.\n',labelmatname);

% where were we?
nmovies = length(movienames);
lastmoviei = length(labeledbehavior);
flieslabeled = [];
for i = 1:length(labeledbehavior{lastmoviei}),
  if ~isempty(labeledbehavior{lastmoviei}(i).starts),
    flieslabeled(end+1) = i;
  end
end
lastflyi = find(ismember(fliestolabel{lastmoviei},flieslabeled),1,'last');

for moviei = 1:nmovies,
  if moviei < lastmoviei,
    fprintf('Done labeling movie %s\n',movienames{moviei});
    continue;
  end
  load(matnames{moviei});
  trx = process_data(trx,matnames{moviei},movienames{moviei});
  for flyi = 1:length(fliestolabel{moviei}),
    if moviei == lastmoviei && flyi <= lastflyi,
      fprintf('Done labeling fly %d for movie %s\n',flyi,movienames{moviei});
      continue;
    end
    fly = fliestolabel{moviei}(flyi);
    t0 = t0tolabel{moviei}(flyi);
    t1 = t1tolabel{moviei}(flyi);
    t0 = t0 + trx(fly).firstframe - 1;
    t1 = t1 + trx(fly).firstframe - 1;
    
    trk = GetPartOfTrack(trx(fly),t0,t1);
    trk.speed = sqrt(diff(trk.x).^2 + diff(trk.y).^2);
    maxspeed = prctile(trk.speed,90);
    [labeledbehavior{moviei}(fly).starts,labeledbehavior{moviei}(fly).ends,labeledbehavior{moviei}(fly).notes] = labelbehaviors(trk,1,movienames{moviei},[],[],[],sprintf('min(%.1f,trk.speed)',maxspeed));
    save('-append',labelmatname,'labeledbehavior');
    fprintf('Done labeling fly %d of movie %s,\n(fly %d / %d for movie %d / %d)\n',fly,movienames{moviei},...
      flyi,length(fliestolabel{moviei}),moviei,nmovies);
  end
  b = questdlg('Keep labeling?','Keep labeling?','Next fly','Quit','Next fly');
  if strcmpi(b,'quit'),
    return;
  end
end

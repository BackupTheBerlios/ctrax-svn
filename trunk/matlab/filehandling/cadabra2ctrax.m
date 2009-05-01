% [succeeded,outmatname,trx] = cadabra2ctrax(featname,roiname,moviename,[outmatname])
% inputs cadabra mat file output names, loads in the trajectories, and
% converts to ctrax trx format. these trx will be saved to the file
% outmatname (user will be prompted for this name if none given). 
% inputs:
% featname: "feat" output by cadabra
% roiname: "roi" output by cadabra
% moviename: name of movie tracked
% outmatname: [optional] mat file to save trx trajectories to. user will be
% prompted for file name if none given
function [succeeded,outmatname,trx] = cadabra2ctrax(featname,roiname,moviename,outmatname)

% initialize outputs
succeeded = false;
if ~exist('outmatname','var'),
  outmatname = '';
end
trx = [];

if ~exist('featname','var'),
  helpmsg = 'Choose CADABRA feat mat file';
  [featname,featpath] = uigetfilehelp('*.mat','Choose feat file','','helpmsg',helpmsg);
  if ~ischar(featname),
    return;
  end
  featname = fullfile(featpath,featname);
end

if ~exist('roiname','var'),
  helpmsg = sprintf('Choose CADABRA roi mat file corresponding to feat file %s',featname);
  roiname = strrep(featname,'feat','roi');
  [roiname,roipath] = uigetfilehelp('*.mat','Choose feat file',roiname,'helpmsg',helpmsg);
  if ~ischar(roiname),
    return;
  end
  roiname = fullfile(roipath,roiname);
end
  
% load in data
feat = load(featname);
roi = load(roiname);

% read fps, moviename is either the name of the movie or an mmreader object
if ~ischar(moviename),
  readerobj = moviename;
  moviename = fullfile(get(readerobj,'Path'),get(readerobj,'Name'));
else
  readerobj = mmreader(moviename);
end
fps = get(readerobj,'FrameRate')+0;

% scale, pxpermm
pxpermm = 1/mean([roi.scale.x,roi.scale.y]);

% allocate
c = cell(1,2);
trx = struct('x',c,'y',c,'theta',c,'a',c,'b',c,...
  'id',c,'moviename',c,'firstframe',c,'arena',c,...
  'nframes',c,'endframe',c,'pxpermm',c,'fps',c,'x_mm',c,...
  'y_mm',c,'a_mm',c,'b_mm',c);
obj = [feat.fly_feat.obj1,feat.fly_feat.obj2];
firstframeoff = feat.fly_feat.frame(1) - 1;
arena = struct('x',nan,'y',nan,'r',nan);
for fly = 1:2,
  
  % frames for which fly is tracked
  
  % (x,y) = (0,0) for untracked frames
  badframes = obj(fly).pos_x == 0 & obj(fly).pos_y == 0;
  lastframe = find(~badframes,1,'last');
  firstframe = find(~badframes,1);
  if isempty(lastframe), 
    firstframe = 1;
    lastframe = 0;
  end
  nframes = lastframe - firstframe + 1;
  
  % allocate
  trx(fly).x = nan(1,nframes);
  trx(fly).y = nan(1,nframes);
  trx(fly).theta = nan(1,nframes);
  trx(fly).a = nan(1,nframes);
  trx(fly).b = nan(1,nframes);
  trx(fly).x_mm = nan(1,nframes);
  trx(fly).y_mm = nan(1,nframes);
  trx(fly).a_mm = nan(1,nframes);
  trx(fly).b_mm = nan(1,nframes);
  
  % store parameters
  trx(fly).moviename = moviename;
  trx(fly).firstframe = firstframe + firstframeoff;
  trx(fly).endframe = lastframe + firstframeoff;
  trx(fly).nframes = nframes;
  trx(fly).pxpermm = pxpermm;
  trx(fly).fps = fps;
  trx(fly).id = fly;
  trx(fly).arena = arena;
  %trx(fly).f2i = @(f) f - trx(fly).firstframe + 1;
  
  % store data
  idx = feat.fly_feat.frame(firstframe:lastframe)-trx(fly).firstframe + 1;
  % we use the quarter major, minor axis length
  trx(fly).a_mm(idx) = obj(fly).FLength/4;
  % maybe area is major/2 * minor/2 * pi, store quarter minor axis length
  trx(fly).b_mm(idx) = (obj(fly).FArea./(obj(fly).FLength/2)/pi)/2;
  % convert degrees to radians
  trx(fly).theta(idx) = obj(fly).headdir*pi/180;
  
  % convert mm to px, incorporate offset
  trx(fly).x(idx) = obj(fly).pos_x*pxpermm + roi.ROI.cols(1) - 1 - 1;
  trx(fly).y(idx) = obj(fly).pos_y*pxpermm + roi.ROI.rows(1) - 1 - 1;
  trx(fly).x_mm = trx(fly).x / pxpermm;
  trx(fly).x_mm = trx(fly).y / pxpermm;
  trx(fly).a = trx(fly).a_mm*pxpermm;
  trx(fly).b = trx(fly).b_mm*pxpermm;
  
end

if isempty(outmatname),
  [pathstr,name] = fileparts(featname);
  outname = strrep(name,'_feat','');
  outmatname = fullfile(pathstr,[outname,'_trx.mat']);
  helpmsg = {};
  helpmsg{1} = 'Choose mat file to save Ctrax version of trajectories loaded from:';
  helpmsg{2} = sprintf('CADABRA feat file: %s',featname);
  helpmsg{3} = sprintf('CADABRA roi file: %s',roiname);
  [outmatname,outmatpath] = uiputfilehelp('*.mat',sprintf('Save ctrax version of %s',name),outmatname,'helpmsg',helpmsg);
  if ~ischar(outmatname),
    outmatname = '';
    return;
  end
  outmatname = fullfile(outmatpath,outmatname);
  fprintf('Saving to %s...\n',outmatname);
end

save(outmatname,'trx');
succeeded = true;
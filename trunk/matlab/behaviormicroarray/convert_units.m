% calibrate tracking results: convert from pixels to mms, frames to seconds
setuppath;

%% set all defaults

pxpermm = 4;
fps = 20;
if ~exist('ISMATNAME','var'),
  ISMATNAME = false;
end
if ~exist('ISMOVIENAME','var'),
  ISMOVIENAME = false;
end
if ~exist('ISAUTOMATIC','var'),
  ISAUTOMATIC = false;
end
if ~ISMATNAME,
  matpath = '';
  matname = '';
end

%% load settings

pathtolearnparams = which('convert_units');
savedsettingsfile = strrep(pathtolearnparams,'convert_units.m','.convertunitsrc.mat');
if exist(savedsettingsfile,'file')
  if ISMATNAME,
    load(savedsettingsfile,'pxpermm','fps');
  else
    load(savedsettingsfile);
  end
end

%% get the mat file

if ~ISMATNAME,
  
fprintf('Choose mat file to convert units for.\n');
matname = [matpath,matname];
[matname,matpath] = uigetfile('*.mat','Choose mat file for which to convert units',matname);
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

end
%% read in the movie
datacurr = load(matname);
istrx = isfield(datacurr,'trx');
if ~istrx,
  if ~all(isfield(datacurr,{'x_pos','y_pos','maj_ax','min_ax','angle'})),
    msgbox(sprintf('Could not find necessary variables in %s to convert to mm, seconds, aborting',matname));
    return;
  end
end

% check to see if it has already been converted
alreadyconverted = false;
if istrx && isfield(datacurr.trx,'pxpermm') && isfield(datacurr.trx,'fps'),
  alreadyconverted = true;
end
if ~istrx && isfield(datacurr,'pxpermm') && isfield(datacurr,'fps'),
  alreadyconverted = true;
end
if alreadyconverted,
  if ~ISAUTOMATIC,
    msgbox(sprintf('Units already converted for file %s, aborting',matname));
  end
  return;
end

%% convert from frames to seconds
prompts = {'Frames per second:'};
while true,
  b = {num2str(fps)};
  b = inputdlg(prompts,'Convert frames to seconds',1,b);
  if isempty(b), return; end
  fpstmp = str2double(b{1});
  if isnan(fpstmp) || fpstmp <= 0,
    fprintf('Illegal value for frames per second -- must be a positive number.\n');
    continue;
  end
  fps = fpstmp;
  break;
end

if ~exist(savedsettingsfile,'file'),
  save(savedsettingsfile,'fps');
else
  save('-append',savedsettingsfile,'fps');
end

%% convert from px to mm
b = questdlg('Enter pixels per mm manually, or compute from landmarks','Pixels to mm','Manual','Compute','Manual');
if strcmpi(b,'manual'),
  
  % compute the mean major axis length
  if istrx,
    meanmaj = mean([datacurr.trx.a])*4;
  else
    meanmaj = mean(datacurr.maj_ax)*4;
  end
  b = {num2str(pxpermm)};
  prompts = {sprintf('Pixels per mm [mean fly length = %.1f px]',meanmaj)};

  while true,
    tmp = inputdlg(prompts,'Pixels per mm',1,b);
    if isempty(tmp),
      return;
    end
    ppmtmp = str2double(tmp{1});
    if isnan(ppmtmp) || ppmtmp <= 0,
      fprintf('Illegal value for pixels per mm -- must be a positive number.\n');
      continue;
    end
    pxpermm = ppmtmp;
    break;
  end
  
else
  
  % get moviename
  
  if ~ISMOVIENAME,
    
  fprintf('Choose the movie corresponding to the mat file %s\n',matname);
  moviename = strrep(matname,'.mat','.fmf');
  if ~exist(moviename,'file'),
    moviename = strrep(matname,'.mat','.sbfmf');
    if ~exist(moviename,'file'),
      moviename = strrep(matname,'.mat','.avi');
      if ~exist(moviename,'file'),
        moviename = matpath;
      end
    end
  end
  [moviename,moviepath] = uigetfile({'*.fmf','*.sbfmf','*.avi'}',...
    sprintf('Choose movie corresponding to %s',matname0),moviename);
  moviename = [moviepath,moviename];
  
  end
  
  % plot one image
  [readframe,nframes,fid] = get_readframe_fcn(moviename);
  im = readframe(round(nframes/2));
  figure(1); clf; imagesc(im); axis image; colormap gray; hax = gca; hold on;

  % make a draggable line
  title('Click to set endpoints of line.');
  fprintf('Draw a line on Figure 1\n');
  hline = imline(hax);
  title({'Drag around line to set landmark distance in pixels.','Double-click on line when done.'});
  position = wait(hline);
  d = sqrt(sum(diff(position).^2));
  delete(hline);
  plot(position(:,1),position(:,2),'r.-');
  b = {''};
  prompts = {sprintf('Length of this line in mm [ = %.1f px]',d)};

  while true,
    tmp = inputdlg(prompts,'Line length',1,b);
    if isempty(tmp),
      return;
    end
    mmtmp = str2double(tmp{1});
    if isnan(mmtmp) || mmtmp <= 0,
      fprintf('Illegal value for line length -- must be a positive number.\n');
      continue;
    end
    mm = mmtmp;
    break;
  end
  
  pxpermm = d / mm;
  fprintf('Pixels per mm set to %f\n',pxpermm);

  if ishandle(1), delete(1); end
end

save('-append',savedsettingsfile,'pxpermm');

%% actually do the conversion now

if istrx,
  pxfns = {'xpred','ypred','dx','dy','v'};
  % these are used for plotting, so we want to keep them in pixels
  pxcpfns = {'x','y','a','b'};
  okfns = {'x','y','theta','a','b','id','moviename','firstframe','arena',...
    'f2i','nframes','endframe','xpred','ypred','thetapred','dx','dy','v'};
  unknownfns = setdiff(fieldnames(datacurr.trx),okfns);
  if ~isempty(unknownfns),
    b = questdlg({'Do not know how to convert the following variables: ',...
      sprintf('%s, ',unknownfns{:}),'Ignore these variables and continue?'},...
      'Unknown Variables','Continue','Abort','Abort');
    if strcmpi(b,'abort'),
      return;
    end
  end
  for ii = 1:length(pxfns),
    fn = pxfns{ii};
    if isfield(datacurr.trx,fn),
      for fly = 1:length(datacurr.trx),
        datacurr.trx(fly).(fn) = datacurr.trx(fly).(fn) / pxpermm;
      end
    end
  end
  for ii = 1:length(pxcpfns),
    fn = pxcpfns{ii};
    if isfield(datacurr.trx,fn),
      for fly = 1:length(datacurr.trx),
        datacurr.trx(fly).([fn,'_mm']) = datacurr.trx(fly).(fn) / pxpermm;
      end
    end
  end

  for fly = 1:length(datacurr.trx),
    datacurr.trx(fly).pxpermm = pxpermm;
    datacurr.trx(fly).fps = fps;
  end
else
  pxfns = {'x_pos','y_pos','maj_ax','min_ax'};
  for ii = 1:length(pxfns),
    fn = pxfns{ii};
    if isfield(datacurr,fn),
      datacurr.([fn,'_mm']) = datacurr.(fn) / pxpermm;
    end
  end
  datacurr.pxpermm = pxpermm;
  datacurr.fps = fps;
end

%% save to file

savename = matname;
while true,
  [savename, savepath] = uiputfile('*.mat', sprintf('Save results for input %s to',matname), savename);
  if ~isnumeric(savename),
    break;
  end
  fprintf('missed\n');
  savename = matname;
end
savename = [savepath,savename];

fns = fieldnames(datacurr);
if strcmpi(matname,savename),
  tmpname = tempname;
  fprintf('Overwriting %s with converted data...\n',savename);
  movefile(matname,tmpname);  
  try
    save(savename,'-struct','datacurr',fns{:});
  catch
    fprintf('Aborting overwriting\n');
    movefile(tmpname,matname);
    return;
  end
  delete(tmpname);
else
  fprintf('Saving converted data to file %s.\nUse this mat file instead of %s in the future\n',savename,matname);
  save(savename,'-struct','datacurr',fns{:});
end
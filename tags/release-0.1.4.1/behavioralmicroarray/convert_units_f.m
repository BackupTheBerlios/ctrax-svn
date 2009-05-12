function [succeeded,savename] = convert_units_f(varargin)

setuppath;

fps = 20;
pxpermm = 4;
% parse inputs
[matname,matpath,moviename,ISAUTOMATIC,trx] = ...
  myparse(varargin,'matname',nan,'matpath',nan,'moviename',nan,'isautomatic',false,'trx',[]);

succeeded = false;
ISTRX = ~isempty(trx);
ISMATNAME = ischar(matname);
ISMATPATH = ischar(matpath);
if ISMATNAME && ~ISMATPATH,
  matpath = '';
  ISMATNAME = true;
end
ISMOVIENAME = ischar(moviename);
if ~ISMATNAME,
  matname = '';
  matpath = '';
end
if ~ISMOVIENAME,
  moviename = '';
end

% load settings
pathtolearnparams = which('convert_units_f');
savedsettingsfile = strrep(pathtolearnparams,'convert_units_f.m','.convertunitsrc.mat');
if exist(savedsettingsfile,'file')
  if ISMATNAME,
    load(savedsettingsfile,'pxpermm','fps');
  else
    load(savedsettingsfile);
  end
end

% get the mat file

if ISTRX,
  if ~ISMATNAME,
    if isfield(trx,'matname'),
      [matpath,matname] = split_path_and_filename(trx(1).matname);
      savename = trx(1).matname;
    else
      matname = '';
      matpath = '';
      savename = '';
    end
  end
  matname0 = matname;
else
  if ~ISMATNAME,
  
    matname = [matpath,matname];
    helpmsg = 'Choose a mat file containing trajectories for which to set the pixel to mm and frame to second conversions';
    [matname,matpath] = uigetfilehelp('*.mat','Choose mat file for which to convert units',matname,'helpmsg',helpmsg);
    if isnumeric(matname) && matname == 0,
      return;
    end
    fprintf('Matfile: %s%s\n\n',matpath,matname);
    
    if exist(savedsettingsfile,'file'),
      save('-append',savedsettingsfile,'matname','matpath');
    else
      save(savedsettingsfile,'matname','matpath');
    end
    
  end

  matname0 = matname;
  matname = [matpath,matname];
  savename = matname;
end

%% read in the movie
inputmatname = matname;
if ~ISTRX,
  [trx,matname,loadsucceeded] = load_tracks(matname);
  if ~loadsucceeded,
    msgbox(sprintf('Could not load trx mat file %s',matname));
  end
end

% check to see if it has already been converted
alreadyconverted = isfield(trx,'pxpermm') && isfield(trx,'fps');
if alreadyconverted,
  pxpermm = trx(1).pxpermm;
  fps = trx(1).fps;
  fprintf('pxpermm = %f, fps = %f\n',pxpermm,fps);
end

%% convert from frames to seconds

if ~alreadyconverted,

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
end
%% convert from px to mm

if ~alreadyconverted,

  b = questdlg('Enter pixels per mm manually, or compute from landmarks','Pixels to mm','Manual','Compute','Manual');
  if strcmpi(b,'manual'),
    
    % compute the mean major axis length
    meanmaj = mean([trx.a])*4;
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
      
      moviename = '';
      if isfield(trx,'moviename'),
        moviename = trx(1).moviename;
      end
      if isempty(moviename) || ~exist(moviename,'file'),
        moviename = strrep(inputmatname,'.mat','.fmf');
        if ~exist(moviename,'file'),
          moviename = strrep(inputmatname,'.mat','.sbfmf');
          if ~exist(moviename,'file'),
            moviename = strrep(inputmatname,'.mat','.avi');
            if ~exist(moviename,'file'),
              moviename = matpath;
            end
          end
        end
      end
      helpmsg = sprintf('Choose the movie corresponding to the mat file %s. We will use a frame from this movie to set landmark locations.',inputmatname);
      [moviename,moviepath] = uigetfilehelp({'*.fmf','*.sbfmf','*.avi'}',...
        sprintf('Choose movie corresponding to %s',matname0),moviename,'helpmsg',helpmsg);
      moviename = [moviepath,moviename];
      
    end
    
    % plot one image
    [readframe,nframes,fid] = get_readframe_fcn(moviename);
    im = readframe(round(nframes/2));
    figure(1); clf; imagesc(im); axis image; colormap gray; hax = gca; hold on;
    
    % make a draggable line
    title('Click to set endpoints of line.');
    fprintf('Draw a line on Figure 1\n');
    % allow users who don't have imline (image processing toolbox)
    if exist('imline','file'),
      try
        hline = imline(hax);
      catch
        % seems that in some releases imline does not accept just one input
        [position,hline] = get2ptline(hax);
        delete(hline);
        hline = imline(hax,position(:,1),position(:,2));
      end
      title({'Drag around line to set landmark distance in pixels.','Double-click on line when done.'});
      position = wait(hline);
    else
      title('Click on two landmarks which you know the distance between');
      [position,hline] = get2ptline(hax);
    end
    
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
end

%% actually do the conversion now

pxfns = {'xpred','ypred','dx','dy','v'};
% these are used for plotting, so we want to keep them in pixels
pxcpfns = {'x','y','a','b'};
okfns = {'x','y','theta','a','b','id','moviename','firstframe','arena',...
  'f2i','nframes','endframe','xpred','ypred','thetapred','dx','dy','v',...
  'a_mm','b_mm','x_mm','y_mm','matname','sex','type'};
unknownfns = setdiff(getperframepropnames(trx),okfns);

if ~alreadyconverted,
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
    if isfield(trx,fn),
      for fly = 1:length(trx),
        trx(fly).(fn) = trx(fly).(fn) / pxpermm;
      end
    end
  end
end

didsomething = false;
for ii = 1:length(pxcpfns),
  fn = pxcpfns{ii};
  newfn = [fn,'_mm'];
  if isfield(trx,fn) && ~isfield(trx,newfn),
    for fly = 1:length(trx),
      trx(fly).(newfn) = trx(fly).(fn) / pxpermm;
      didsomething = true;
    end
  end
end

if ~alreadyconverted,
  for fly = 1:length(trx),
    trx(fly).pxpermm = pxpermm;
    trx(fly).fps = fps;
  end
end

%% save to file

if ~alreadyconverted || didsomething,

  nmissed = 0;
  while true,
    helpmsg = sprintf('Choose the file to which to save the trx from %s augmented with the pixel to mm and frame to second conversions',inputmatname);
    [savename, savepath] = uiputfilehelp('*.mat', sprintf('Save results for input %s to',matname0), savename,'helpmsg',helpmsg);
    if ~isnumeric(savename),
      break;
    end
    fprintf('missed\n');
    nmissed = nmissed + 1;
    if nmissed > 1,
      return;
    end
    savename = inputmatname;
  end
  savename = [savepath,savename];
  
  for i = 1:length(trx),
    trx(i).matname = savename;
  end
  
  if strcmpi(matname,savename),
    tmpname = tempname;
    fprintf('Overwriting %s with converted data...\n',savename);
    movefile(matname,tmpname);
    try
      save('-append',savename,'trx');
    catch
      fprintf('Aborting overwriting\n');
      movefile(tmpname,matname);
      return;
    end
    delete(tmpname);
  else
    fprintf('Saving converted data to file %s.\nUse this mat file instead of %s in the future\n',savename,inputmatname);
    copyfile(matname,savename);
    save('-append',savename,'trx');
  end
end
succeeded = true;
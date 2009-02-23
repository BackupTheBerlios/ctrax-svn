% script that prompts user for mat, annotation, and movie files, parameters
% for computing suspicious frames, then computes suspicious frames, then
% brings up the fixerrors gui

%% set all defaults

moviename = '';
moviepath = '';
setuppath;

%% read last settings
pathtofixerrors = which('fixerrors');
savedsettingsfile = strrep(pathtofixerrors,'fixerrors.m','.fixerrorsrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% choose movie

fprintf('Choose a movie to fix errors in\n');
movieexts = {'*.fmf','*.sbfmf','*.avi'}';
[moviename,moviepath,filterindex] = uigetfile(movieexts,'Choose movie file',moviename);
if isnumeric(moviename) && moviename == 0, 
  return;
end

fprintf('Choose the corresponding mat file containing the trajectories.\n');
matname = [moviepath,strrep(moviename,movieexts{filterindex}(2:end),'.mat')];
[matname,matpath] = uigetfile({'.mat'},'Choose mat file',matname);
if isnumeric(matname) && matname == 0, 
  return;
end

annname = [matpath,moviename,'.ann'];
fprintf('Choose the corresponding annotation file containing the trajectories.\n');
[annname,annpath] = uigetfile({'.ann'},'Choose ann file',annname);
if isnumeric(annname) && annname == 0,
  return;
end

moviename = [moviepath,moviename];
matname = [matpath,matname];
annname = [annpath,annname];

[readframe,nframes,fid] = get_readframe_fcn(moviename);

if exist('savedsettingsfile','file'),
  save('-append',savedsettingsfile,'moviename','moviepath');
else
  save(savedsettingsfile,'moviename','moviepath');
end

%% convert to px, seconds

ISAUTOMATIC = true;
ISMATNAME = true;
ISMOVIENAME= true;
savedsettingsfile0 = savedsettingsfile;
convert_units;
if ~alreadyconverted,
  if ~ischar(savename),
    error('Savename not set in convert_units');
    return;
  end
  matname = savename;
end
savedsettingsfile = savedsettingsfile0;

%% see if we should restart

tag = moviename(length(moviepath)+1:end-length(movieexts{filterindex})+1);
loadname = sprintf('tmpfixed_%s.mat',tag);
DORESTART = false;
  
if exist(loadname,'file'),
  tmp = dir(loadname);
  prompt = ['An old file saved by fixerrors was found for movie ',...
    sprintf('%s.%s, last modified %s. ',...
    tag,movieexts{filterindex}(2:end),tmp(1).date),...
    'Would you like to load the saved results and restart?'];
  button = questdlg(prompt,'Restart?');
  if strcmpi(button,'yes'),
    DORESTART = true;
  elseif strcmpi(button,'cancel'),
    return
  end
  
end

%% set parameters for detecting suspicious frames

if ~DORESTART,

[max_jump,maxmajor,meanmajor,arena_radius] = ...
  read_ann(annname,'max_jump','maxmajor','meanmajor','arena_radius');
meanmajor = meanmajor * 4;
maxmajor = maxmajor * 4;

px2mm = pxpermm;
save('-append',savedsettingsfile,'px2mm');

max_jump = max_jump / px2mm;
maxmajor = maxmajor / px2mm;
meanmajor = meanmajor / px2mm;

% set default values
if ~exist('minerrjump','var')
  minerrjump = .2*max_jump;
end
if ~exist('minorientchange','var'),
  minorientchange = 45;
end
if ~exist('largemajor','var'),
  largemajor = meanmajor + 2/3*(maxmajor-meanmajor);
end
if ~exist('minanglediff','var'),
  minanglediff = 90;
end
if ~exist('minwalkvel','var'),
  minwalkvel = 1 / 4;
end
if ~exist('matcherrclose','var'),
  matcherrclose = 10/4^2;
end
tmp = [minerrjump,minorientchange,largemajor,minanglediff,minwalkvel,matcherrclose];
defaultv = cell(size(tmp));
for i = 1:length(tmp),
  defaultv{i} = num2str(tmp(i));
end

shortdescr = cell(1,6);
descr = cell(1,6);
relev = cell(1,6);
shortdescr{1} = 'Minimum suspicious prediction-detection error (mm)';
descr{1} = ['All sequences in which the error between the constant velocity ',...
  'prediction and measured position is greater than the given value ',...
  'will be flagged. '];
relev{1} = sprintf('Max jump error: %.1f (mm)',max_jump);
shortdescr{2} = 'Minimum suspicious orientation change (deg)';
descr{2} = ['All sequences in which the change in orientation is greater ',...
  'than the given value will be flagged.'];
relev{2} = '';
shortdescr{3} = 'Minimum suspiciously large major axis (mm)';
descr{3} = ['All sequences in which the major axis length is greater than ',...
  'the given value will be flagged.'];
relev{3} = sprintf('Mean major axis length: %.2f mm, max major axis length: %.2f mm',meanmajor,maxmajor);
shortdescr{4} = 'Minimum suspicious orientation-velocity direction mismatch (deg): ';
descr{4} = '';
relev{4} = '';
shortdescr{5} = 'Minimum walking speed (mm/frame)';
descr{5} = ['All sequences in which the fly is walking (has speed greater than ',...
  'the given value) and the orientation and velocity differ by the given ',...
  'value will be flagged.'];
relev{5} = '';
shortdescr{6} = 'Maximum ambiguous error (mm^2)';
descr{6} = ['All sequences in which the increase in error for swapping ',...
  'a pair of identities is less than the given value will be flagged.'];
relev{6} = '';
prompts = cell(size(shortdescr));
for i = 1:length(shortdescr),
  prompts{i} = sprintf('**%s**: %s. ',shortdescr{i},descr{i});
  if ~isempty(relev{i}),
    prompts{i} = [prompts{i},sprintf('. [Relevant quantities: %s]',relev{i})];
  end
end
title1 = 'Suspiciousness Parameters';
tmp = inputdlg(prompts,title1,1,defaultv,'on');

if isempty(tmp),
  return;
end

minerrjump = str2double(tmp{1});
minorientchange = str2double(tmp{2});
largemajor = str2double(tmp{3});
minanglediff = str2double(tmp{4});
minwalkvel = str2double(tmp{5});
matcherrclose = str2double(tmp{6});

save('-append',savedsettingsfile,...
  'minerrjump','minorientchange','largemajor','minanglediff',...
  'minwalkvel','matcherrclose');

end

%% convert to the units expected by suspicious_sequences

if ~DORESTART,

minerrjumpfrac = minerrjump / max_jump;
minorientchange = minorientchange*pi/180;
maxmajorfrac = (largemajor - meanmajor)/(maxmajor - meanmajor);
minwalkvel = minwalkvel*px2mm;
matcherrclose = matcherrclose*px2mm^2;
minanglediff = minanglediff*pi/180;
[seqs,trx0,params] = suspicious_sequences(matname,annname,...
  'minerrjumpfrac',minerrjumpfrac,'minorientchange',minorientchange,...
  'maxmajorfrac',maxmajorfrac,'minwalkvel',minwalkvel,...
  'matcherrclose',matcherrclose,'minanglediff',minanglediff);

end

%% call the fixerrors gui

fprintf('Movie: %s\n',moviename);
fprintf('Mat: %s\n',matname);
fprintf('Annname: %s\n',annname);
fprintf('Temporary file created at: %s\n',loadname);

if ~DORESTART,
  trx = fixerrorsgui(seqs,moviename,trx0,annname,params,loadname);
else
  load(loadname);
  trx0 = trx;
  trx = fixerrorsgui(seqs,moviename,trx0,annname,params,loadname);
end

%% save

for i = 1:length(trx),
  trx(i).pxpermm = pxpermm;
  trx(i).fps = fps;
end

fprintf('Save fixed tracks to a file\n');
savename = sprintf('fixed_%s.mat',tag);
[savename, savepath] = uiputfile('*.mat', 'Save results?', savename);
if isnumeric(savename) && savename == 0,
  return;
end
trx = rmfield(trx,{'xpred','ypred','thetapred','dx','dy','v','f2i'});
if ~isempty(savename),
  save(savename,'trx');
else
  trx = rmfield(trx,{'xpred','ypred','thetapred','dx','dy','v','f2i'});
  tmpsavename = sprintf('backupfixed_movie%s.mat',tag);
  save(tmpsavename,'trx');
  msgbox(sprintf('saving trx to file %s\n',tmpsavename));
end

%% histogram per-frame stats
setuppath;

%% set all defaults

matname = '';

%% load settings

pathtohistogramproperties = which('histogramproperties');
savedsettingsfile = strrep(pathtohistogramproperties,'histogramproperties.m','.histogrampropertiesscriptrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% choose a mat file to analyze
fprintf('Choose per-frame stats mat file to plot.\n');
[matname,matpath] = uigetfile('*.mat','Choose mat file to analyze',matname);
if isnumeric(matname) && matname == 0,
  return;
end
matname = [matpath,matname];
fprintf('Matfile: %s\n\n',matname);

if exist(savedsettingsfile,'file'),
  save('-append',savedsettingsfile,'matname','matpath');
else
  save(savedsettingsfile,'matname','matpath');
end

datacurr = load(matname);
if ~isfield(datacurr,'trx') || ~isfield(datacurr.trx,'units'),
  fprintf('No variable trx or trx.units required for histogramming\n');
  return;
end

fprintf('Starting GUI histogramproperties, manipulate plot using histogramproperties figure\n');
histogramproperties(datacurr.trx);
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
[trx,matname,loadsucceeded] = load_tracks();
if ~loadsucceeded,
  return;
end
fprintf('Matfile: %s\n\n',matname);

if exist(savedsettingsfile,'file'),
  save('-append',savedsettingsfile,'matname','matpath');
else
  save(savedsettingsfile,'matname','matpath');
end

fprintf('Starting GUI histogramproperties, manipulate plot using histogramproperties figure\n');
histogramproperties(trx);
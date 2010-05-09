function succeeded = learn_area_classifier()

succeeded = false;

% load defaults
matname = '';
donormalize = true;
normalizematname = '';

pathtolearnareaclassifier = which('learn_area_classifier');
savedsettingsfile = strrep(pathtoclassifybyarea,'learn_area_classifier.m','.learnareaclassifierrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end


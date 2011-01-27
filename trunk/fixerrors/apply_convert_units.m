function trx = apply_convert_units(trx,pxpermm,fps)

pxpermminput = exist('pxpermm','var');
fpsinput = exist('fps','var');

if ~pxpermminput
  pxpermm = trx(1).pxpermm;
end
if ~fpsinput
  fps = trx(1).fps;
end

%% actually do the conversion now

pxfns = {'xpred','ypred','dx','dy','v'};
% these are used for plotting, so we want to keep them in pixels
pxcpfns = {'x','y','a','b'};
okfns = {'x','y','theta','a','b','id','moviename','firstframe','arena',...
  'f2i','nframes','endframe','xpred','ypred','thetapred','dx','dy','v',...
  'a_mm','b_mm','x_mm','y_mm','matname','sex','type','timestamps'};
unknownfns = setdiff(getperframepropnames(trx),okfns);

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

for ii = 1:length(pxcpfns),
  fn = pxcpfns{ii};
  newfn = [fn,'_mm'];
  if isfield(trx,fn),
    for fly = 1:length(trx),
      trx(fly).(newfn) = trx(fly).(fn) / pxpermm;
    end
  end
end

for fly = 1:length(trx),
  if pxpermminput,
    trx(fly).pxpermm = pxpermm;
  end
  if fpsinput,
    trx(fly).fps = fps;
  end
end

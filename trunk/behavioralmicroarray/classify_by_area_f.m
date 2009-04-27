function [succeeded,savematnames,classifiermatname,areathresh] = classify_by_area(varargin)

setuppath;
succeeded = false;
savematnames = {};
classifiermatname = '';
areathresh = nan;

% parse inputs
[matnames,donormalize,normalizematname,classifiermatname] = ...
  myparse(varargin,'matnames',{},'donormalize',nan,'normalizematname','',...
  'classifiermatname','');
ismatnames = ~isempty(matnames);
isdonormalize = ~isnan(donormalize);
isnormalizematname = ~isempty(normalizematname);
isclassifiermatname = ~isempty(classifiermatname);

% load defaults
pathtoclassifybyarea = which('classify_by_area');
savedsettingsfile = strrep(pathtoclassifybyarea,'classify_by_area.m','.classifybyarearc.mat');
if exist(savedsettingsfile,'file')
  defaults = load(savedsettingsfile);
  if ~isdonormalize && isfield(defaults,'donormalize'),
    donormalize = defaults.donormalize;
  end
  if ~ismatnames && isfield(defaults,'matname'),
    matnames = {defaults.matname};
  end
  if ~isnormalizematname && isfield(defaults,'normalizematname'),
    normalizematname = defaults.normalizematname;
  end
  if ~isclassifiermatname && isfield(defaults,'classifiermatname'),
    classifiermatname = defaults.classifiermatname;
  end
else
  if ~ismatnames,
    matnames = {''};
  end
  if ~isdonormalize,
    donormalize = false;
  end
end

% get trajectory file names
if ~ismatnames,
  fprintf('Choose trajectory mat file(s) for which to classify flies based on area.\n');
  [matname0s,matpath] = uigetfile('*.mat','Choose trajectory mat file(s)',matnames{end},'multiselect','on');
  if iscell(matname0s),
    matnames = cell(size(matname0s));
    for i = 1:length(matname0s),
      matnames{i} = [matpath,matname0s{i}];
    end
  else
    if ~ischar(matname0s),
      return;
    end
    matnames = {[matpath,matname0s]};
  end
  matname = matnames{end};
end

if isclassifiermatname,
  doloadclassifier = true;
else
  b = questdlg('Load in an existing area-based classifier, or create a new classifier?',...
    'Load classifier?','Load','Create','Cancel','Create');
  if strcmpi(b,'cancel'),
    return;
  end
  doloadclassifier = strcmpi(b,'load');
end

if doloadclassifier,

  while true,

    fprintf('Choose mat file containing area-based classifier\n');
    [classifiermatname,classifiermatpath] = uigetfile('*.mat','Choose classifier mat file',classifiermatname);

    if ~ischar(classifiermatname),
      return;
    end

    classifiermatname = [classifiermatpath,classifiermatname];
    tmp = load(classifiermatname);
    if ~all(isfield(tmp,{'areathresh','area_coeffs','normalize_area','typename','type0','type1'})),
      msgbox(sprintf('Mat file %s does not contain the correct variables',classifiermatname));
      continue;
    end

    areathresh = tmp.areathresh;
    area_coeffs = tmp.area_coeffs;
    normalize_area = tmp.normalize_area;
    typename = tmp.typename;
    type0 = tmp.type0;
    type1 = tmp.type1;

    break;

  end

else

  % will we normalize by location in the image?
  if ~isdonormalize,
    if donormalize,
      defaults = 'Yes';
    else
      defaults = 'No';
    end
    b = questdlg('Normalize the area to account for differences in lighting in different parts of the arena?',...
      'Normalize area?','Yes','No','Cancel',defaults);
    if strcmpi(b,'cancel'),
      return;
    end
    donormalize = strcmpi(b,'yes');
  end

  % get normalization data
  if donormalize,

    if isnormalizematname,
      doload = true;
    else
      b = questdlg('Load in normalization terms, or calculate normalization terms from trajectories of flies of one type?',...
        'Load in normalization terms?','Load','Calculate','Cancel','Calculate');
      if strcmpi(b,'cancel'),
        return;
      end
      doload = strcmpi(b,'load');
    end

    % load normalization data from a mat file
    if doload,
      while true,
        [normalizematname0,normalizematpath] = uigetfile('*.mat','Load Normalization Terms',normalizematname);
        if ~ischar(normalizematname0),
          return;
        end
        normalizematname1 = [normalizematpath,normalizematname0];
        if ~exist(normalizematname1,'file'),
          msgbox(sprintf('File %s does not exist',normalizematname1));
          continue;
        end
        normalizationdata = load(normalizematname1);
        if ~isfield(normalizationdata,'area_coeffs'),
          msgbox(sprintf('Invalid normalization file: %s does not contain the variable area_coeffs',normalizematname1));
          continue;
        end
        normalizematname = normalizematname1;
        break;
      end
      normalizematname = normalizematname1;
      area_coeffs = normalizationdata.area_coeffs;
      minx = normalizationdata.minx;
      maxx = normalizationdata.maxx;
      miny = normalizationdata.miny;
      maxy = normalizationdata.maxy;
      hmatnames = setdiff(normalizationdata.matnamesnorm,matnames);
    else
      % compute normalization terms
      
      % input more trajectories, if desired
      b = questdlg('Do you want to enter other mat files of trajectories from which to learn the normalization function?');
      if strcmpi(b,'cancel'),
        return;
      end
      if strcmpi(b,'yes'),
        fprintf('Choose mat file(s) of trajectories for use in normalization.\n');
        % load in homogeneious matnames
        [hmatnames0,hmatpath] = uigetfile('*.mat','Choose trajectory mat file(s)',matname,'multiselect','on');
        if iscell(hmatnames0),
          hmatnames = cell(size(hmatnames0));
          for i = 1:length(hmatnames0),
            hmatnames{i} = [hmatpath,hmatnames0{i}];
          end
        else
          if ~ischar(hmatnames0),
            return;
          end
          hmatnames = {[hmatpath,hmatnames0]};
        end
        b = questdlg('In the extra mat files to use for normalization, are all the flies described by a single mat file of the same type (e.g. all male wild-type in one mat file, all female wild-type in another)');
        if strcmpi(b,'cancel'),
          return;
        end
        ishomogeneous = strcmpi(b,'yes');
      else
        ishomogeneous = false;
        hmatnames = {};
      end
      
      if ishomogeneous,
        matnamesnorm = hmatnames;
      else
        matnamesnorm = [hmatnames,matnames];
      end
      
      % collect x, y, and area for all frames and all flies
      nmovies = length(matnamesnorm);
      X = [];
      Y = [];
      AREA = [];
      for i = 1:nmovies,
        fprintf('Loading movie %s to compute normalization function.\n',matnamesnorm{i});
        [trx,matnamesnorm{i}] = load_tracks(matnamesnorm{i});
        X = [X,[trx.x]];
        Y = [Y,[trx.y]];
        if ishomogeneous,
          areacurr = [trx.a].*[trx.b];
          medianarea = median(areacurr);
          AREA = [AREA,areacurr/medianarea];
        else
          nflies = length(trx);
          for fly = 1:nflies,
            areacurr = trx(fly).a.*trx(fly).b;
            medianarea = median(areacurr);
            AREA = [AREA,areacurr/medianarea];
          end
        end
      end
      fprintf('Performing regression to compute normalization function\n');
      ndata = length(X);
      in = [X.^2;Y.^2;X;Y;ones(1,ndata)];
      area_coeffs = regress(AREA',in');
      
      minx = min(X);
      maxx = max(X);
      miny = min(Y);
      maxy = max(Y);
      
      fprintf('Choose file to save area normalization function to.\n');
      [savenormname,savenormpath] = uiputfile('*.mat','Save area normalization function','');
      if ischar(savenormname),
        savenormname = [savenormpath,savenormname];
        save(savenormname,'area_coeffs','minx','maxx','miny','maxy','matnamesnorm');
      end
      
      hmatnames = matnamesnorm(1:length(hmatnames));
      if ~ishomogeneous,
        matnames = matnamesnorm(length(hmatnames)+1:end);
      end
      
    end % doload
    
    % compute normalization functions
    predict_areafactor = @(x,y) reshape([x(:).^2,y(:).^2,x(:),y(:),ones(length(x(:)),1)]*area_coeffs,size(x));
    normalize_area = @(x,y,a) a ./ predict_areafactor(x,y);
    
    [xgrid,ygrid] = meshgrid(linspace(minx,maxx,50),linspace(miny,maxy,50));
    areafactorpredicted = predict_areafactor(xgrid,ygrid);

    figure;
    imagesc([minx,maxx],[miny,maxy],areafactorpredicted);
    colorbar;
    colormap jet;
    title('Predicted multiple of median area based on position in the arena');
    xlabel('x-position (pixels)');
    ylabel('y-position (pixels)');
    
  else

    hmatnames = {};
    area_coeffs = zeros(5,1); area_coeffs(end) = 1;
    normalize_area = @(x,y,a) a;
    
  end

  matnamesall = [matnames,hmatnames];
  nmovies = length(matnamesall);
  
  area = [];
  moviei = [];
  prcts = [5,25];
  prcts = [prcts,50,100-fliplr(prcts)];
  medprct = ceil(length(prcts)/2);
  areaprcts = zeros(length(prcts),0);

  for i = 1:nmovies,
    fprintf('Processing movie %s\n',matnamesall{i});
    [trx,matnamesall{i}] = load_tracks(matnamesall{i});
    nflies = length(trx);
    moviei(end+1:end+nflies) = i;
    for fly = 1:nflies,
      areacurr = normalize_area(trx(fly).x,trx(fly).y,trx(fly).a.*trx(fly).b*4*pi);
      areaprctscurr = prctile(areacurr,prcts);
      areaprcts(:,end+1) = areaprctscurr;
    end
  end
  area = areaprcts(medprct,:);
  [sortedarea,order] = sort(area);
  sortedmoviei = moviei(order);
  sortedareaprcts = areaprcts(:,order);

  meanareapermovie = zeros(1,nmovies);
  for i = 1:nmovies,
    meanareapermovie(i) = mean(area(moviei==i));
  end
  
  [tmp,movieorder] = sort(meanareapermovie);
  [tmp,moviecolororder] = sort(movieorder);
  
  fig = figure;
  clf;
  hold on;
  hax = gca;

  hprcts = [];
  colorprcts = 1-repmat((.5+abs(50-prcts)'/50)/2,[1,3]);
  colorsprcts = flipud(gray(ceil(length(prcts)/2)));
  for i = 1:length(prcts),
    if i == medprct, continue; end
    hprcts(end+1) = plot(sortedareaprcts(i,:),'-','color',colorprcts(i,:));
    text(length(area)+1,sortedareaprcts(i,end),sprintf('%d%%',prcts(i)),'horizontalalignment','left','color','w');
  end

  h = zeros(1,nmovies);
  % sort movies by 
  colors = jet(nmovies);
  colors = colors(moviecolororder,:);

  for i = nmovies:-1:1,
    idx = sortedmoviei == i;
    h(i) = plot(find(idx),sortedarea(idx),'o','color',colors(i,:),'markerfacecolor',colors(i,:));
  end
  xlabel('Fly identity');
  ylabel('Area (px^2)');
  legends = cell(1,nmovies);
  for i = 1:nmovies,
    [spath,sname] = split_path_and_filename(matnamesall{i});
    sname = strrep(sname,'_','\_');
    if i > length(matnames),
      legends{i} = sprintf('Median area %s',sname);
    else
      legends{i} = sprintf('Median area %s',sname);
    end
  end
  legend(h(movieorder),legends(movieorder));
  axisalmosttight;
  set(hax,'color','k');

  if exist('imline','file'),
  
    % cluster the areas into two groups to initialize threshold
    [mu,cost,idx,areathresh] = onedimkmeans(sortedarea,2,'issorted',true);
    areathresh = areathresh(2);
    flythresh = find(idx ~= idx(1),1,'first');
    flythresh = flythresh - .5;

    ax = axis;
    threshy = ax(3:4)+100*[-1,1];
    threshx = ax(1:2)+100*[-1,1];
    hthresh = imline(hax,flythresh+[0,0],threshy);
    setColor(hthresh,'r');
    hareathresh = imline(hax,threshx,areathresh+[0,0]);
    setColor(hareathresh,'r');
    setPositionConstraintFcn(hthresh,@(x) [mean(x(1:2,1))+[0;0],threshy']);
    setPositionConstraintFcn(hareathresh,@(x) [threshx',mean(x(1:2,2))+[0;0]]);
    addNewPositionCallback(hthresh,@(x) hthreshCallback(x,sortedarea,hareathresh,threshx));
    addNewPositionCallback(hareathresh,@(x) hareathreshCallback(x,sortedarea,hthresh,threshy));
    
    realclosereq = get(fig,'closeRequestFcn');
    set(fig,'closerequestfcn','');
    input('Drag around the red lines to set the area threshold. Hit ENTER when done: ');
    
    x = getPosition(hareathresh);
    areathresh = mean(x(:,2));
    
    set(fig,'closerequestfcn',realclosereq);

  else
    
    % no imline function
    ax = axis;
    tmpx = [ax(1),1:length(sortedarea),ax(2)];
    if length(sortedarea) > 1,
      tmpy = [2*sortedarea(1) - sortedarea(2), sortedarea, 2*sortedarea(end)-sortedarea(end-1)];
    else
      tmpy = [sortedarea(1)-1,sortedarea,sortedarea(end)+1];
    end
    plot(tmpx,tmpy,'r-');
    title('Click on the red line to set the area threshold. Note that only the y-position of the clicked point is used.');
    while true,
      [x,areathresh] = ginput(1);
      hthresh = plot(ax(1:2),areathresh+[0,0],'m-','linewidth',2);
      b = questdlg('Use area threshold shown in magenta?','Okay?','Use','Redo','Use');
      if strcmpi(b,'use'),
        break;
      end
      delete(hthresh);
    end
  end
    
  while true,
    v = inputdlg({'Name of fly type classification (e.g. "sex")',...
      'Name of smaller type (e.g. "M" for male)',...
      'Name of larger type (e.g. "F" for female)'},...
      'Classification names',1,...
      {'sex','M','F'});
    if isempty(v),
      return;
    end
    typename = v{1};
    type0 = v{2};
    type1 = v{3};
    if ~isvarname(typename),
      msgbox('Classification type must be a valid variable name, e.g. sex');
      continue;
    end
    break;
  end

end % end doloadclassifier

% classify the flies
savematnames = cell(size(matnames));
[savematpath,savematname] = split_path_and_filename(matnames{1});
for i = 1:length(matnames),
  fprintf('Classifying area for movie %s\n',matnames{i});
  [trx,matnames{i}] = load_tracks(matnames{i});
  nflies = length(trx);
  for fly = 1:nflies,
    area = median(normalize_area(trx(fly).x,trx(fly).y,trx(fly).a.*trx(fly).b*4*pi));
    if area <= areathresh,
      trx(fly).(typename) = type0;
    else
      trx(fly).(typename) = type1;
    end
  end

  % save the results
  fprintf('Choose output mat file for classified trajectories\n');
  [tmp,savematname] = split_path_and_filename(matnames{i});
  [savematnames{i},savematpath] = uiputfile('*.mat',...
    sprintf('Save classified trajectories for %s',savematname),...
    [savematpath,savematname]);
  if ~ischar(savematnames{i}),
    savematnames{i} = '';
    continue;
  end
  savematnames{i} = [savematpath,savematnames{i}];
  if strcmp(savematnames{i},matnames{i}),
    save('-append',savematnames{i},'trx');
  else
    save(savematnames{i},'trx');
  end
end

if ~doloadclassifier,
  % save the classifier
  fprintf('Choose mat file to save classifier parameters to.\n');
  [classifiermatname,savepath] = uiputfile('*.mat','Save area-based classifier');
  if ischar(classifiermatname),
    classifiermatname = [savepath,classifiermatname];
    save(classifiermatname,'areathresh','area_coeffs','normalize_area','typename','type0','type1');
  end
end

% save defaults
matname = matnames{end};
save(savedsettingsfile,'donormalize','matname','normalizematname','classifiermatname');

succeeded = true;

function hthreshCallback(x,sortedarea,hareathresh,threshx)

flythresh = mean(x(1:2,1));
if flythresh <= 1,
  areathresh = sortedarea(1);
elseif flythresh >= length(sortedarea),
  areathresh = sortedarea(end);
else
  areathresh = (flythresh - floor(flythresh))*sortedarea(ceil(flythresh)) + ...
    (ceil(flythresh) - flythresh)*sortedarea(floor(flythresh));
end
setPosition(hareathresh,threshx,areathresh+[0,0]);


function hareathreshCallback(x,sortedarea,hthresh,threshy)

areathresh = mean(x(1:2,2));
flythresh = find(sortedarea>areathresh,1,'first');
if isempty(flythresh),
  flythresh = length(sortedarea);
elseif flythresh > 1,
  off1 = sortedarea(flythresh) - areathresh;
  Z = sortedarea(flythresh) - sortedarea(flythresh-1);
  flythresh = flythresh - off1/Z;
end
setPosition(hthresh,flythresh+[0,0],threshy);
function succeeded = classify_by_area(varargin)

succeeded = false;

% parse inputs
[matnames,donormalize,normalizematname] = ...
  myparse(varargin,'matnames',{},'donormalize',nan,'normalizematname','');
ismatnames = ~isempty(matnames);
isdonormalize = ~isnan(donormalize);
isnormalizematname = ~isempty(normalizematname);

% load defaults
pathtoclassifybyarea = which('classify_by_area');
savedsettingsfile = strrep(pathtoclassifybyarea,'classify_by_area.m','.classifybyarearc.mat');
if exist(savedsettingsfile,'file')
  defaults = load(savedsettingsfile);
  if ~isdonormalize,
    donormalize = defaults.donormalize;
  end
  if ~ismatnames,
    matnames = {defaults.matname};
  end
  if ~isnormalizematname,
    normalizematname = defaults.normalizematname;
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
  [matpath,matname0s] = uigetfile('*.mat','Choose trajectory mat file(s)',matnames{end},'multiselect','on');
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
      if isfield(normalizationdata,'area_coeffs'),
        msgbox(sprintf('Invalid normalization file: %s does not contain the variable area_coeffs',normalizematname1));
        continue;
      end
      break;
    end
    normalizematname = normalizematname1;
    area_coeffs = normalizationdata.area_coeffs;
    minx = normalizationdata.minx;
    maxx = normalizationdata.maxx;
    miny = normalizationdata.miny;
    maxy = normalizationdata.maxy;
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
      b = questdlg('Are all the flies described by a single mat file of the same type (e.g. all male wild-type in one mat file, all female wild-type in another)');
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
  
  normalize_area = @(x,y,a) a;

end

matnamesall = [matnames,hmatnames];
nmovies = length(matnamesall);

area = [];
moviei = [];
prcts = [1,25,50,75,99];
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

figure;
if nmovies <= 7,
  colors = lines(nmovies);
else
  colors = [lines(7);jet(nmovies-7)];
end

h = zeros(1,nmovies);

for i = nmovies:-1:1,
  idx = sortedmoviei == i;
  h(i) = plot(find(idx),sortedarea(idx),'o','color',colors(i,:),'markerfacecolor',colors(i,:));
  hold on;
end
xlabel('Fly identity');
ylabel('Median area');
legends = cell(1,nmovies);
for i = 1:nmovies,
  [spath,sname] = split_path_and_filename(matnamesall{i});
  sname = strrep(sname,'_','\_');
  if i > length(matnames),
    legends{i} = sprintf('Extra movie %s',sname);
  else
    legends{i} = sprintf('Movie %s',sname);
  end
end
legend(h,legends);
axisalmosttight;
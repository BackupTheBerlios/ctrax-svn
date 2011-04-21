function [succeeded,figs] = plot_detectbehaviors(trx,seg,issocial)

succeeded = false;
figs = [];

if ~exist('trx','var'),
  helpmsg = 'Choose trx file to load in';
  [matname,matpath] = uigetfilehelp('*.mat','Choose trx mat file to load in','','helpmsg',helpmsg);
  if ~ischar(matname),
    return;
  end
  matname = [matpath,matname];
  [trx,loadmatname,loadsucceeded] = load_tracks(matname);
  if ~loadsucceeded,
    msgbox(sprintf('Could not load trx from %s',matname));
    return;
  end
else
  if isfield(trx,'matname'),
    matname = trx(1).matname;
    [matpath,tmp] = split_path_and_filename(matname);
  else
    matname = '';
    matpath = '';
  end
end

nflies = length(trx);

if ~exist('seg','var'),
  helpmsg = sprintf('Choose seg file corresponding to trx file %s\n',matname);
  [segname,segpath] = uigetfilehelp('*.mat','Choose trx mat file to load in',matpath,'helpmsg',helpmsg);
  if ~ischar(segname),
    return;
  end
  segname = [segpath,segname];
  tmp = load(segname);
  if ~isfield(tmp,'seg'),
    msgbox(sprintf('seg is not defined in %s',segname));
    return;
  end
  if isfield(tmp,'trxname') && ~strcmp(matname,tmp.trxname),
    fprintf('seg in %s computed for matfile %s, but loaded in file %s\n',segname,tmp.trxname,matname);
  end
  seg = tmp.seg;
  if length(seg) ~= nflies,
    msgbox(sprintf('Number of flies in %s does not match %s',segname,matname));
    return;
  end
end

if ~exist('issocial','var'),
  issocial = isfield(trx,'fly1');
end
  
if ~isfield(trx,'x_mm'),
  for fly = 1:nflies,
    trx(fly).x_mm = trx(fly).x;
    trx(fly).y_mm = trx(fly).y;
  end
end   

% show at most maxn1 x maxn2 plots per figure
maxn1 = 3;
maxn2 = 4;
nperpage = maxn1*maxn2;
npages = ceil(nflies/nperpage);
nlastpage = mod(nflies,nperpage);
if nlastpage == 0, nlastpage = nperpage; end

page = 1;
plotcurr = 1;
if npages == 1,
  nperpagecurr = nlastpage;
else
  nperpagecurr = nperpage;
end

for fly = 1:nflies,

  % go to the next page
  if plotcurr == 1,
    figs(end+1) = figure('position',[360   402   886   519],'units','pixels');
    if page == npages,
      n2 = ceil(sqrt(nlastpage));
      n1 = ceil(nlastpage / n2);
    else
      n1 = maxn1;
      n2 = maxn2;
    end
    hax = createsubplots(n1,n2,.01);
  end
  
  axes(hax(plotcurr));
  hold on;
  
  % plot detections
  for i = 1:length(seg(fly).t1),
    hlabel = plot(trx(fly).x_mm(seg(fly).t1(i):seg(fly).t2(i)),...
      trx(fly).y_mm(seg(fly).t1(i):seg(fly).t2(i)),...
      'm-','linewidth',8);
    plot(trx(fly).x_mm(seg(fly).t1(i)),...
      trx(fly).y_mm(seg(fly).t1(i)),'go','markerfacecolor','g');
    plot(trx(fly).x_mm(seg(fly).t2(i)),...
      trx(fly).y_mm(seg(fly).t2(i)),'ro','markerfacecolor','r');
  end
  
  % plot the other fly's position
  if issocial,
    % get frames in and around detections
    idxdetect = false(1,trx(fly).nframes);
    for i = 1:length(seg(fly).t1),
      idxdetect(seg(fly).t1(i):seg(fly).t2(i)) = true;
    end
    idxplotother = oned_binary_imdilate(idxdetect,ones(1,5));    
    x = trx(fly).x_mm + trx(fly).dcenter.*cos(trx(fly).theta+trx(fly).anglefrom1to2);
    y = trx(fly).y_mm + trx(fly).dcenter.*sin(trx(fly).theta+trx(fly).anglefrom1to2);
    
    [idxstarts,idxends] = get_interval_ends(idxplotother);
    idxplotconnector = false(1,trx(fly).nframes);
    for i = 1:length(idxstarts),
      idxplotconnector([idxstarts(i):5:idxends(i)-2,idxends(i)-1]) = true;
    end
    hconnect = plot([x(idxplotconnector);trx(fly).x_mm(idxplotconnector)],...
      [y(idxplotconnector);trx(fly).y_mm(idxplotconnector)],'--','color',[1,1,0]*.7);
    
    % plot other fly's position
    hother = plot(x(idxplotother),y(idxplotother),'g.-');

  end
  
  % plot the fly's position
  hmain = plot(trx(fly).x_mm,trx(fly).y_mm,'k.-');
    
  axisalmosttight;
  axis equal;

  set(hax(plotcurr),'xtick',[],'ytick',[]);
  
  % increment plot
  plotcurr = plotcurr + 1;
  if plotcurr > nperpagecurr && page ~= npages,
    plotcurr = 1;
    page = page + 1;
    if npages == page,
      nperpagecurr = nlastpage;
    else
      nperpagecurr = nperpage;
    end
  end
  
end

if plotcurr > 1,
  delete(hax(plotcurr:end));
end

fprintf('Illustration of detected behaviors\n');
fprintf('The fly''s position is plotted in black\n');
fprintf('We highlight frames in which the behavior is detected in magenta.\n');
fprintf('At the start and end of each detected behavior, we plot\n');
fprintf('a green and red circle, resp. There is a subplot for labeled fly.\n');
if issocial,
  fprintf('During and near detected behaviors, we plot the position of the other fly\n');
  fprintf('in green. We connect corresponding frames for the main and other fly in yellow.\n');
end

succeeded = true;
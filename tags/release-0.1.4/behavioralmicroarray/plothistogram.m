function plotstuff = plothistogram(histstuff,varargin)
%fighandles has member variables {'fig','haxes'}
%plotmode in {'Mean Fraction per Fly','Mean Count per Fly','Total Count'}
%plotindivs is true/false
%plotstd in {'off','standard deviation','standard error'}
%flycolors

[plotmode,plotstuff,plotindivs,plotstd,logplot] = ...
  myparse(varargin,'plotmode','Mean Fraction per Fly','fighandles',struct,...
  'plotindivs',false,'plotstd','standard deviation','logplot',false);

% total count has no single-fly stuff
if strcmpi(plotmode,'total count') || strcmpi(plotmode,'total fraction'),
  plotindivs = false;
  plotstd = 'off';
end

nprops = length(histstuff.centers);
if nprops == 1,
  nflies = size(histstuff.countsperfly,1);
  ndata = size(histstuff.countsperfly,3);
else
  nflies = size(histstuff.countsperfly,3);
  ndata = size(histstuff.countsperfly,4);
end
if ~isfield(plotstuff,'fig') || ~ishandle(plotstuff.fig),
  plotstuff.fig = figure;
  tmp = get(plotstuff.fig,'position');
  desiredwidth = 800;
  if tmp(3) < desiredwidth,
    tmp(1) = tmp(1) - (desiredwidth-tmp(3))/2;
    tmp(3) = desiredwidth;
  end
  if tmp(4) < desiredwidth,
    tmp(2) = tmp(2) - (desiredwidth-tmp(4))/2;
    tmp(4) = desiredwidth;
  end
  set(plotstuff.fig,'position',tmp);
else
  figure(plotstuff.fig);
  clf;
end

switch lower(plotmode),
  case 'mean count per fly',
    perflydata = histstuff.countsperfly;
    meandata = histstuff.countsmean;
    stddata = histstuff.countsstd;
    stderrdata = histstuff.countsstderr;
  case 'mean fraction per fly',
    perflydata = histstuff.fracperfly;
    meandata = histstuff.fracmean;
    stddata = histstuff.fracstd;
    stderrdata = histstuff.fracstderr;
  case 'total count',
    perflydata = nan;
    meandata = histstuff.countstotal;
    stddata = nan;
    stderrdata = nan;
  case 'total fraction',
    perflydata = nan;
    meandata = histstuff.fractotal;
    stddata = nan;
    stderrdata = nan;
end

% plus and minus one standard deviation or standard error
if strcmpi(plotstd,'standard error')
  above = meandata + stderrdata;
  below = meandata - stderrdata;
elseif strcmpi(plotstd,'standard deviation')
  above = meandata + stddata;
  below = meandata - stddata;
end
if logplot,
  epsilon = min(meandata(meandata > 0));
  if ~strcmpi(plotstd,'off'),
    epsilon = min(epsilon,min(min([below(below>0);above(above>0)])));
  end
  if plotindivs,
    epsilon = min(epsilon,min(perflydata(perflydata>0)));
  end
  if isempty(epsilon),
    meandata(:) = 1;
    if ~strcmpi(plotstd,'off'),
      below(:) = 1;
      above(:) = 1;
    end
    if plotindivs,
      perflydata(:) = 1;
    end
  else
    meandata = max(meandata,epsilon);
    if ~strcmpi(plotstd,'off'),
      below = max(below,epsilon);
      above = max(above,epsilon);
    end
    if plotindivs,
      perflydata = max(perflydata,epsilon);
    end
  end
end

% plot one property histogram
if nprops == 1,
  
  maincolors = zeros(ndata,3);
  if ndata > 1,
    maincolors(2:end,:) = lines(ndata-1);
    % take gray and make it black
    maincolors(all(maincolors==.25,2),:) = 0;
  end
  sigcolors = (2+maincolors)/3;

  plotstuff.haxes = axes;
  hold on;
  %plotstuff.haxes = zeros(ndata,1);
  %for datai = 1:ndata,
  %  plotstuff.haxes(datai) = subplot(ndata,1,datai);
  %  hold on;
  %end

  % plot standard deviation/standard error
  if ~strcmpi(plotstd,'off'),
    plotstuff.hsig = nan(1,ndata);
    px = [histstuff.centers{1},fliplr(histstuff.centers{1})];
    axes(plotstuff.haxes);
    for datai = 1:ndata,
      py = [above(:,:,datai),fliplr(below(:,:,datai))];
      plotstuff.hsig(datai) = patch(px,py,sigcolors(datai,:),'linestyle','none');
    end
  end
  
  % plot individual flies
  if plotindivs,
    
    % choose colors for flies
    % find projection on principal component
    %tmp = reshape(perflydata,[size(perflydata,1),size(perflydata,2)*size(perflydata,3)]);
    %coeff = princomp(tmp');
    %coeff = coeff(:,1);
    % sort flies by this coefficient
    %[tmp,order] = sort(coeff);
    %[tmp,order] = sort(order);
    if ndata == 1,
      % use different colors for each fly if there is only one curve
      plotstuff.flycolors = (3+jet(nflies))/4*.7;
    else
      plotstuff.flycolors = zeros([nflies,3,ndata]);
      for datai = 1:ndata,
        nfliescurr = nnz(histstuff.istype(datai,:));
        if all(maincolors(datai,:)==0),
          % if black, then make gray
          plotstuff.flycolors(histstuff.istype(datai,:),:,datai) = ...
              repmat(linspace(.25,.75,nfliescurr)',[1,3]);
        else
          % mess with the saturation
          tmp = rgb2hsv(maincolors(datai,:));
          news = linspace(.25,.75,nfliescurr)';
          plotstuff.flycolors(histstuff.istype(datai,:),:,datai) = ...
              hsv2rgb([repmat(tmp(:,1),[nfliescurr,1]),news,repmat(tmp(:,3),[nfliescurr,1])]);
        end
      end
    end
    % sort
    %plotstuff.flycolors = plotstuff.flycolors(order,:,:);
    plotstuff.hindiv = zeros(ndata,nflies);
    for datai = 1:ndata,
      %axes(plotstuff.haxes(datai));
      for fly = 1:nflies,
        if ~histstuff.istype(datai,fly),
          continue;
        end
        plotstuff.hindiv(datai,fly) = ...
            plot(histstuff.centers{1},perflydata(fly,:,datai),'-','color',...
                 plotstuff.flycolors(fly,:,datai));
      end
    end
  end
  
  % plot mean
  plotstuff.hmu = zeros(1,ndata);
  for datai = 1:ndata,
    %axes(plotstuff.haxes(datai));
    plotstuff.hmu(datai) = plot(histstuff.centers{1},meandata(:,:,datai),'o-',...
                                'linewidth',2,'markerfacecolor',maincolors(datai,:),'markersize',6,...
                                'color',maincolors(datai,:));
  end
  axisalmosttight;
  
  if logplot,
    ylim = get(plotstuff.haxes,'ylim');
    ylim(1) = max(ylim(1),epsilon);
    set(plotstuff.haxes,'ylim',ylim);
    set(plotstuff.haxes,'yscale','log');
  else
    set(plotstuff.haxes,'yscale','linear');
  end
  
  xlabel(strrep(histstuff.propnames{1},'_','\_'));
  ylabel(strrep(plotmode,'_','\_'));
  legend(plotstuff.hmu,histstuff.dataname{:});
  
else
  % two properties
  
  % create axes

  % if we're plotting individuals, compute a good number of rows, columns
  plotstuff.haxes = struct;
  if plotindivs,
    nindivrow = zeros(1,ndata);
    nindivcol = zeros(1,ndata);
    pos = get(plotstuff.fig,'position');
    w = pos(3); h = pos(4);
    for datai = 1:ndata,
      nfliescurr = nnz(histstuff.istype(datai,:));
      nindivcol(datai) = max(1,round(sqrt(2*nfliescurr*w/h)));
      nindivrow(datai) = ceil(nfliescurr/nindivcol(datai));
    end
  end
  if strcmpi(plotstd,'off') && ~plotindivs,
    % mean only
    plotstuff.haxes.mean = zeros(ndata,1);
    for datai = 1:ndata,
      plotstuff.haxes.mean(datai) = subplot(ndata,1,datai);
    end
  elseif ~strcmpi(plotstd,'off') && ~plotindivs,
    % mean, plus and minus std
    plotstuff.haxes.minus = zeros(ndata,1);
    plotstuff.haxes.mean = zeros(ndata,1);
    plotstuff.haxes.plus = zeros(ndata,1);
    for datai = 1:ndata,
      plotstuff.haxes.minus(datai) = subplot(ndata,3,sub2ind([3,ndata],1,datai));
      plotstuff.haxes.mean(datai) = subplot(ndata,3,sub2ind([3,ndata],2,datai));
      plotstuff.haxes.plus(datai) = subplot(ndata,3,sub2ind([3,ndata],3,datai));
    end
  elseif strcmpi(plotstd,'off') && plotindivs,
    % mean, indivs
    plotstuff.haxes.mean = zeros(ndata,1);
    plotstuff.haxes.indiv = zeros(ndata,nflies);
    for datai = 1:ndata,
      % first row is all the means
      plotstuff.haxes.mean(datai) = subplot(ndata+1,ndata,datai);
      % next ndata rows are for individual flies
      nfliescurr = nnz(histstuff.istype(datai,:));
      for i = 1:nfliescurr,
        plotstuff.haxes.indiv(datai,i) = subplot(nindivrow(datai)*(ndata+1),...
                                           nindivcol(datai),...
                                           nindivrow(datai)*nindivcol(datai)*datai+i);
      end
    end
  elseif ~strcmpi(plotstd,'off') && plotindivs,
    % mean, plus and minus std, indivs
    plotstuff.haxes.minus = zeros(ndata,1);
    plotstuff.haxes.mean = zeros(ndata,1);
    plotstuff.haxes.plus = zeros(ndata,1);
    plotstuff.haxes.indiv = zeros(ndata,nflies);
    for datai = 1:ndata,
      plotstuff.haxes.minus(datai) = subplot(2*ndata,3,2*3*(datai-1)+1);
      plotstuff.haxes.mean(datai) = subplot(2*ndata,3,2*3*(datai-1)+2);
      plotstuff.haxes.plus(datai) = subplot(2*ndata,3,2*3*(datai-1)+3);
      nfliescurr = nnz(histstuff.istype(datai,:));
      for i = 1:nfliescurr,
        plotstuff.haxes.indiv(datai,i) = subplot(nindivrow(datai)*2*ndata,...
                                                 nindivcol(datai),...
                                                 nindivrow(datai)*nindivcol(datai)*(2*datai-1)+i);
      end
    end
  end
  
  % plot log of counts/frac
  if logplot,
    meandata = log(meandata);
    if ~strcmpi(plotstd,'off'),
      below = log(below);
      above = log(above);
    end
    if plotindivs,
      perflydata = log(perflydata);
    end
  end
  
  % plot the mean
  x = [histstuff.centers{2}(1),histstuff.centers{2}(end)];
  y = [histstuff.centers{1}(1),histstuff.centers{1}(end)];
  maxv = max(meandata(:));
  minv = min(meandata(:));
  for datai = 1:ndata,
    axes(plotstuff.haxes.mean(datai));
    plotstuff.him.mean = imagesc(x,y,meandata(:,:,:,datai),[minv,maxv]);
    if datai == ndata,
      xlabel(strrep(histstuff.propnames{2},'_','\_'));
    end
    colormap jet;
    if logplot,
      if ndata > 1,
        title(strrep(sprintf('Log(%s), %s',plotmode,histstuff.dataname{datai}),'_','\_'));
      else
        title(strrep(sprintf('Log(%s), %s',plotmode,histstuff.dataname{datai}),'_','\_'));
      end
    else
      if ndata > 1,
        title(strrep(sprintf('%s, %s',plotmode,histstuff.dataname{datai}),'_','\_'));
      else
        title(strrep(plotmode,'_','\_'));
      end
    end
  end
  if ~strcmpi(plotstd,'off'),
    minv = min(minv,min(below(:)));
    maxv = max(maxv,max(above(:)));
    for datai = 1:ndata,
      axes(plotstuff.haxes.minus(datai));
      plotstuff.him.minus(datai) = imagesc(x,y,below(:,:,:,datai),[minv,maxv]);
      if datai == 1,
        plotstuff.hcolorbar = colorbar('eastoutside');
      end
      if datai == ndata,
        ylabel(strrep(histstuff.propnames{1},'_','\_'));
      end
      if datai == 1,
        if logplot,
          title(strrep(sprintf('Log(-1 %s)',plotstd),'_','\_'));
        else
          title(strrep(['-1 ',plotstd],'_','\_'));
        end
      end
      axes(plotstuff.haxes.plus(datai));
      plotstuff.him.plus = imagesc(x,y,above(:,:,datai),[minv,maxv]);
      set(plotstuff.haxes.mean(datai),'clim',[minv,maxv]);
      if datai == 1,
        if logplot,
          title(strrep(sprintf('Log(-1 %s)',plotstd),'_','\_'));
        else
          title(strrep(['+1 ',plotstd],'_','\_'));
        end
      end
    end
  else
    ylabel(strrep(histstuff.propnames{1},'_','\_'));
    plotstuff.hcolorbar = colorbar;
  end
  set(plotstuff.hcolorbar,'edgecolor','w');
  if plotindivs,
    for datai = 1:ndata,
      nfliescurr = nnz(histstuff.istype(datai,:));
      for i = 1:nfliescurr,
        axes(plotstuff.haxes.indiv(datai,i));
        imagesc(x,y,perflydata(:,:,i,datai),[minv,maxv]);
        if i == 1 && ndata > 1,
          title(sprintf('Fly %d, %s',i,histstuff.dataname{datai}));
        else
          title(sprintf('Fly %d',i));
        end
        axis off;
      end
    end
  end
  % link the color maps of all the axes
  allaxes = struct2cell(plotstuff.haxes);
  allaxes = [allaxes{:}];
  plotstuff.hlink = linkprop(allaxes(:)','clim');
  
  % add a button for changing colormap
  tmp = get(plotstuff.fig,'position');
  plotstuff.hcolormapbutton = uicontrol('style','pushbutton','string',...
    'Edit Colormap...','position',[tmp(3)-128     7   107    26],...
    'callback',@(hObject, eventdata, handles) colormapeditor);
  
end
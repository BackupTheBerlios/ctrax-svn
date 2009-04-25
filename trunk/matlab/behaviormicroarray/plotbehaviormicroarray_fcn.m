%inputs:
%propsperfly: nflies x nprops
%types: binary nflies x ntypes
%optional:
%propnames: 1 x nprops cell
%typenames: 1 x ntypes cell
%plotindivs: true/false
%normalize: true/false
%errorbars: {'none','std dev','std err'}
%hfig: figure handle for population means
%hfig2: figure handle for microarray

function plotstuff = plotbehaviormicroarray_fcn(propsperfly,types,varargin)

plotstuff = struct;
[nflies,nprops] = size(propsperfly);
[nflies2,ntypes] = size(types);
if nflies2 ~= nflies,
  warning('Number of flies in propsperfly = %d, in types = %d\n',nflies,nflies2);
  return;
end

[propnames,typenames,plotindivs,normalize,errorbars,hfig,hfig2] = ...
    myparse(varargin,'propnames',{},'typenames',{},'plotindivs',false,...
            'normalize',true,'errorbars','none','hfig',nan,'hfig2',nan);

if isempty(propnames),
  propnames = cell(1,nprops);
  for i = 1:nprops,
    propnames{i} = sprintf('%d',i);
  end
end

if isempty(typenames),
  typenames = cell(1,ntypes);
  for i = 1:ntypes,
    typenames{i} = sprintf('Type %d',i);
  end
end

if ~ishandle(hfig),
  hfig = figure;
  tmp = get(hfig,'position');
  desiredwidth = 1000;
  if normalize,
    desiredheight = 600;
  else
    desiredheight = 400;
  end
  if tmp(3) < desiredwidth,
    tmp(1) = tmp(1) - (desiredwidth-tmp(3))/2;
    tmp(3) = desiredwidth;
  end
  if tmp(4) < desiredheight,
    tmp(2) = tmp(2) - (desiredheight-tmp(4))/2;
    tmp(4) = desiredheight;
  end
  set(hfig,'position',tmp);
end
plotstuff.hfig = hfig;

figure(hfig);
clf;

% normalize data, if necessary
if normalize,
  mu_all = zeros(1,nprops);
  sig_all = zeros(1,nprops);
  nflies_all = zeros(1,nprops);
  propsperfly0 = propsperfly;
  for i = 1:nprops,
    ignoreidx = isnan(propsperfly(:,i));
    mu_all(i) = mean(propsperfly(~ignoreidx,i));
    sig_all(i) = std(propsperfly(~ignoreidx,i),1);
    propsperfly(:,i) = (propsperfly(:,i)-mu_all(i))/sig_all(i);
    nflies_all(i) = nnz(~ignoreidx);
  end
  plotstuff.mu_all = mu_all;
  plotstuff.sig_all = sig_all;
  plotstuff.nflies_all = nflies_all;

  % plot the means for all flies
  plotstuff.hall_axes = subplot(3,1,1);
  hold on;

  % individuals, lines
  if plotindivs,
    plotstuff.hall_indivs_lines = plot(1:nprops,propsperfly0,'-');
  end
  
  % errorbars
  if strcmpi(errorbars,'std dev'),
    plotstuff.hall_errorbars = errorbar(1:nprops,mu_all,sig_all,'ko');
  elseif strcmpi(errorbars,'std err'),
    plotstuff.stderr_all = sig_all./sqrt(nflies_all);
    plotstuff.hall_errorbars = errorbar(1:nprops,mu_all,plotstuff.stderr_all,'ko');
  end
    
  % means
  plotstuff.hall_means = plot(1:nprops,mu_all,'ko','markerfacecolor','k');

  % individuals
  if plotindivs,
    plotstuff.hall_indivs = plot(1:nprops,propsperfly0,'.');
  end
  tmp = propnames;
  [tmp{1:2:end}] = deal('');
  set(plotstuff.hall_axes,'xtick',1:nprops,'xticklabel',tmp);
  title('Unnormalized properties for all flies');  

  axisalmosttight;
  ax = axis;
  ax(1:2) = [.5,nprops+.5];
  axis(ax);

  
  % axis for normalized data
  plotstuff.hmain_axes = subplot(3,1,2:3);
  
else
  
  plotstuff.hmain_axes = gca;
  
end

hold on;

% compute the means, stddevs per type
nflies_type = zeros(ntypes,nprops);
mu = zeros(ntypes,nprops);
sig = zeros(ntypes,nprops);
for propi = 1:nprops,
  ignoreidx = isnan(propsperfly(:,propi));
  for typei = 1:ntypes,
    idx = ~ignoreidx & types(:,typei);
    nflies_type(typei,propi) = nnz(idx);
    mu(typei,propi) = mean(propsperfly(idx,propi));
    sig(typei,propi) = std(propsperfly(idx,propi),1);
  end
end

% individuals, lines
if plotindivs,
  plotstuff.hindivs_lines = zeros(nflies,ntypes);
  for typei = 1:ntypes,
    plotstuff.hindivs_lines(types(:,typei),typei) = plot(1:nprops,propsperfly(types(:,typei),:),'-');
  end
end

% no need to plot means, stddevs
skipmeans = normalize && plotindivs && ntypes == 1;

if ~skipmeans,
  
  colors = lines(nflies);
  
  % errorbars
  if strcmpi(errorbars,'std dev'),
    plotstuff.herrorbars = zeros(1,ntypes);
    for i = 1:ntypes,
      plotstuff.herrorbars(i) = errorbar(1:nprops,mu(i,:),sig(i,:),'o','color',colors(i,:));
    end
  elseif strcmpi(errorbars,'std err'),
    plotstuff.stderr = sig./sqrt(nflies_type);
    plotstuff.herrorbars = zeros(1,ntypes);
    for i = 1:ntypes,
      plotstuff.herrorbars(i) = errorbar(1:nprops,mu(i,:),plotstuff.stderr(i,:),'o',...
                                         'color',colors(i,:));
    end
  end
    
  % means
  plotstuff.hmeans = zeros(1,ntypes);
  for i = 1:ntypes,
    plotstuff.hmeans(i) = plot(1:nprops,mu(i,:),'o',...
                               'color',colors(i,:),'markerfacecolor',colors(i,:));
  end
  
  % individuals
  if plotindivs,
    plotstuff.hindivs = zeros(nflies,ntypes);
    for typei = 1:ntypes,
      c = (colors(typei,:)+2)/3;
      plotstuff.hindivs(types(:,typei),typei) = ...
          plot(1:nprops,propsperfly(types(:,typei),:),'.','color',c);
      set(plotstuff.hindivs_lines(types(:,typei),typei),'color',c);
    end
  end
else
  plotstuff.hindivs = plot(1:nprops,propsperfly,'o-');
  for i = 1:nflies,
    set(plotstuff.hindivs(i),'markerfacecolor',get(plotstuff.hindivs(i),'color'));
  end    
end

axisalmosttight;
ax = axis;
ax(1:2) = [.5,nprops+.5];
axis(ax);

if normalize,
  ylabel('Std devs');
end

if normalize,
  tmp = propnames;
  [tmp{2:2:end}] = deal('');
  set(plotstuff.hmain_axes,'xtick',1:nprops,'xticklabel',tmp);  
else
  set(plotstuff.hmain_axes,'xtick',1:nprops,'xticklabel',propnames);
end
fprintf('Properties, from left to right:\n');
fprintf('%s\n',propnames{:});
if normalize,
  linkaxes([plotstuff.hmain_axes,plotstuff.hall_axes],'x');
end

if skipmeans,
  s = cell(1,nflies);
  for i = 1:nflies,
    s{i} = sprintf('Fly %d',i);
  end
  plotstuff.hlegend = legend(plotstuff.hindivs,s);  
else
  plotstuff.hlegend = legend(plotstuff.hmeans,typenames);
end

% show the microarray

% set up the figure, if necessary
if ~ishandle(hfig2),
  hfig2 = figure;
  tmp = get(hfig2,'position');
  desiredwidth = 1000;
  desiredheight = 400;
  if tmp(3) < desiredwidth,
    tmp(1) = tmp(1) - (desiredwidth-tmp(3))/2;
    tmp(3) = desiredwidth;
  end
  if tmp(4) < desiredheight,
    tmp(2) = tmp(2) - (desiredheight-tmp(4))/2;
    tmp(4) = desiredheight;
  end
  set(hfig2,'position',tmp);
end
plotstuff.hfig2 = hfig2;
figure(hfig2);
clf;

% how many flies of each type
% note that some flies might be plotted more than once if they are included
% in more than one type
nfliespertype = zeros(1,ntypes);
for typei = 1:ntypes,
  nfliespertype(typei) = nnz(types(:,typei));
end

% set up the subplots to be of a width proportional to the number of flies
% of each type
p_type = zeros(ntypes,4);
hax_type = zeros(1,ntypes);
% first type
p_type(1,:) = [.1,.05,.8*nfliespertype(1)/sum(nfliespertype),.9];
% all other types
for typei = 2:ntypes,
  p_type(typei,:) = [.01 + p_type(typei-1,1) + p_type(typei-1,3),.05,.8*nfliespertype(typei)/sum(nfliespertype),.9];
end

him_type = zeros(1,ntypes);
for typei = 1:ntypes,
  hax_type(typei) = axes('position',p_type(typei,:));
  him_type(typei) = imagesc(propsperfly(types(:,typei),:)');
  set(hax_type(typei),'xtick',[],'ytick',[]);
  title(typenames{typei});
end
set(hax_type(1),'ytick',1:nprops,'yticklabel',propnames);
linkprop(hax_type,'clim');
linkaxes(hax_type,'y');
hcol = colorbar;
set(hcol,'units','normalized','position',[ 0.925    0.05    0.0236    0.9])
if normalize,
  set(get(hcol,'ylabel'),'string','Stds');
end
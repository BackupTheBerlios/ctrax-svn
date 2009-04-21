% params = systematic_learn_params2(data,labels,minxfns,maxxfns,minxclosefns,maxxclosefns,
%            minsumxfns,maxsumxfns,...)
% optional parameters:
% minminx, minmaxx, minminxclose, minmaxxclose, minminsumx, minmaxsumx,
% maxminx, maxmaxx, maxminxclose, maxmaxxclose, maxminsumx, maxmaxsumx,
% forceminx, forcemaxx, forceminxclose, forcemaxxclose, forceminsumx, forcemaxsumx,
% minr, maxr, isminseqlength, ismaxseqlength,
% costparams, rho, pctinit, minseqlengthorder, maxseqlengthorder, nsamples

function params = systematic_learn_params2(data,labels,minxfns,maxxfns,minxclosefns,maxxclosefns,...
  minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns,fps,varargin)

costparams.prepend = [0,.2];
costparams.append = [0,.2];
costparams.connect = [1,.2];
costparams.spurious = [1,.2];
costparams.preshorten = [0,.2];
costparams.postshorten = [0,.2];
costparams.lost = [1,.2];
costparams.split = [1,.2];
costparams.maxextend = 5;

nflies = length(data);

isevent0 = cell(1,nflies);
for fly = 1:nflies,
  isevent0{fly} = true(1,data(fly).nframes);
end

[lb.minx,lb.maxx,lb.minxclose,lb.maxxclose,lb.minsumx,lb.maxsumx,lb.minmeanx,lb.maxmeanx,...
  ub.minx,ub.maxx,ub.minxclose,ub.maxxclose,ub.minsumx,ub.maxsumx,ub.minmeanx,ub.maxmeanx,...
  force.minx,force.maxx,force.minxclose,force.maxxclose,force.minsumx,force.maxsumx,force.minmeanx,force.maxmeanx,...
  minr,maxr,costparams,rho,pct,minseqlengthorder,nsamples,...
  isminseqlength,maxseqlengthorder,ismaxseqlength,nitersnochange,isevent0,...
  initminxc,initmaxxc,initminxclosec,initmaxxclosec,initminsumxc,initmaxsumxc,initminmeanxc,initmaxmeanxc] = ...
  myparse(varargin,'minminx',{},'minmaxx',{},'minminxclose',{},'minmaxxclose',{},...
  'minminsumx',{},'minmaxsumx',{},'minminmeanx',{},'minmaxmeanx',{},...
  'maxminx',{},'maxmaxx',{},'maxminxclose',{},'maxmaxxclose',{},...
  'maxminsumx',{},'maxmaxsumx',{},'maxminmeanx',{},'maxmaxmeanx',{},...
  'forceminx',{},'forcemaxx',{},'forceminxclose',{},'forcemaxxclose',{},...
  'forceminsumx',{},'forcemaxsumx',{},'forceminmeanx',{},'forcemaxmeanx',{},...
  'minr',0,'maxr',10,'costparams',costparams,...
  'rho',.75,'pctinit',30,'minseqlengthorder',2,'nsamples',1000,...
  'isminseqlength',true,'maxseqlengthorder',1,'ismaxseqlengthorder',false,...
  'nitersnochange',5,'isevent0',isevent0,...
  'initminx',{},'initmaxx',{},'initminxclose',{},'initmaxxclose',{},...
  'initminsumx',{},'initmaxsumx',{},'initminmeanx',{},'initmaxmeanx',{});

for s = {'minx','maxx','minxclose','maxxclose','minsumx','maxsumx','minmeanx','maxmeanx'},
  eval(sprintf('init%s = struct(init%sc{:});',s{1},s{1}));
end
  
fns = fieldnames(lb);
for fn = fns',
  lb.(fn{1}) = struct(lb.(fn{1}){:});
end
fns = fieldnames(ub);
for fn = fns',
  ub.(fn{1}) = struct(ub.(fn{1}){:});
end
fns = fieldnames(force);
for fn = fns',
  force.(fn{1}) = struct(force.(fn{1}){:});
end

% collect data
datafns = union(minxfns,union(maxxfns,union(minxclosefns,union(maxxclosefns,union(minsumxfns,union(maxsumxfns,union(minmeanxfns,maxmeanxfns)))))));
for fn = datafns,
  pos.(fn{1}) = [];
end

isevent = cell(1,nflies);
for fly = 1:length(labels),
  isevent{fly} = false(1,data(fly).nframes-1);
  for i = 1:length(labels(fly).starts),
    isevent{fly}(labels(fly).starts(i):labels(fly).ends(i)-1) = true;
  end
  for fn = datafns,
    pos.(fn{1}) = [pos.(fn{1}),data(fly).(fn{1})(isevent{fly})];
  end
end

% compute minima
minpos = struct;
for fn = minxfns,
  minpos.(fn{1}) = min(pos.(fn{1}));
end

% compute maxima
maxpos = struct;
for fn = maxxfns,
  maxpos.(fn{1}) = max(pos.(fn{1}));
end

% sequence length bounds
seqlengths = getstructarrayfield(labels,'ends')-getstructarrayfield(labels,'starts');
if isminseqlength
  minseqlength = prctile(seqlengths,minseqlengthorder);
else
  minseqlength = 1;
end
if ismaxseqlength
  maxseqlength = prctile(seqlengths,100-minseqlengthorder);
else
  maxseqlength = inf;
end

% compute sums, means for positives
meanxfns = union(minmeanxfns,maxmeanxfns);
fns = union(minsumxfns,union(maxsumxfns,meanxfns));

for fn = fns,
  sumx.(fn{1}) = [];
  meanx.(fn{1}) = [];
end
ii = 1;
for fly = 1:length(labels),
  if isempty(labels(fly).starts), continue; end
  for i = 1:length(labels(fly).starts),
    t0 = labels(fly).starts(i);
    t1 = labels(fly).ends(i)-1;
    doflip = 0;
    for fn = fns,
      if strcmpi(fn{1},'absdtheta'),
        sumx.(fn{1})(ii) = modrange(data(fly).theta(t1+1)-data(fly).theta(t0),-pi,pi)*fps;
        doflip = 2*double(sumx.(fn{1})(ii) < 0)-1;
        sumx.(fn{1})(ii) = abs(sumx.(fn{1})(ii));
      elseif strcmpi(fn{1},'dtheta'),
        sumx.(fn{1})(ii) = modrange(data(fly).theta(t1+1)-data(fly).theta(t0),-pi,pi)*fps;
      elseif strcmpi(fn{1},'smoothdtheta'),
        sumx.(fn{1})(ii) = modrange(data(fly).smooththeta(t1+1)-data(fly).smooththeta(t0),-pi,pi)*fps;
      elseif strcmpi(fn{1},'abssmoothdtheta'),
        sumx.(fn{1})(ii) = abs(modrange(data(fly).smooththeta(t1+1)-data(fly).smooththeta(t0),-pi,pi))*fps;
      elseif strcmpi(fn{1},'corfrac'),
        sumx.(fn{1})(ii) = sum(data(fly).corfrac(1,t0:t1));
      elseif strcmpi(fn{1},'flipdv_cor'),
        sumx.(fn{1}) = sum(data(fly).dv_cor(t0:t1));
      elseif ~isempty(strfind(fn{1},'abs')) && ~strcmpi(fn{1},'absdv_cor'),
        fnraw = strrep(fn{1},'abs','');
        sumx.(fn{1})(ii) = abs(sum(data(fly).(fnraw)(t0:t1)));
      else
        sumx.(fn{1})(ii) = sum(data(fly).(fn{1})(t0:t1));
      end
    end
    if ismember('flipdv_cor',fns),
      if doflip == 0,
        sumx.flipdv_cor = abs(sumx.flipdv_cor);
      elseif doflip == 1,
        sumx.flipdv_cor = -sumx.flipdv_cor;
      end
    end
    for fn = fns,
      if ismember(fn{1},meanxfns),
        meanx.(fn{1})(ii) = sumx.(fn{1})(ii)/(t1-t0+1);
      end
    end
    ii = ii + 1;
  end
end


% compute minima and maxima for sums
minsumpos = struct;
for fn = minsumxfns,
  minsumpos.(fn{1}) = min(sumx.(fn{1}));
end
maxsumpos = struct;
for fn = maxsumxfns,
  maxsumpos.(fn{1}) = max(sumx.(fn{1}));
end

% compute minima and maxima for means
minmeanpos = struct;
for fn = minmeanxfns,
  minmeanpos.(fn{1}) = min(meanx.(fn{1}));
end
maxmeanpos = struct;
for fn = maxmeanxfns,
  maxmeanpos.(fn{1}) = max(meanx.(fn{1}));
end


% choose mu and sig
mu = struct;
sig = struct;

% for lower bounds
for fn = minxfns,
  if isfield(force.minx,fn{1}),
    mu.minx.(fn{1}) = force.minx.(fn{1});
    sig.minx.(fn{1}) = 0;
    continue;
  end
  lb1 = minpos.(fn{1});
  ub1 = prctile(pos.(fn{1}),pct);
  if isfield(lb.minx,fn{1}),
    lb1 = max(lb1,lb.minx.(fn{1}));
    ub1 = max(ub1,lb.minx.(fn{1}));
  end
  if isfield(ub.minx,fn{1}),
    ub1 = min(ub1,ub.minx.(fn{1}));
    lb1 = min(lb1,ub.minx.(fn{1}));
  end
  mu.minx.(fn{1}) = (lb1+ub1)/2;
  sig.minx.(fn{1}) = (ub1-lb1)/4;
end
% for upper bounds
for fn = maxxfns,
  if isfield(force.maxx,fn{1}),
    mu.maxx.(fn{1}) = force.maxx.(fn{1});
    sig.maxx.(fn{1}) = 0;
    continue;
  end
  lb1 = prctile(pos.(fn{1}),100-pct);
  ub1 = maxpos.(fn{1});
  if isfield(lb.maxx,fn{1}),
    lb1 = max(lb1,lb.maxx.(fn{1}));
    ub1 = max(ub1,lb.maxx.(fn{1}));
  end
  if isfield(ub.maxx,fn{1}),
    ub1 = min(ub1,ub.maxx.(fn{1}));
    lb1 = min(lb1,ub.maxx.(fn{1}));
  end
  mu.maxx.(fn{1}) = (lb1+ub1)/2;
  sig.maxx.(fn{1}) = (ub1-lb1)/4;
end

% for close lower bounds
for fn = minxclosefns,
  if isfield(force.minxclose,fn{1}),
    mu.minxclose.(fn{1}) = force.minxclose.(fn{1});
    sig.minxclose.(fn{1}) = 0;
    continue;
  end
  lb1 = minpos.(fn{1});
  ub1 = prctile(pos.(fn{1}),pct);
  if isfield(lb.minxclose,fn{1}),
    lb1 = max(lb1,lb.minxclose.(fn{1}));
    ub1 = max(ub1,lb.minxclose.(fn{1}));
  end
  if isfield(ub.minxclose,fn{1}),
    ub1 = min(ub1,ub.minxclose.(fn{1}));
    lb1 = min(lb1,ub.minxclose.(fn{1}));
  end
  mu.minxclose.(fn{1}) = (3*lb1+ub1)/4;
  sig.minxclose.(fn{1}) = (ub1-lb1)/8;
end
% for close upper bounds
for fn = maxxclosefns,
  if isfield(force.maxxclose,fn{1}),
    mu.maxxclose.(fn{1}) = force.maxxclose.(fn{1});
    sig.maxxclose.(fn{1}) = 0;
    continue;
  end
  lb1 = prctile(pos.(fn{1}),100-pct);
  ub1 = maxpos.(fn{1});
  if isfield(lb.maxxclose,fn{1}),
    lb1 = max(lb1,lb.maxxclose.(fn{1}));
    ub1 = max(ub1,lb.maxxclose.(fn{1}));
  end
  if isfield(ub.maxx,fn{1}),
    ub1 = min(ub1,ub.maxxclose.(fn{1}));
    lb1 = min(lb1,ub.maxxclose.(fn{1}));
  end
  mu.maxxclose.(fn{1}) = (lb1+3*ub1)/4;
  sig.maxxclose.(fn{1}) = (ub1-lb1)/4;
end

% for lower sum bounds
for fn = minsumxfns,
  if isfield(force.minsumx,fn{1}),
    mu.minsumx.(fn{1}) = force.minsumx.(fn{1});
    sig.minsumx.(fn{1}) = 0;
    continue;
  end
  mid = minsumpos.(fn{1});
  ub1 = prctile(sumx.(fn{1}),pct);
  lb1 = 2*mid - ub1;
  if isfield(lb.minsumx,fn{1}),
    lb1 = max(lb1,lb.minsumx.(fn{1}));
    ub1 = max(ub1,lb.minsumx.(fn{1}));
  end
  if isfield(ub.minsumx,fn{1}),
    ub1 = min(ub1,ub.minsumx.(fn{1}));
    lb1 = min(lb1,ub.minsumx.(fn{1}));
  end
  mu.minsumx.(fn{1}) = (lb1+ub1)/2;
  sig.minsumx.(fn{1}) = (ub1-lb1)/2;
end
% for upper sum bounds
for fn = maxsumxfns,
  if isfield(force.maxsumx,fn{1}),
    mu.maxsumx.(fn{1}) = force.maxsumx.(fn{1});
    sig.maxsumx.(fn{1}) = 0;
    continue;
  end
  mid = maxsumpos.(fn{1});
  lb1 = prctile(sumx.(fn{1}),100-pct);
  ub1 = 2*mid - lb1;
  if isfield(lb.maxsumx,fn{1}),
    lb1 = max(lb1,lb.maxsumx.(fn{1}));
    ub1 = max(ub1,lb.maxsumx.(fn{1}));
  end
  if isfield(ub.maxsumx,fn{1}),
    ub1 = min(ub1,ub.maxsumx.(fn{1}));
    lb1 = min(lb1,ub.maxsumx.(fn{1}));
  end
  mu.maxsumx.(fn{1}) = (lb1+ub1)/2;
  sig.maxsumx.(fn{1}) = (ub1-lb1)/2;
end

% for lower mean bounds
for fn = minmeanxfns,
  if isfield(force.minmeanx,fn{1}),
    mu.minmeanx.(fn{1}) = force.minmeanx.(fn{1});
    sig.minmeanx.(fn{1}) = 0;
    continue;
  end
  mid = minmeanpos.(fn{1});
  ub1 = prctile(meanx.(fn{1}),pct);
  lb1 = 2*mid - ub1;
  if isfield(lb.minmeanx,fn{1}),
    lb1 = max(lb1,lb.minmeanx.(fn{1}));
    ub1 = max(ub1,lb.minmeanx.(fn{1}));
  end
  if isfield(ub.minmeanx,fn{1}),
    ub1 = min(ub1,ub.minmeanx.(fn{1}));
    lb1 = min(lb1,ub.minmeanx.(fn{1}));
  end
  mu.minmeanx.(fn{1}) = (lb1+ub1)/2;
  sig.minmeanx.(fn{1}) = (ub1-lb1)/4;
end
% for upper mean bounds
for fn = maxmeanxfns,
  if isfield(force.maxmeanx,fn{1}),
    mu.maxmeanx.(fn{1}) = force.maxmeanx.(fn{1});
    sig.maxmeanx.(fn{1}) = 0;
    continue;
  end
  mid = maxmeanpos.(fn{1});
  lb1 = prctile(meanx.(fn{1}),100-pct);
  ub1 = 2*mid - lb1;
  if isfield(lb.maxmeanx,fn{1}),
    lb1 = max(lb1,lb.maxmeanx.(fn{1}));
    ub1 = max(ub1,lb.maxmeanx.(fn{1}));
  end
  if isfield(ub.maxmeanx,fn{1}),
    ub1 = min(ub1,ub.maxmeanx.(fn{1}));
    lb1 = min(lb1,ub.maxmeanx.(fn{1}));
  end
  mu.maxmeanx.(fn{1}) = (lb1+ub1)/2;
  sig.maxmeanx.(fn{1}) = (ub1-lb1)/4;
end

% % minimum sequence length
% if isminseqlength,
% 
%   if isfield(force,'minseqlength'),
%     mu.minseqlength.(fn{1}) = force.minseqlength;
%     sig.minseqlength.(fn{1}) = 0;
%   else
%     mid = seqlengths(minseqlengthorder);
%     ub1 = prctile(seqlengths(minseqlengthorder+1:end),pct);
%     lb1 = 2*mid - ub1;
%     if isfield(lb,'minseqlength')
%       lb1 = max(lb1,lb.minseqlength);
%       ub1 = max(ub1,lb.minseqlength);
%     end
%     if isfield(ub,'minseqlength'),
%       ub1 = min(ub1,ub.minseqlength);
%       lb1 = min(lb1,ub.minseqlength);
%     end
%     mu.minseqlength = (lb1+ub1)/2;
%     sig.minseqlength = (ub1-lb1)/4;
%   end
% 
% end
  
% if ismaxseqlength,
% 
%   if isfield(force,'maxseqlength'),
%     mu.maxseqlength.(fn{1}) = force.maxseqlength;
%     sig.maxseqlength.(fn{1}) = 0;
%   else
%     mid = seqlengths(end-maxseqlengthorder+1);
%     lb1 = prctile(seqlengths(1:end-maxseqlengthorder),100-pct);
%     ub1 = 2*mid - lb1;
%     if isfield(lb,'maxseqlength')
%       lb1 = max(lb1,lb.maxseqlength);
%       ub1 = max(ub1,lb.maxseqlength);
%     end
%     if isfield(ub.maxseqlength),
%       ub1 = min(ub1,ub.maxseqlength);
%       lb1 = min(lb1,ub.maxseqlength);
%     end
%     mu.maxseqlength = (lb1+ub1)/2;
%     sig.maxseqlength = (ub1-lb1)/4;
%   end
% 
% end

% set initial values
fnstmp = fieldnames(initminx);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.minx.(fn) = initminx.(fn)(1);
  sig.minx.(fn) = initminx.(fn)(2);
end
fnstmp = fieldnames(initmaxx);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.maxx.(fn) = initmaxx.(fn)(1);
  sig.maxx.(fn) = initmaxx.(fn)(2);
end
fnstmp = fieldnames(initminxclose);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.minxclose.(fn) = initminxclose.(fn)(1);
  sig.minxclose.(fn) = initminxclose.(fn)(2);
end
fnstmp = fieldnames(initmaxxclose);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.maxxclose.(fn) = initmaxxclose.(fn)(1);
  sig.maxxclose.(fn) = initmaxxclose.(fn)(2);
end
fnstmp = fieldnames(initminsumx);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.minsumx.(fn) = initminsumx.(fn)(1);
  sig.minsumx.(fn) = initminsumx.(fn)(2);
end
fnstmp = fieldnames(initmaxsumx);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.maxsumx.(fn) = initmaxsumx.(fn)(1);
  sig.maxsumx.(fn) = initmaxsumx.(fn)(2);
end
fnstmp = fieldnames(initminmeanx);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.minmeanx.(fn) = initminmeanx.(fn)(1);
  sig.minmeanx.(fn) = initminmeanx.(fn)(2);
end
fnstmp = fieldnames(initmaxmeanx);
for i = 1:length(fnstmp),
  fn = fnstmp{i};
  mu.maxmeanx.(fn) = initmaxmeanx.(fn)(1);
  sig.maxmeanx.(fn) = initmaxmeanx.(fn)(2);
end


r_try = [minr,maxr];
r_mu = mean(r_try);
r_sig = max(abs(r_try - r_mu))/2;

vinit = [];
for fn = minxfns,
  vinit(end+1) = mu.minx.(fn{1});
  vinit(end+1) = sig.minx.(fn{1});
end
for fn = maxxfns,
  vinit(end+1) = mu.maxx.(fn{1});
  vinit(end+1) = sig.maxx.(fn{1});
end
for fn = minxclosefns,
  vinit(end+1) = mu.minxclose.(fn{1});
  vinit(end+1) = sig.minxclose.(fn{1});
end
for fn = maxxfns,
  vinit(end+1) = mu.maxxclose.(fn{1});
  vinit(end+1) = sig.maxxclose.(fn{1});
end
for fn = minsumxfns,
  vinit(end+1) = mu.minsumx.(fn{1});
  vinit(end+1) = sig.minsumx.(fn{1});
end
for fn = maxsumxfns,
  vinit(end+1) = mu.maxsumx.(fn{1});
  vinit(end+1) = sig.maxsumx.(fn{1});
end
for fn = minmeanxfns,
  vinit(end+1) = mu.minmeanx.(fn{1});
  vinit(end+1) = sig.minmeanx.(fn{1});
end
for fn = maxmeanxfns,
  vinit(end+1) = mu.maxmeanx.(fn{1});
  vinit(end+1) = sig.maxmeanx.(fn{1});
end

vinit(end+1) = r_mu;
vinit(end+1) = r_sig;

vinit = vinit';

if any(isnan(vinit)),
  fprintf('vinit is nan\n');
  keyboard;
end

scoref = @(x) systematic_score_params2(x,data,isevent,isevent0,minseqlength,maxseqlength,costparams,...
  minxfns,maxxfns,minxclosefns,maxxclosefns,minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns,fps);
generatef = @(x,N) systematic_gen_params2(x,N,lb,ub,minr,maxr,minxfns,maxxfns,...
  minxclosefns,maxxclosefns,minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns);
maximizef = @systematic_ml_distr;

% set up plotting
SHOWFIGURE = 1234;
figure(SHOWFIGURE);
clf;
xlabel('Iteration');
ylabel('Score');
title('Scores of samples generated ... refreshed every 30s');
hax = gca;
hold on;
hbutton = uicontrol('style','pushbutton','string','Quit','units','normalized',...
  'backgroundcolor',[.7,0,0],'foregroundcolor',[1,1,1],...
  'userdata',false);
pos = get(hbutton,'position');
pos(1) = 0.8786;
set(hbutton,'position',pos);
set(hax,'userdata',0);
set(hbutton,'callback',@(hObject,eventdata)set(hObject,'UserData',true));
realclosereq = get(SHOWFIGURE,'CloseRequestFcn');
set(SHOWFIGURE,'CloseRequestFcn','');

% optimize
[v,xbest,sbest] = cross_entropy_method(scoref,generatef,maximizef,vinit,rho,'nsamples',nsamples,...
 'niters no change',nitersnochange );
fprintf('Best Score: %f\n',sbest);
i = 1;
params = struct('minx',struct,'maxx',struct','minxclose',struct,'maxxclose',struct,...
  'minsumx',struct,'maxsumx',struct,'minmeanx',struct,'maxmeanx',struct,...
  'r',0,'minseqlength',0,'maxseqlength',0);

for fn = minxfns,
  params.minx.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = maxxfns,
  params.maxx.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = minxclosefns,
  params.minxclose.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = maxxclosefns,
  params.maxxclose.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = minsumxfns,
  params.minsumx.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = maxsumxfns,
  params.maxsumx.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = minmeanxfns,
  params.minmeanx.(fn{1}) = xbest(i);
  i = i + 1;
end
for fn = maxmeanxfns,
  params.maxmeanx.(fn{1}) = xbest(i);
  i = i + 1;
end

params.r = round(xbest(end));
params.minseqlength = minseqlength;
params.maxseqlength = maxseqlength;
params.fps = fps;

delete(hbutton);
set(SHOWFIGURE,'closerequestfcn',realclosereq);
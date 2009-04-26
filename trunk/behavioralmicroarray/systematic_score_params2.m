function score = systematic_score_params2(x,data,isevent,isevent0,minseqlength,maxseqlength,costparams,...
  minxfns,maxxfns,minxclosefns,maxxclosefns,minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns,fps)

DEBUG = 1;
SHOWFIGURE = 1234;

N = size(x,2);
score = zeros(N,1);

% crop data to labeled portions
for fly = 1:length(isevent),
  nlabels = length(isevent{fly});
  if nlabels == 0, continue; end
  datacrop(fly) = GetPartOfTrack(data(fly),1,nlabels+1);    
  isevent0{fly} = isevent0{fly}(1:nlabels);
end

% score each sample
tic;
lastplotted = 0;
chil = get(SHOWFIGURE,'children');
isax = strcmpi(get(chil,'type'),'axes');
hax = chil(isax);
hbutton = chil(~isax);
iter = get(hax,'userdata');
iter = iter+1;
set(hax,'userdata',iter);
bestscore = -inf;

for n = 1:N,

  minx = struct;
  maxx = struct;
  minxclose = struct;
  maxxclose = struct;
  minsumx = struct;
  maxsumx = struct;
  minmeanx = struct;
  maxmeanx = struct;
  
  % parse current parameter x
  i = 1;
  for fn = minxfns,
    minx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxxfns,
    maxx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = minxclosefns,
    minxclose.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxxclosefns,
    maxxclose.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = minsumxfns,
    minsumx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxsumxfns,
    maxsumx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = minmeanxfns,
    minmeanx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxmeanxfns,
    maxmeanx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  r = x(end,n);

  err = 0;
  lasterr = 0;
  
  % for printing out results
  errcounts.spurious_maxextend = 0;
  errcounts.prepend = 0;
  errcounts.append = 0;
  errcounts.connect = 0;
  errcounts.spurious = 0;
  errcounts.lost_maxextend = 0;
  errcounts.preshorten = 0;
  errcounts.postshorten = 0;
  errcounts.split = 0;
  errcounts.lost = 0;
  tmpfns = fieldnames(errcounts);
  for tmpi = 1:length(tmpfns),
    errcounts.frames.(tmpfns{tmpi}) = 0;
  end
  
  for fly = 1:length(isevent),
    nlabels = length(isevent{fly});
    if nlabels == 0, continue; end

    seg = systematic_detect_event4(datacrop(fly),'event',isevent0{fly},minx,maxx,minxclose,maxxclose,minsumx,maxsumx,...
      minmeanx,maxmeanx,r,minseqlength,maxseqlength,fps);
    if DEBUG == 2,
      global currentparams;
      currentparams = {datacrop(fly),'event',isevent0{fly},minx,maxx,minxclose,maxxclose,minsumx,maxsumx,...
        minmeanx,maxmeanx,r,minseqlength,maxseqlength};
    end
    isevent_detected = false(1,nlabels);
    for i = 1:length(seg.t1),
      isevent_detected(seg.t1(i):seg.t2(i)-1) = true;
    end
    
    isboth = isevent_detected & isevent{fly};
    
    % extra detections
    isfalsepos = isevent_detected & ~isevent{fly};
    [starts,ends] = get_interval_ends(isfalsepos);
    ends = ends - 1;
    
    for i = 1:length(starts),
      w = ends(i)-starts(i)+1;
      % spurious
      if w > costparams.maxextend,
        err = err + costparams.spurious(1) + w*costparams.spurious(2);
        errcounts.spurious_maxextend = errcounts.spurious_maxextend + 1;
        errcounts.frames.spurious_maxextend = errcounts.frames.spurious_maxextend + w;
      elseif (starts(i) > 1) && isboth(starts(i)-1) && (ends(i) == nlabels || ~isevent_detected(ends(i)+1)),
        % prepend
        err = err + costparams.prepend(1) + w*costparams.prepend(2);
        errcounts.prepend = errcounts.prepend + 1;
        errcounts.frames.prepend = errcounts.frames.prepend + w;
      elseif (ends(i) < nlabels) && isboth(ends(i)+1) && (starts(i) == 1 || ~isevent_detected(starts(i)-1)),
        % append
        err = err + costparams.append(1) + w*costparams.append(2);
        errcounts.append = errcounts.append + 1;
        errcounts.frames.append = errcounts.frames.append + w;
      elseif (starts(i) > 1) && (ends(i) < nlabels) && isboth(starts(i)-1) && isboth(ends(i)+1),
        % connect
        err = err + costparams.connect(1) + w*costparams.connect(2);
        errcounts.connect = errcounts.connect + 1;
        errcounts.frames.connect = errcounts.frames.connect + w;
      else
        % spurious
        err = err + costparams.spurious(1) + w*costparams.spurious(2);
        errcounts.spurious = errcounts.spurious + 1;
        errcounts.frames.spurious = errcounts.frames.spurious + w;
      end
    end
    
    % missed detections
    isfalseneg = ~isevent_detected & isevent{fly};
    [starts,ends] = get_interval_ends(isfalseneg);
    ends = ends - 1;
   
    for i = 1:length(starts),
      w = ends(i)-starts(i)+1;
      if w > costparams.maxextend,
        % lost
        err = err + costparams.lost(1) + w*costparams.lost(2);
        errcounts.lost_maxextend = errcounts.lost_maxextend + 1;
        errcounts.frames.lost_maxextend = errcounts.frames.lost_maxextend + w;
      elseif (starts(i) > 1) && isboth(starts(i)-1) && (ends(i) == nlabels || ~isevent{fly}(ends(i)+1)),
        % preshorten
        err = err + costparams.preshorten(1) + w*costparams.preshorten(2);
        errcounts.preshorten = errcounts.preshorten + 1;
        errcounts.frames.preshorten = errcounts.frames.preshorten + w;
      elseif (ends(i) < nlabels) && isboth(ends(i)+1) && (starts(i) == 1 || ~isevent{fly}(starts(i)-1)),
        % postshorten
        err = err + costparams.postshorten(1) + w*costparams.postshorten(2);
        errcounts.postshorten = errcounts.postshorten + 1;
        errcounts.frames.postshorten = errcounts.frames.postshorten + w;
      elseif (starts(i) > 1) && (ends(i) < nlabels) && isboth(starts(i)-1) && isboth(ends(i)+1),
        % split
        err = err + costparams.split(1) + w*costparams.split(2);
        errcounts.split = errcounts.split + 1;
        errcounts.frames.split = errcounts.frames.split + w;
      else
        % lost
        err = err + costparams.lost(1) + w*costparams.lost(2);
        errcounts.lost = errcounts.lost + 1;
        errcounts.frames.lost = errcounts.frames.lost + w;
      end
    end
        
    if DEBUG == 2,
      
      sumxfns = unionmany(minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns);
      
      clf;
      fns = unionmany(minxfns,maxxfns,minxclosefns,maxxclosefns,minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns);
      nfn = length(fns);
      hax = createsubplots(nfn,1,[.03,.07]);
      ax = cell(1,nfn);
      for axi = 1:nfn,
        fn = fns{axi};
        if ~isempty(strfind(fn,'theta')),
          fact = 180/pi;
        else
          fact = 1;
        end
        axes(hax(axi));
        hold on;
        plot(fact*data(fly).(fn)(1:nlabels),'k.-');
        plot(find(isboth),fact*data(fly).(fn)(isboth),'go','markerfacecolor','g');
        plot(find(isfalsepos),fact*data(fly).(fn)(isfalsepos),'ro','markerfacecolor','r');
        plot(find(isfalseneg),fact*data(fly).(fn)(isfalseneg),'bo','markerfacecolor','b');
        if ismember(fn,minxfns),
          plot([1,nlabels],[minx.(fn),minx.(fn)]*fact,'m-');
        end
        if ismember(fn,maxxfns),
          plot([1,nlabels],[maxx.(fn),maxx.(fn)]*fact,'c-');
        end
        if ismember(fn,minxclosefns),
          plot([1,nlabels],[minxclose.(fn),minxclose.(fn)]*fact,'-','color',[1,0,1]/2);
        end
        if ismember(fn,maxxclosefns),
          plot([1,nlabels],[maxxclose.(fn),maxxclose.(fn)]*fact,'-','color',[0,1,1]/2);
        end
        
        ax{axi} = [1,nlabels,-inf,inf];
        if ismember(fn,minxclosefns),
          ax{axi}(3) = minxclose.(fn);
        elseif ismember(fn,minxfns),
          ax{axi}(3) = minx.(fn);
        end
        ax{axi} = max(ax{axi},min(data(fly).(fn)(1:nlabels)));

        if ismember(fn,maxxclosefns),
          ax{axi}(4) = maxxclose.(fn);
        elseif ismember(fn,maxxfns),
          ax{axi}(4) = maxx.(fn);
        end
        ax{axi}(4) = min(ax{axi}(4),max(data(fly).(fn)(1:nlabels)));

        tmp = ax{axi}(4)-ax{axi}(3);
        ax{axi}(3) = ax{axi}(3) - .1*tmp; 
        ax{axi}(4) = ax{axi}(4) + .1*tmp;
        
        ax{axi}(3:4) = ax{axi}(3:4)*fact;
        
        s = [fn,': '];
        if ismember(fn,minsumxfns) && ~ismember(fn,maxsumxfns),
          s = [s,sprintf('%f <= sumx, ',minsumx.(fn)*fact)];
        elseif ismember(fn,minsumxfns) && ismember(fn,maxsumxfns),
          s = [s,sprintf('%f <= sumx <= %f, ',minsumx.(fn)*fact,fact*maxsumx.(fn))];
        elseif ~ismember(fn,minsumxfns) && ismember(fn,maxsumxfns),
          s = [s,sprintf('sumx <= %f, ',maxsumx.(fn)*fact)];
        end
        if ismember(fn,minmeanxfns) && ~ismember(fn,maxmeanxfns),
          s = [s,sprintf('%f <= meanx',fact*minmeanx.(fn))];
        elseif ismember(fn,minmeanxfns) && ismember(fn,maxmeanxfns),
          s = [s,sprintf('%f <= meanx <= %f',fact*minmeanx.(fn),fact*maxmeanx.(fn))];
        elseif ~ismember(fn,minmeanxfns) && ismember(fn,maxmeanxfns),
          s = [s,sprintf('meanx <= %f',maxmeanx.(fn)*fact)];
        end
        title(s);
        
        axis(ax{axi});
      end
      for k = 1:length(seg.t1),
        t0 = seg.t1(k);
        t1 = seg.t2(k)-1;
        w = t1 - t0 + 1;
        tmid = (t0 + t1+1)/2;
        
        sumx = struct;
        doflip = 0;
        for i = 1:length(sumxfns),
          fn = sumxfns{i};
          if strcmpi(fn,'absdtheta'),
            sumx.(fn) = mod(data(fly).theta(t1+1)-data(fly).theta(t0)+pi,2*pi)-pi;
            doflip = 2*double(sumx.(fn) < 0)-1;
            sumx.(fn) = abs(sumx.(fn));
          elseif strcmpi(fn,'dtheta'),
            sumx.(fn) = mod(data(fly).theta(t1+1)-data(fly).theta(t0)+pi,2*pi)-pi;
          elseif strcmpi(fn,'smoothdtheta'),
            sumx.(fn) = mod(data(fly).smooththeta(t1+1)-data(fly).smooththeta(t0)+pi,2*pi)-pi;
          elseif strcmpi(fn,'abssmoothdtheta'),
            sumx.(fn)(i) = abs(modrange(data(fly).smooththeta(t1+1)-data(fly).smooththeta(t0),-pi,pi));
          elseif strcmpi(fn,'flipdv_cor'),
            sumx.(fn) = sum(data(fly).dv_cor(t0:t1));
          elseif ~isempty(strfind(fn,'abs')),
            fnraw = strrep(fn,'abs','');
            sumx.(fn) = abs(sum(data(fly).(fnraw)(t0:t1)));
          else
            sumx.(fn) = sum(data(fly).(fn)(t0:t1));
          end          
        end
        if ismember('flipdv_cor',sumxfns),
          if doflip == 0,
            sumx.flipdv_cor = abs(sumx.flipdv_cor);
          elseif doflip == 1,
            sumx.flipdv_cor = -sumx.flipdv_cor;
          end
        end
        for i = 1:length(sumxfns),
          fn = sumxfns{i};
          if ~isempty(strfind(fn,'theta')),
            fact = 180/pi;
          else
            fact = 1;
          end
          j = find(strcmpi(fns,fn));
          axes(hax(j));
          text(tmid,ax{j}(3)+(ax{j}(4)-ax{j}(3))*.9,{sprintf('sum = %.1f',fact*sumx.(fn)),sprintf('mu = %.1f',fact*sumx.(fn)/w)},'horizontalalignment','center','backgroundcolor','w');
        end
      end
      
      [startstmp,endstmp] = get_interval_ends(isevent{fly});
      endstmp = endstmp - 1;
      for k = 1:length(startstmp),
        t0 = startstmp(k);
        t1 = endstmp(k)-1;
        w = t1 - t0 + 1;
        tmid = (t0 + t1+1)/2;

        sumx = struct;
        doflip = 0;
        for i = 1:length(sumxfns),
          fn = sumxfns{i};
          if strcmpi(fn,'absdtheta'),
            sumx.(fn) = mod(data(fly).theta(t1+1)-data(fly).theta(t0)+pi,2*pi)-pi;
            doflip = 2*double(sumx.(fn) < 0)-1;
            sumx.(fn) = abs(sumx.(fn));
          elseif strcmpi(fn,'dtheta'),
            sumx.(fn) = mod(data(fly).theta(t1+1)-data(fly).theta(t0)+pi,2*pi)-pi;
          elseif strcmpi(fn,'smoothdtheta'),
            sumx.(fn) = mod(data(fly).smooththeta(t1+1)-data(fly).smooththeta(t0)+pi,2*pi)-pi;
          elseif strcmpi(fn,'abssmoothdtheta'),
            sumx.(fn) = abs(mod(data(fly).smooththeta(t1+1)-data(fly).smooththeta(t0)+pi,2*pi)-pi);
          elseif strcmpi(fn,'flipdv_cor'),
            sumx.(fn) = sum(data(fly).dv_cor(t0:t1));
          elseif ~isempty(strfind(fn,'abs')),
            fnraw = strrep(fn,'abs','');
            sumx.(fn) = abs(sum(data(fly).(fnraw)(t0:t1)));
          else
            sumx.(fn) = sum(data(fly).(fn)(t0:t1));
          end          
        end
        if ismember('flipdv_cor',sumxfns),
          if doflip == 0,
            sumx.flipdv_cor = abs(sumx.flipdv_cor);
          elseif doflip == 1,
            sumx.flipdv_cor = -sumx.flipdv_cor;
          end
        end
        for i = 1:length(sumxfns),
          fn = sumxfns{i};
          if ~isempty(strfind(fn,'theta')),
            fact = 180/pi;
          else
            fact = 1;
          end
          j = find(strcmpi(fns,fn));
          axes(hax(j));
          text(tmid,ax{j}(3)+(ax{j}(4)-ax{j}(3))*.1,{sprintf('sum = %f',fact*sumx.(fn)),sprintf('mu = %f',fact*sumx.(fn)/w)},'horizontalalignment','center','backgroundcolor','w');
        end
      end

      linkaxes(hax,'x');
      
 
      currerr = err - lasterr;
      lasterr = err;
      input(sprintf('Sample %d, fly = %d, err = %f',n,fly,currerr));
      
    end
    
  end

  score(n) = -err;
  if score(n) > bestscore,
    besterrcounts = errcounts;
    bestscore = score(n);
  end
  
  if toc > 30 || n == N,
    fprintf('scoring sample %d. best score so far: %f.\n',n,max(score(1:n)));
    fprintf('Best score error counts:\n');
    fprintf('Spurious sequences [extended more than %d frames]: %d, mean nframes extended: %.1f\n',...
      costparams.maxextend,besterrcounts.spurious_maxextend,...
      besterrcounts.frames.spurious_maxextend/besterrcounts.spurious_maxextend);
    fprintf('Prepended sequences: %d, mean nframes prepended: %.1f\n',...
      besterrcounts.prepend,besterrcounts.frames.prepend/besterrcounts.prepend);
    fprintf('Appended sequences: %d, mean nframes appended: %.1f\n',...
      besterrcounts.append,besterrcounts.frames.append/besterrcounts.append);
    fprintf('Connected sequences: %d, mean nframes connect: %.1f\n',...
      besterrcounts.connect,besterrcounts.frames.connect/besterrcounts.connect);
    fprintf('Spurious: %d, mean nframes spurious: %.1f\n',...
      besterrcounts.spurious,besterrcounts.frames.spurious/besterrcounts.spurious);
    fprintf('Lost sequences [shortened more than %d frames]: %d, mean nframes shortened: %.1f\n',...
      costparams.maxextend,besterrcounts.lost_maxextend,...
      besterrcounts.frames.lost_maxextend/besterrcounts.lost_maxextend);
    fprintf('Preshortened sequences: %d, mean nframes preshortened: %.1f\n',...
      besterrcounts.preshorten,besterrcounts.frames.preshorten/besterrcounts.preshorten);
    fprintf('Postshortened sequences: %d, mean nframes postshortened: %.1f\n',...
      besterrcounts.postshorten,besterrcounts.frames.postshorten/besterrcounts.postshorten);
    fprintf('Split sequences: %d, mean nframes split: %.1f\n',...
      besterrcounts.split,besterrcounts.frames.split/besterrcounts.split);
    fprintf('Lost: %d, mean nframes lost: %.1f\n',...
      besterrcounts.lost,besterrcounts.frames.lost/besterrcounts.lost);
    fprintf('\n\n');
    tmpfig = get(0,'CurrentFigure');
    set(0,'CurrentFigure',SHOWFIGURE);
    plot(repmat(iter,[1,n-lastplotted]),score(lastplotted+1:n),'k.');
    lastplotted = n;
    axisalmosttight;
    set(0,'CurrentFigure',tmpfig);
    drawnow;
    tic;
  end
  drawnow;
  if get(hbutton,'userdata'),
    % abort
    ME = MException('CrossEntropyMethod:Abort','User ended computation');
    throw(ME);
  end
  
end

if DEBUG,
  fprintf('Best Sample with score %f:\n',max(score));
  n = argmax(score);
  i = 1;
  for fn = minxfns,
    minx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxxfns,
    maxx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = minxclosefns,
    minxclose.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxxclosefns,
    maxxclose.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = minsumxfns,
    minsumx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxsumxfns,
    maxsumx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = minmeanxfns,
    minmeanx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  for fn = maxmeanxfns,
    maxmeanx.(fn{1}) = x(i,n);
    i = i + 1;
  end
  r = x(end,n);
  allfns = union(minxfns,union(maxxfns,union(minxclosefns,union(maxxclosefns,union(minsumxfns,union(maxsumxfns,union(minmeanxfns,maxmeanxfns)))))));
  fprintf('\nvariable   minxclose    minx         maxx         maxxclose    minsumx      maxsumx      minmeanx     maxmeanx\n');
  tmp = {minxclose,minx,maxx,maxxclose,minsumx,maxsumx,minmeanx,maxmeanx};
  for i = 1:length(allfns),
    fn = allfns{i};
    l = length(fn);
    fprintf([fn,repmat(' ',[1,max(0,13-l)])]);
    if ~isempty(strfind(fn,'theta')) || ~isempty(strfind(fn,'angle')) || ...
        ~isempty(strfind(fn,'phi')) || ~isempty(strfind(fn,'yaw')),
      isangle = true;
    else
      isangle = false;
    end

    for j = 1:length(tmp),
      paramcurr = tmp{j};
      if isfield(paramcurr,fn),
        if isangle,
          fprintf('%-9.2f    ',paramcurr.(fn)*180/pi);
        else
          fprintf('%-9.2f    ',paramcurr.(fn));
        end
      else
        fprintf('             ');
      end
    end
    fprintf('\n');
  end
  fprintf('r = %d, minseqlength = %d, maxseqlength = %d\n',r,minseqlength,maxseqlength);
end
% types of false positives:
% true:     0 0 0
% detected: 0 1 0    spurious
% true:     0 0 1
% detected: 0 1 0    spurious
% true:     0 0 1
% detected: 0 1 1    extend
% true:     1 0 0
% detected: 0 1 0    spurious
% true:     1 0 0
% detected: 1 1 0    extend
% true:     1 0 1
% detected: 0 1 0    spurious
% true:     1 0 1
% detected: 1 1 0    extend
% true:     1 0 1
% detected: 0 1 1    extend
% true:     1 0 1
% detected: 1 1 1    connect

% types of false negatives:
% true: 0 1 0    lost
% dete: 0 0 0
% true: 0 1 0    lost
% dete: 0 0 1
% true: 0 1 0    lost
% dete: 1 0 0
% true: 0 1 0    lost
% dete: 1 0 1
% true: 0 1 1    shortened
% dete: 0 0 1
% true: 0 1 1    shortened
% dete: 1 0 1
% true: 1 1 0    shortened
% dete: 1 0 0
% true: 1 1 0    shortened
% dete: 1 0 1
% true: 1 1 1    split
% dete: 1 0 1

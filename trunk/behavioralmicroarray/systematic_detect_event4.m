function seg = systematic_detect_event4(trk,eventname,isevent0,minx,maxx,minxclose,maxxclose,minsumx,maxsumx,minmeanx,maxmeanx,r,minseqlength,maxseqlength,fps)

% initialize
%minseqlength = max(1,minseqlength - 1);
maxnfailsinarow = 10;
nfails = 0;
minxfns = fieldnames(minx)';
maxxfns = fieldnames(maxx)';
minxclosefns = fieldnames(minxclose)';
maxxclosefns = fieldnames(maxxclose)';
minsumxfns = fieldnames(minsumx)';
maxsumxfns = fieldnames(maxsumx)';
minmeanxfns = fieldnames(minmeanx)';
maxmeanxfns = fieldnames(maxmeanx)';
meanxfns = union(minmeanxfns,maxmeanxfns);
sumxfns = union(minsumxfns,union(maxsumxfns,meanxfns));
issumflipdvcor = ismember('flipdv_cor',sumxfns);
sumxfnissmoothdtheta = strcmpi(sumxfns,'smoothdtheta');
sumxfnisabssmoothdtheta = strcmpi(sumxfns,'abssmoothdtheta');
sumxfnisabsdtheta = strcmpi(sumxfns,'absdtheta');
sumxfnisdtheta = strcmpi(sumxfns,'dtheta');
nsumxfns = length(sumxfns);
sumxfnisabs = false(1,nsumxfns);
for i = 1:nsumxfns,
  fn = sumxfns{i};
  sumxfnisabs(i) = ~isempty(strfind(fn,'abs'));
end
nclose = 0;
se = ones(1,2*r+1);
sumx = struct;
score = 0;
dobreak = false;
didfit = false;
params = {};
w = 0;

% threshold for smaller cube
isnotinsmallcube = false(1,trk.nframes-1);
for fn1 = minxfns,
  if strcmpi(fn1{1},'corfrac'),
    isnotinsmallcube = isnotinsmallcube | trk.(fn1{1})(1,1:trk.nframes-1) < minx.(fn1{1});
  else
    isnotinsmallcube = isnotinsmallcube | trk.(fn1{1})(1:trk.nframes-1) < minx.(fn1{1});
  end
end
for fn1 = maxxfns,
  if strcmpi(fn1{1},'corfrac'),
    isnotinsmallcube = isnotinsmallcube | trk.(fn1{1})(1,1:trk.nframes-1) > maxx.(fn1{1});
  else
    isnotinsmallcube = isnotinsmallcube | trk.(fn1{1})(1:trk.nframes-1) > maxx.(fn1{1});
  end
end

% threshold for larger cube
isnotinlargecube = false(1,trk.nframes-1);
for fn1 = minxclosefns,
  if strcmpi(fn1{1},'corfrac'),
    isnotinlargecube = isnotinlargecube | trk.(fn1{1})(1,1:trk.nframes-1) < minxclose.(fn1{1});
  else
    isnotinlargecube = isnotinlargecube | trk.(fn1{1})(1:trk.nframes-1) < minxclose.(fn1{1});
  end
end
for fn1 = maxxclosefns,
  if strcmpi(fn1{1},'corfrac'),
    isnotinlargecube = isnotinlargecube | trk.(fn1{1})(1,1:trk.nframes-1) > maxxclose.(fn1{1});
  else
    isnotinlargecube = isnotinlargecube | trk.(fn1{1})(1:trk.nframes-1) > maxxclose.(fn1{1});
  end
end

isinsmallcube = ~isnotinsmallcube;
dilateisinsmallcube = oned_binary_imdilate(~isnotinsmallcube,se);

% compute cumulative sums
cumsumx = struct;
for fn1 = sumxfns,
  fn = fn1{1};
  if strcmpi(fn,'corfrac'),
    cumsumx.(fn) = cumsum(trk.(fn)(1,:));
  elseif ~isempty(strfind(fn,'dtheta')),
    continue;
  elseif strcmpi(fn,'flipdv_cor'),
    cumsumx.(fn) = cumsum(trk.dv_cor);
  elseif ~isempty(strfind(fn,'abs')),
    fnraw = strrep(fn,'abs','');
    cumsumx.(fn) = cumsum(trk.(fnraw));
  else      
    cumsumx.(fn) = cumsum(trk.(fn));
  end
end

scoref = @(trk,t0,t1) scorefun(trk,t0,t1,false);
checkf = @(trk,t0,t1) scorefun(trk,t0,t1,true);

% initialize dynamic programming table
dpT = trk.nframes;
dpscore = zeros(1,dpT+1);
dpdidfitonce = false(1,dpT);
dpscore(1) = 0;
% prev(t+1) is the first frame of the last segment in the best segmentation
% of 1:t
dpprev = 0:dpT;

for dpt2 = 1:dpT,
  
  dpprev(dpt2+1) = dpt2-1;
  for dpt1 = dpt2:-1:1,
  
    % compute the score for the segment from t1 to t2
    [dpscorecurr,dpdobreak] = scoref(trk,dpt1,dpt2);
    dpdidfitonce(dpt2) = dpdidfitonce(dpt2) || dpscorecurr > 0;
    
    % if dobreak == true, then this and all future iterations will have
    % scorecurr == 0
    if dpdobreak,
      % if we never fit, then set prev to be t2-1
      if ~dpdidfitonce(dpt2),
        dpprev(dpt2+1) = dpt2-1;
        dpscore(dpt2+1) = dpscore(dpt2);
      end
      break;
    end
    
    % compute the total score for segmenting up to t2 and segmenting at t1
    % score(t1) is the score for the segment ending at t1 - 1, which is
    % stored in t1, so this is really (t1-1) + 1
    dpnewscore = dpscore(dpt1) + dpscorecurr;

    % store if this is the best so far
    if dpnewscore >= dpscore(dpt2+1),
      dpscore(dpt2+1) = dpnewscore;
      dpprev(dpt2+1) = dpt1-1;
    end
    
  end
  
end

% backtrack to get the segmentation
seg = struct('t1',{[]},'t2',{[]},'score',{[]},'params',{{}});

dpj = 1;
dpt2 = dpT;
while 1,
  if dpt2 == 0,
    break;
  end
  dpt1 = dpprev(dpt2+1)+1;
  [dpscorecurr,dpdobreak,dpdidfit,dpparamscurr] = checkf(trk,dpt1,dpt2);
  if dpdidfit,
    seg.t1(dpj) = dpt1;
    seg.t2(dpj) = dpt2;
    seg.score(dpj) = dpscorecurr;
    seg.params{dpj} = dpparamscurr;
    dpj = dpj+1;
  end
  if dpt1-1 >= dpt2,
    keyboard;
  end
  dpt2 = dpt1-1;
end
seg.t1 = fliplr(seg.t1);
seg.t2 = fliplr(seg.t2);
seg.score = fliplr(seg.score);
seg.params = fliplr(seg.params);

function [score,dobreak,didfit,params] = scorefun(trk,t0,t1,docheck)

score = 0;
dobreak = false;
didfit = false;
params = {};


% different scheme for indexing
t1 = t1 - 1;

w = t1 - t0 + 1;

if w < minseqlength,
  return;
end

if w > maxseqlength,
  dobreak = true;
  return;
end

% all frames have to be in the larger cube
if w-1 < minseqlength || w == 1 || docheck,
  nfails = 0;
  if any(isnotinlargecube(t0:t1) | ~isevent0(t0:t1)),
    dobreak = true;
    return;
  end
else
  if isnotinlargecube(t0) || ~isevent0(t0),
    dobreak = true;
    return;
  end
end

% check whether all frames are within smaller cube or can be covered by
% dilating

if w == minseqlength || w == 1 || docheck,

  if w <= r,
    
    % nclose is the number of frames outside the small cube at the
    % beginning of the sequence
    nclose = find(isinsmallcube(t0:t1),1);
    if isempty(nclose), 
      nclose = w;
      nfails = nfails + 1;
      if nfails > maxnfailsinarow,
        dobreak = true;
      end
      return;
    else
      nclose = nclose - 1;
    end
    
  else
    
    % make sure gaps in small cube are not too big
    
    % find the first frame in the small cube
    nclose = find(isinsmallcube(t0:t1),1);

    % if no frames in the small cube, then there is a gap of size > r
    if isempty(nclose),
      dobreak = true;
      return;
    end
    nclose = nclose - 1;

    % make sure that all frames after first frame in small cube can be
    % covered by dilating
    if ~all(dilateisinsmallcube(t0+nclose:t1)),
      dobreak = true;
      return;
    end

    % make sure first frame is not more than r
    if nclose > r,
      dobreak = 0;
      return;
    end
  end
    
else

  % check to see if we've made a hole larger than r in this frame
  if isnotinsmallcube(t0),
    nclose = nclose + 1;
    if nclose > r,
      dobreak = true;
      return;
    end
    if nclose == w,
      return;
    end
  else
    nclose = 0;
  end
  
end

% compute the sum
doflip = 0;
for i = 1:nsumxfns,
  fn = sumxfns{i}; 
  if sumxfnisabsdtheta(i),
    sumx.(fn) = (mod(trk.theta(t1)-trk.theta(t0)+pi,2*pi)-pi)*fps;
    doflip = 2*double(sumx.(fn) < 0)-1;
    sumx.(fn) = abs(sumx.(fn));
  elseif sumxfnisdtheta(i),
    sumx.(fn) = (mod(trk.theta(t1)-trk.theta(t0)+pi,2*pi)-pi)*fps;
  elseif sumxfnissmoothdtheta(i),
    sumx.(fn) = (mod(trk.smooththeta(t1)-trk.smooththeta(t0)+pi,2*pi)-pi)*fps;
  elseif sumxfnisabssmoothdtheta(i),
    sumx.(fn) = abs(mod(trk.smooththeta(t1)-trk.smooththeta(t0)+pi,2*pi)-pi)*fps;
  else
    if t0 > 1,
      sumx.(fn) = cumsumx.(fn)(t1) - cumsumx.(fn)(t0-1);
    else
      sumx.(fn) = cumsumx.(fn)(t1);
    end
    if sumxfnisabs(i),
      sumx.(fn) = abs(sumx.(fn));
    end
  end
end
if issumflipdvcor,
  if doflip == 0,
    sumx.flipdv_cor = abs(sumx.flipdv_cor);
  elseif doflip == 1,
    sumx.flipdv_cor = -sumx.flipdv_cor;
  end
end

% threshold the sum
for fn1 = minsumxfns,
  fn = fn1{1};
  if sumx.(fn) < minsumx.(fn),
    nfails = nfails + 1;
    if nfails > maxnfailsinarow,
      dobreak = true;
    end
    return;
  end
end
for fn1 = maxsumxfns,
  fn = fn1{1};
  if sumx.(fn) > maxsumx.(fn),
    nfails = nfails + 1;
    if nfails > maxnfailsinarow,
      dobreak = true;
    end
    return;
  end
end

% compute the mean
for fn1 = meanxfns,
  fn = fn1{1};
  meanx.(fn) = sumx.(fn) / w;
end

% threshold the mean
for fn1 = minmeanxfns,
  fn = fn1{1};
  if meanx.(fn) < minmeanx.(fn),
    nfails = nfails + 1;
    if nfails > maxnfailsinarow,
      dobreak = true;
    end
    return;
  end
end
for fn1 = maxmeanxfns,
  fn = fn1{1};
  if meanx.(fn) > maxmeanx.(fn),
    nfails = nfails + 1;
    if nfails > maxnfailsinarow,
      dobreak = true;
    end
    return;
  end
end

score = w.^2;
didfit = true;
if nargout > 3,
  params = {};%eventname,mean(trk.du_cor(t0:t1)),...
    %mean(trk.dv_cor(t0:t1)),mean(trk.dtheta(t0:t1)),...
    %mean(trk.corfrac(1,t0:t1))};
end

end

end
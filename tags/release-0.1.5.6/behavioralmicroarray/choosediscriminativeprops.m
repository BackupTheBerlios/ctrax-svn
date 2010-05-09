function chosenprops = choosediscriminativeprops(X,y,nchoose,varargin)

chosenprops = [];
MINSIG = .00001;

[nflies,ntypes] = size(y);
[nflies2,nprops] = size(X);
if nflies ~= nflies2,
  warning('Number of flies in y does not match X');
  return;
end

if nflies <= 1,
  warning('Must be more than one fly');
  return;
end

isindep = myparse(varargin,'isindependent',false);

nchoose = min(nchoose,nprops);
if nchoose == nprops,
  chosenprops = 1:nprops;
  return;
end

badidx = any(isnan(X),2);
if any(badidx),
  warning('Some flies have nan properties, excluding');
  X(badidx,:) = [];
  y(badidx,:) = [];
end

if ntypes == 1,

  allowedsplits1 = ceil(nflies*.3);
  allowedsplits2 = floor(nflies*.7);
  
  if isindep,

    cost = inf(1,nprops);
    n1 = (allowedsplits1:allowedsplits2)';
    n2 = nflies - n1;
    for propi = 1:nprops,
      % z-score
      sig = std(X(:,propi),1);
      if sig < MINSIG,
        continue;
      end
      mu = mean(X(:,propi));
      Xz = sort((X(:,propi) - mu)/sig);
      
      % prepare to compute the mean and standard deviation quickly using
      % the cumulative sum
      Xz2 = Xz.^2;
      csXz = cumsum(Xz(allowedsplits1:allowedsplits2))+sum(Xz(1:allowedsplits1-1));
      csXz2 = cumsum(Xz2(allowedsplits1:allowedsplits2))+sum(Xz2(1:allowedsplits1-1));
      sXz = sum(Xz);
      sXz2 = sum(Xz2);

      % compute the mean and variance for each allowed split
      mu1 = csXz./n1;
      mu2 = (sXz-csXz)./n2;
      sig1 = csXz2./n1 - mu1.^2;
      sig2 = (sXz2-csXz2)./n2 - mu2.^2;
      
      % compute the between class scatter
      bc = (mu2-mu1).^2;
      % compute the within class scatter
      wc = sig1.*n1 + sig2.*n2;
      
      % choose the best
      cost(propi) = min(wc ./ bc);

    end
    
    % sort by splitting cost
    [tmp,order] = sort(cost);
    chosenprops = order(1:nchoose);
    
  else
    
    chosenprops = zeros(1,nchoose);

    Xcurr = X;
    n1 = (allowedsplits1:allowedsplits2)';
    n2 = nflies - n1;
    ischosen = false(1,nprops);
    hwait = waitbar(0,'Choosing properties...');
    for choosei = 1:nchoose,
          
      waitbar((choosei-1)/nchoose,hwait,sprintf('Choosing property %d',choosei))
      bestcost = inf;
      for propi = 1:nprops,
        
        if ischosen(propi), continue; end
        
        % z-score
        sig = std(Xcurr(:,propi),1);
        if sig < MINSIG,
          continue;
        end
        mu = mean(Xcurr(:,propi));
        Xz = sort((Xcurr(:,propi) - mu)/sig);
        
        % prepare to compute the mean and standard deviation quickly using
        % the cumulative sum
        Xz2 = Xz.^2;
        csXz = cumsum(Xz(allowedsplits1:allowedsplits2))+sum(Xz(1:allowedsplits1-1));
        csXz2 = cumsum(Xz2(allowedsplits1:allowedsplits2))+sum(Xz2(1:allowedsplits1-1));
        sXz = sum(Xz);
        sXz2 = sum(Xz2);
        
        % compute the mean and variance for each allowed split
        mu1 = csXz./n1;
        mu2 = (sXz-csXz)./n2;
        sig1 = csXz2./n1 - mu1.^2;
        sig2 = (sXz2-csXz2)./n2 - mu2.^2;
      
        % compute the between class scatter
        bc = (mu2-mu1).^2;
        % compute the within class scatter
        wc = sig1.*n1 + sig2.*n2;
      
        % choose the best
        cost = min(wc ./ bc);
        
        % store
        if cost < bestcost,
          bestcost = cost;
          bestpropi = propi;
        end
        
      end

      % store
      chosenprops(choosei) = bestpropi;
      ischosen(bestpropi) = true;
      
      % remove correlation
      xchosen = [X(:,ischosen),ones(nflies,1)];
      Xcurr = X;
      for propi = 1:nprops,
        if ischosen(propi),
          Xcurr(:,propi) = 0;
        else
          [b,bint,r] = regress(X(:,propi),xchosen);
          Xcurr(:,propi) = r;
        end
      end
      
    end
    
    if ishandle(hwait),
      delete(hwait);
    end
    
  end
  
  
else
  
  Z = sum(y,1);
  
  if isindep,

    cost = inf(1,nprops);
    
    for propi = 1:nprops,
      % z-score
      sig = std(X(:,propi),1);
      if sig < MINSIG,
        continue;
      end
      mu = mean(X(:,propi));
      Xz = (X(:,propi) - mu)/sig;
      X2 = Xz.^2;

      mu = zeros(1,ntypes);
      sig = zeros(1,ntypes);
      for typei = 1:ntypes,
        mu(typei) = sum(Xz.*y(:,typei))/Z(typei);
        sig(typei) = sum(X2.*y(:,typei))/Z(typei) - mu(typei).^2;
      end
      
      % compute the between class scatter
      bc = sum(pdist(mu','euclidean').^2);
      % compute the within class scatter
      wc = sum(sig.*Z);
      
      % store
      cost(propi) = wc / bc;
      
    end

    % sort by splitting cost
    [tmp,order] = sort(cost);
    chosenprops = order(1:nchoose);
    
  else

    hwait = waitbar(0,'Choosing property 1');
    chosenprops = zeros(1,nchoose);

    Xcurr = X;
    ischosen = false(1,nprops);
    for choosei = 1:nchoose,
          
      waitbar((choosei-1)/nchoose,hwait,sprintf('Choosing property %d',choosei))
      bestcost = inf;
      for propi = 1:nprops,
        
        
        if ischosen(propi), continue; end
        
        % z-score
        sig = std(Xcurr(:,propi),1);
        if sig < MINSIG,
          continue;
        end
        mu = mean(Xcurr(:,propi));
        Xz = (Xcurr(:,propi) - mu)/sig;
        X2 = Xz.^2;

        mu = zeros(1,ntypes);
        sig = zeros(1,ntypes);
        for typei = 1:ntypes,
          mu(typei) = sum(Xz.*y(:,typei))/Z(typei);
          sig(typei) = sum(X2.*y(:,typei))/Z(typei) - mu(typei).^2;
        end
        
        % compute the between class scatter
        bc = sum(pdist(mu','euclidean').^2);
        % compute the within class scatter
        wc = sum(sig.*Z);
        % compute Fisher's linear criterion
        cost = wc / bc;
     
        % store
        if cost < bestcost,
          bestcost = cost;
          bestpropi = propi;
        end
        
      end

      % store
      chosenprops(choosei) = bestpropi;
      ischosen(bestpropi) = true;
      
      % remove correlation
      xchosen = [X(:,ischosen),ones(nflies,1)];
      Xcurr = X;
      for propi = 1:nprops,
        if ischosen(propi),
          Xcurr(:,propi) = 0;
        else
          [b,bint,r] = regress(X(:,propi),xchosen);
          Xcurr(:,propi) = r;
        end
      end
      
    end
    
    if ishandle(hwait),
      delete(hwait);
    end
    
  end
        
end

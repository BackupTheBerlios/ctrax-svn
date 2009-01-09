function chosenprops = choosediscriminativeprops(X,y,nchoose,varargin)

chosenprops = [];

[nflies,ntypes] = size(y);
[nflies2,nprops] = size(X);
if nflies ~= nflies2,
  warning('Number of flies in y does not match X');
  return;
end

isindep = myparse(varargin,'isindependent',false);

nchoose = min(nchoose,nprops);

badidx = any(isnan(X),2);
if any(badidx),
  warning('Some flies have nan properties, excluding');
  X(badidx,:) = [];
  y(badidx,:) = [];
end

if isindep,

  if ntypes == 1,

    % z-score
    for i = 1:nprops,      
      sig = std(X(:,i),1);
      mu = mean(X(:,i));
      X(:,i) = (X(:,i) - mu)/sig;
    end
    
    cost = zeros(1,nprops);
    for propi = 1:nprops,
      [mu,cost(propi)] = onedimkmeans(X(:,propi),2);
    end
    [tmp,order] = sort(cost);
    chosenprops = order(1:nchoose);
    
  else
    
    cost = zeros(1,nprops);
    for propi = 1:nprops,
      [thresh,cost(propi)] = onedim2kclass(X(:,propi),y);
    end
    [tmp,order] = sort(cost);
    chosenprops = order(1:nchoose);
    
  end
  
else
  
  if ntypes == 1,
    
    % z-score
    for i = 1:nprops,      
      keep = ~isnan(X(:,i));
      sig = std(X(keep,i),1);
      mu = mean(X(keep,i));
      X(:,i) = (X(:,i) - mu)/sig;
    end
    
    % cluster using GMM into two groups
    fprintf('Clustering...\n');
    [mu,S,priors,post] = mygmm(X,2,'replicates',100,'covartype','diag');
    [tmp,yy] = max(post,[],2);
    
    fprintf('Selecting features...\n');
    % criterion for evaluating a subset of features is the number 
    % misclassification errors of LDA
    opts = statset('display','iter');
    fun = @(XT,yT)...
      (sum(~strcmp(yT,classify(XT,XT,yT,'linear'))));
    % input X should be n x d, y should be n x 1
    [fs,history] = sequentialfs(fun,X,yy,'cv','none',...
                                'nfeatures',nchoose,...
                                'direction','forward',...
                                'options',opts);
    
    % order features from first-chosen to last
    chosenprops = [];
    for i = 1:size(history.In,1),
      chosenprops = [chosenprops,find(history.In(i,:))];
    end
    
  else

    % criterion for evaluating a subset of features is the number 
    % misclassification errors of LDA
    opts = statset('display','iter');
    fun = @(XT,yT,Xt,yt)...
      (sum(~strcmp(yt,classify(Xt,XT,yT,'linear'))));
    % input X should be n x d, y should be n x 1
    [fs,history] = sequentialfs(fun,X,y,'cv','resubstition',...
                                'nfeatures',nchoose,...
                                'direction','forward',...
                                'options',opts);
    
    % order features from first-chosen to last
    chosenprops = [];
    for i = 1:size(history.In,1),
      chosenprops = [chosenprops,find(history.In(i,:))];
    end
    
  end
  
end
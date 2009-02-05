function x = systematic_gen_params2(v,N,lb,ub,minr,maxr,minxfns,maxxfns,...
  minxclosefns,maxxclosefns,minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns)

DEBUG = true;

% parse current distribution parameters v
i = 1;
minx_mu = struct;
minx_sig = struct;
maxx_mu = struct;
maxx_sig = struct;
minxclose_mu = struct;
minxclose_sig = struct;
maxxclose_mu = struct;
maxxclose_sig = struct;
minsumx_mu = struct;
minsumx_sig = struct;
maxsumx_mu = struct;
maxsumx_sig = struct;
minmeanx_mu = struct;
minmeanx_sig = struct;
maxmeanx_mu = struct;
maxmeanx_sig = struct;
for fn = minxfns,
  minx_mu.(fn{1}) = v(i);
  minx_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = maxxfns,
  maxx_mu.(fn{1}) = v(i);
  maxx_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = minxclosefns,
  minxclose_mu.(fn{1}) = v(i);
  minxclose_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = maxxclosefns,
  maxxclose_mu.(fn{1}) = v(i);
  maxxclose_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = minsumxfns,
  minsumx_mu.(fn{1}) = v(i);
  minsumx_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = maxsumxfns,
  maxsumx_mu.(fn{1}) = v(i);
  maxsumx_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = minmeanxfns,
  minmeanx_mu.(fn{1}) = v(i);
  minmeanx_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
for fn = maxmeanxfns,
  maxmeanx_mu.(fn{1}) = v(i);
  maxmeanx_sig.(fn{1}) = v(i+1);
  i = i + 2;
end
r_mu = v(end-1);
r_sig = v(end);

% generate lower bound
[mu,sig,lb1,ub1] = convert2vector(minx_mu,minx_sig,lb.minx,ub.minx);
minx_gen = generate(mu,sig,repmat(lb1,[1,N]),repmat(ub1,[1,N]));
% generate upper bound
[mu,sig,lb1,ub1,minx_gen_bnd] = convert2vector(maxx_mu,maxx_sig,lb.maxx,ub.maxx,minxfns,minx_gen);
maxx_gen = generate(mu,sig,max(repmat(lb1,[1,N]),minx_gen_bnd),repmat(ub1,[1,N]));

% generate lower lower bound
[mu,sig,lb1,ub1,minx_gen_bnd] = convert2vector(minxclose_mu,minxclose_sig,lb.minxclose,ub.minxclose,minxfns,minx_gen);
minxclose_gen = generate(mu,sig,repmat(lb1,[1,N]),min(minx_gen_bnd,repmat(ub1,[1,N])));
% generate upper upper bound
[mu,sig,lb1,ub1,maxx_gen_bnd] = convert2vector(maxxclose_mu,maxxclose_sig,lb.maxxclose,ub.maxxclose,maxxfns,maxx_gen);
maxxclose_gen = generate(mu,sig,max(repmat(lb1,[1,N]),maxx_gen_bnd),repmat(ub1,[1,N]));

% generate lower bound on sum
[mu,sig,lb1,ub1] = convert2vector(minsumx_mu,minsumx_sig,lb.minsumx,ub.minsumx);
minsumx_gen = generate(mu,sig,repmat(lb1,[1,N]),repmat(ub1,[1,N]));
% generate upper bound on sum
[mu,sig,lb1,ub1,minsumx_gen_bnd] = convert2vector(maxsumx_mu,maxsumx_sig,lb.maxsumx,ub.maxsumx,minsumxfns,minsumx_gen);
maxsumx_gen = generate(mu,sig,max(repmat(lb1,[1,N]),minsumx_gen_bnd),repmat(ub1,[1,N]));

% generate lower bound on mean
[mu,sig,lb1,ub1] = convert2vector(minmeanx_mu,minmeanx_sig,lb.minmeanx,ub.minmeanx);
minmeanx_gen = generate(mu,sig,repmat(lb1,[1,N]),repmat(ub1,[1,N]));
% generate upper bound on mean
[mu,sig,lb1,ub1,minmeanx_gen_bnd] = convert2vector(maxmeanx_mu,maxmeanx_sig,lb.maxmeanx,ub.maxmeanx,minmeanxfns,minmeanx_gen);
maxmeanx_gen = generate(mu,sig,max(repmat(lb1,[1,N]),minmeanx_gen_bnd),repmat(ub1,[1,N]));


% closing radius
r = round(randn([1,N])*r_sig + r_mu);
for iter = 1:3,
  idx = r < minr | r > maxr;
  if ~any(idx),
    break;
  end
  r(idx) = round(randn([1,nnz(idx)])*r_sig + r_mu);
end
r(r < minr) = minr;
r(r > maxr) = maxr;

% when r is 0, equiv to having same cubes
idx = r == 0;

for i = 1:length(minxclosefns),
  fn = minxclosefns{i};
  j = find(strcmpi(fn,minxfns));
  if ~isempty(j),
    minxclose_gen(i,idx) = minx_gen(j,idx);
  end
end
for i = 1:length(maxxclosefns),
  fn = maxxclosefns{i};
  j = find(strcmpi(fn,maxxfns));
  if ~isempty(j),
    maxxclose_gen(i,idx) = maxx_gen(j,idx);
  end
end

x = [minx_gen;maxx_gen;minxclose_gen;maxxclose_gen;minsumx_gen;maxsumx_gen;minmeanx_gen;maxmeanx_gen;r];

if DEBUG,
  
  fprintf('\nCurrent Distribution:\n');
  allfns = union(minxfns,union(maxxfns,union(minxclosefns,union(maxxclosefns,union(minsumxfns,union(maxsumxfns,union(minmeanxfns,maxmeanxfns)))))));
  fprintf('variable       minxclose      minx           maxx           maxxclose      minsumx        maxsumx        minmeanx       maxmeanx\n');
  tmp_mu = {minxclose_mu,minx_mu,maxx_mu,maxxclose_mu,minsumx_mu,maxsumx_mu,minmeanx_mu,maxmeanx_mu};
  tmp_sig = {minxclose_sig,minx_sig,maxx_sig,maxxclose_sig,minsumx_sig,maxsumx_sig,minmeanx_sig,maxmeanx_sig};
  for i = 1:length(allfns),
    fn = allfns{i};
    l = length(fn);
    fprintf([fn,repmat(' ',[1,max(0,15-l)])]);
    if ~isempty(strfind(fn,'theta')) || ~isempty(strfind(fn,'angle')) || ...
        ~isempty(strfind(fn,'phi')) || ~isempty(strfind(fn,'yaw')),
      isangle = true;
    else
      isangle = false;
    end
    for j = 1:length(tmp_mu),
      mucurr = tmp_mu{j};
      sigcurr = tmp_sig{j};
      if isfield(mucurr,fn),
        if isangle,
          fprintf('%-5.2f+%5.2f    ',mucurr.(fn)*180/pi,sigcurr.(fn)*180/pi);
        else
          fprintf('%-5.2f+%5.2f    ',mucurr.(fn),sigcurr.(fn));
        end
      else
        fprintf('               ');
      end
    end
    fprintf('\n');
  end
  fprintf('r_mu = %f, r_sig = %f\n',r_mu,r_sig);
  
  fprintf('\nSample Distribution:\n');
  fprintf('variable       minxclose      minx           maxx           maxxclose      minsumx        maxsumx        minmeanx       maxmeanx\n');
  tmp = {minxclose_gen,minx_gen,maxx_gen,maxxclose_gen,minsumx_gen,maxsumx_gen,minmeanx_gen,maxmeanx_gen};
  tmpfns = {minxclosefns,minxfns,maxxfns,maxxclosefns,minsumxfns,maxsumxfns,minmeanxfns,maxmeanxfns};
  for i = 1:length(allfns),
    fn = allfns{i};
    l = length(fn);
    fprintf([fn,repmat(' ',[1,max(0,15-l)])]);
    if ~isempty(strfind(fn,'theta')) || ~isempty(strfind(fn,'angle')) || ...
        ~isempty(strfind(fn,'phi')) || ~isempty(strfind(fn,'yaw')),
      isangle = true;
    else
      isangle = false;
    end
    for j = 1:length(tmp),
      fnscurr = tmpfns{j};
      gencurr = tmp{j};
      k = find(strcmpi(fn,fnscurr),1);
      if ~isempty(k),
        gencurr = gencurr(k,:);
        mucurr = mean(gencurr);
        sigcurr = std(gencurr);
        if isangle,
          fprintf('%-5.2f+%5.2f    ',mucurr*180/pi,sigcurr*180/pi);
        else
          fprintf('%-5.2f+%5.2f    ',mucurr,sigcurr);
        end
      else
        fprintf('               ');
      end
    end
    fprintf('\n');
  end
  fprintf('r_mu = %f, r_sig = %f\n',mean(r),std(r));

end

function [mu,sig,lb,ub,varargout] = convert2vector(mus,sigs,lbs,ubs,varargin)

fns = fieldnames(mus);
n = length(fns);
mu = zeros(n,1);
sig = zeros(n,1);
lb = zeros(n,1);
ub = zeros(n,1);
m = length(varargin)/2;

for i = 1:m,
  varargout{i} = zeros(n,size(varargin{2*i},2));
end

for i = 1:n,
  fn = fns{i};
  mu(i) = mus.(fn);
  sig(i) = sigs.(fn);
  if isfield(lbs,fn),
    lb(i) = lbs.(fn);
  else
    lb(i) = -inf;
  end
  if isfield(ubs,fn),
    ub(i) = ubs.(fn);
  else
    ub(i) = inf;
  end
  for j = 1:m,
    fns2 = varargin{2*(j-1)+1};
    k = find(strcmpi(fn,fns2));
    if isempty(k),
      varargout{j}(i,:) = nan;
    else
      varargout{j}(i,:) = varargin{2*j}(k,:);
    end
  end
end

function x = generate(mu,sig,minx,maxx)

[n,N] = size(minx);
x = randn([n,N]).*repmat(sig,[1,N]) + repmat(mu,[1,N]);
for iter = 1:3,
  idx = x > maxx | x < minx;
  if ~any(idx(:)), break; end
  for i = 1:n,
    if ~any(idx(i,:)), continue; end
    x(i,idx(i,:)) = randn([1,nnz(idx(i,:))])*sig(i) + mu(i);
  end
end
idx = x > maxx;
if any(idx(:)),
  for i = 1:n,
    if ~any(idx(i,:)), continue; end
    x(i,idx(i,:)) = maxx(i,idx(i,:));
  end
end
idx = x < minx;
if any(idx(:)),
  for i = 1:n,
    if ~any(idx(i,:)), continue; end
    x(i,idx(i,:)) = minx(i,idx(i,:));
  end
end

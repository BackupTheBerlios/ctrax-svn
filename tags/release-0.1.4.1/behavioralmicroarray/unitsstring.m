function s = unitsstring(units)

s = '';
num = unique(units.num);
if isempty(num),
  s = '1';
else
  s = '( ';
  for i = 1:length(num),
    n = nnz(strcmpi(num{i},units.num));
    if n > 1,
      s = [s,num{i},'^',num2str(n),' '];
    else
      s = [s,num{i},' '];
    end
  end
  s = [s,')'];
end

den = unique(units.den);
if ~isempty(den),
  s = [s,' / ( '];
  for i = 1:length(den),
    n = nnz(strcmpi(den{i},units.den));
    if n > 1,
      s = [s,den{i},'^',num2str(n),' '];
    else
      s = [s,den{i},' '];
    end
  end
  s = [s,')'];
end
function v = isanglename(prop)

v = ~isempty(strfind(prop,'theta')) || ~isempty(strfind(prop,'angle')) || ...
  ~isempty(strfind(prop,'phi')) || ~isempty(strfind(prop,'yaw'));
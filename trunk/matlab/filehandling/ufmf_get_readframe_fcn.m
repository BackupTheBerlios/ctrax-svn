% readframe = ufmf_get_readframe_fcn(header)
%
% Returns the handle to a function which inputs the frame number and
% outputs the corresponding frame for a UFMF file described by the input
% 'header'. 
%
% Note that this 'header' is modified by ufmf_read_frame and that header
% seems to be stored in the function handle instance somehow. I don't know
% how fragile this is. 
%
function readframe = ufmf_get_readframe_fcn(header)

readframe = @(f) ufmfreadframe(f);

  function [im,timestamp,cachedidx,bb,mu] = ufmfreadframe(f)
    % this 'header' is stored in the function handle
    [im,header,timestamp,bb,mu] = ufmf_read_frame(header,f);
    cachedidx = header.cachedmeans_idx;
  end

end
% [sbfmfnames,rois] = avitosbfmf(...)
%
% Inputs an AVI movie and converts to an SBFMF movie.
% Also allows user to specify regions of interest, 
% in which case each region of interest is saved to a 
% different SBFMF. 
%
% Usage:
% The parameters to avitosbfmf can be specified as pairs 
% consisting of a string specifying the parameter name 
% followed by the parameter value. If a given parameter is 
% not specifed, then either the user is prompted for this
% value or default values are chosen.
%
% Parameters:
%
% 'aviname': Name of the AVI file to convert. If not 
%   specified, the user is prompted for a name. 
% 'sbfmfname': Name of the SBFMF file to output. If more
%   than one file is to be output, names are generated by 
%   appending '_box%d' to each name, e.g. movie.sbfmf will 
%   result in SBFMF files movie_box1.sbfmf, movie_box2.sbfmf.
%   If not specified, then the user is prompted for the name.
% 'nrois': Number of regions of interest. If not specified
%   but rois is specified, it is set to the number of rows in
%   rois. Otherwise, it is set to 0 and the whole frame is 
%   stored to one SBFMF. 
% 'rois': Regions of interest to crop and store. This should be a 
%   nrois x 4 matrix, where rois(i,:) is of the form 
%   [xmin,xmax,ymin,ymax]. If there are less rois specified than 
%   'nrois', more regions of interest are probpted. 
% 'startframe': first frame to write out to file. 1 by default.
% 'endframe': last frame to write out to file. nframes by default.
% 'differencemode': direction used in background subtraction. One of
%   'Dark flies on a light background',
%   'Light flies on a dark background',
%   'Other'
%   This is 'Other' by default.
% 'bgthresh': Background subtraction threshold for choosing which
%   pixels to store per frame. 1 by default. This should be between 1 
%   and 255.
% 'bgalgorithm': Algorithm to use for computing background center
%   and standard deviation. Options are 'median' and 'mean'. 'median'
%   by default.
% 'bgnframes': Number of frames sampled to estimate the background 
%   model. 200 by default.
% 'bgstartframe': start frame of interval used in background 
%   model estimation. 1 by default.
% 'bgendframe': end frame of interval used in background model 
%   estimation. nframes by default.
% 'maxmem': maximum amount of memory to use when computing the 
%   median backgroud model.
% 'verbose': how much diagnostic print statements should be executed.
%   0 means none, 1 means some, 2 means all. 1 by default.
%
% Outputs:
% sbfmfnames: Cell array containing the names of the SBFMF files saved. 
% rois: Regions of interest selected. 
%
% Example Usage:
%
% [sbfmfnames,rois] = avitosbfmf('ncropboxes',2,'bgthresh',3,...
%   'firstframe',1,'endframe',inf,'bgalgorithm','median')
%
% Relies On:
% compute_bg_model, compute_bg_mean, compute_bg_median,
% myparse, sbfmf_write_header, sbfmf_write_index, sbfmf_write_frame
%
% KB 01/13/2010
%
function [sbfmfnames,rois] = avitosbfmf(varargin)

version = '0.3b';
sbfmfnames = {};
rois = []; %#ok<NASGU>

% maxmem: maximum amount of memory to use when computing median

% parse inputs, set defaults as necessary
[aviname,sbfmfname,nrois,rois,startframe,endframe,...
  differencemode,bgthresh,bgalgorithm,bgnframes,bgstartframe,bgendframe,...
  maxmem,verbose] = ...
  myparse(varargin,...
  'aviname','','sbfmfname','','nrois',0,'rois',zeros(0,4),...
  'startframe',1,'endframe',inf,...
  'differencemode','Other','bgthresh',1,...
  'bgalgorithm','median','bgnframes',200,...
  'bgstartframe',1,'bgendframe',inf,...
  'maxmem',2^20*50,...
  'verbose',1);
differencemode = lower(differencemode);

% check for input errors
if startframe < 1 || endframe < startframe,
  error('1 <= startframe <= endframe must hold');
end

% set nrois if necessary
nroisin = size(rois,1);
nrois = max(nrois,nroisin);

% get filenames
if isempty(aviname),
  global AVITOSBFMF_OPENPATH; %#ok<TLEV>
  [aviname,avipath] = uigetfile('*.avi','Choose input AVI file',AVITOSBFMF_OPENPATH);
  if ~ischar(aviname), return; end
  AVITOSBFMF_OPENPATH = avipath;
  aviname = [avipath,aviname];
end
if isempty(sbfmfname),
  sbfmfname = strrep(aviname,'.avi','.sbfmf');
  [sbfmfname,sbfmfpath] = uiputfile(sbfmfname,'Choose output SBFMF file.');
  sbfmfname = [sbfmfpath,sbfmfname];
end

% open avi file for reading
readerobj = mmreader(aviname);
nframes = get(readerobj,'NumberOfFrames');
if isempty(nframes),
  % approximate nframes from duration
  nframes = get(readerobj,'Duration')*get(readerobj,'FrameRate');
end

% set endframes to be at most last frame
endframe = min(endframe,nframes);
bgendframe = min(bgendframe,nframes);

% create readframe function
readframe = @(f) flipdim(read(readerobj,f),1);
headerinfo = get(readerobj);
headerinfo.type = 'avi';
fps = headerinfo.FrameRate;

% read in the first frame to get the size
im = readframe(startframe);
nr = size(im,1);
nc = size(im,2);
ncolors = size(im,3);

if verbose >= 1,
  fprintf('AVI File Info:\n');
  fprintf('Name = %s\n',aviname);
  fprintf('Frame size = %d rows x %d columns x %d colors',nr,nc,ncolors);
  fprintf('Number of frames = %d\n',nframes);
  fprintf('Frame rate = %d\n',fps);
  fprintf('VideoFormat = %s\n',headerinfo.VideoFormat);
  fprintf('BitsPerPixel = %d\n',headerinfo.BitsPerPixel);
end

% set one roi for whole image of no ROIs specified
if nrois == 0,
  nrois = 1;
  rois = [1,nc,1,nr];
end

% make sure roi sizes are allowed, threshold if not
rois(:,1) = max(min(rois(:,1),nc),1);
rois(:,3) = max(min(rois(:,3),nr),1);
rois(:,2) = max(min(rois(:,2),nc),rois(:,1));
rois(:,4) = max(min(rois(:,4),nr),rois(:,3));

% get rois from user if not all specified
if nroisin < nrois,
  clf;
  fig = gcf;
  image(im);
  axis image;
  hold on;
  for i = 1:nroisin,
    plot(rois(i,[1,1,2,2,1]),rois(i,[3,4,4,3,3]),'r-');
  end
  rois = [rois;nan(nrois-nroisin,4)];
  for i = nroisin+1:nrois,
    title(sprintf('Select rectangular region of interest %d/%d',i,nrois));
    rect = getrect(fig);
    rois(i,:) = round([rect(1),rect(1)+rect(3),rect(2),rect(2)+rect(4)]);
    plot(rois(i,[1,1,2,2,1]),rois(i,[3,4,4,3,3]),'r-');
  end    

end

% make filenames for each roi
if nrois > 1,

  fprintf('ROIs will be saved to files:\n');
  sbfmfnames = cell(1,nrois);
  for i = 1:nrois,
    sbfmfnames{i} = strrep(sbfmfname,'.sbfmf',sprintf('_box%d.sbfmf',i));
    fprintf('%s\n',sbfmfnames{i});
  end  
else
  sbfmfnames = {sbfmfname};  
end

if verbose >= 1,
  fprintf('\nSaving SBFMFs from frame %d to %d\n',startframe,endframe);
  for i = 1:nrois,
    fprintf('ROI %d stored in %s\n',i,sbfmfnames{i});
    fprintf('ROI spans y = %d:%d, x = %d:%d, nrows = %d, ncolumns = %d\n',...
      rois(i,3),rois(i,4),...
      rois(i,1),rois(i,2),...
      rois(i,4)-rois(i,3)+1,...
      rois(i,2)-rois(i,1)+1);
  end
end


% compute bg model
bgcenter = cell(1,nrois);
bgstd = cell(1,nrois);
if verbose >= 1,
  fprintf('\nBackground modeling, subtraction parameters:\n');
  fprintf('Algorithm: %s\n',bgalgorithm);
  fprintf('%d frames sampled between frames %d and %d\n',bgnframes,bgstartframe,bgendframe);
  fprintf('Threshold = %f\n',bgthresh);
end
for i = 1:nrois,
  if verbose >= 1,
    fprintf('Computing background model for ROI %d/%d.\n',i,nrois);
  end
  [bgcenter{i},bgstd{i}] = compute_bg_model(readframe,nframes,...
                                            'bgalgorithm',bgalgorithm,...
                                            'bgstartframe',bgstartframe,...
                                            'bgendframe',bgendframe,...
                                            'bgnframes',bgnframes,...
                                            'roi',rois(i,:),...
                                            'maxmem',maxmem);
  bgcenter{i} = bgcenter{i}';
  bgstd{i} = bgstd{i}';
end

% open the sbfmfs
fp = zeros(1,nrois);
for i = 1:nrois,
  fp(i) = fopen(sbfmfnames{i},'w');
end

% where the index location will be stored
indexlocpr = nan(1,nrois);
% how many frames we are actually writing
nframeswrite = endframe - startframe + 1;
% write headers
for i = 1:nrois,
  % already taken the transpose
  indexlocpr(i) = sbfmf_write_header(fp(i),version,...
    nframeswrite,differencemode,bgcenter{i},bgstd{i});
end

% store the locations of the frames in the file
framesloc = nan(nrois,endframe-startframe+1);
% make stamps based on frame rate
stamp = (startframe-1)/fps;

if verbose >= 1,
  estsize = zeros(1,nrois);
end
for t = startframe:endframe,
    
  % read in the current frame
  im = readframe(t);
  
  % we'll need this for indexing into framesloc
  j = t - startframe + 1;
  
  if verbose >= 1 && mod(t,100) == 0,
    fprintf('Compressing frame %d of %d\n',j,endframe-startframe+1);
    if j > 1,
      for i = 1:nrois,
        fprintf('Average number of points stored for ROI %d = %d\n',...
                i,round(estsize(i)/(j-1)));
      end
    end
  end
  
  % write the current frame for each roi
  for i = 1:nrois,
    
    % convert to double gray scale, transpose
    im1 = double(rgb2gray(im(rois(i,3):rois(i,4),rois(i,1):rois(i,2),:)))';
    
    % subtract and threshold
    switch differencemode,
      case 'dark flies on a light background',
        idx = find(imsubtract(bgcenter{i},im1) >= bgthresh);
      case 'light flies on a dark backgroud',
        idx = find(imsubtract(im1,bgcenter{i}) >= bgthresh);
      otherwise,
        idx = find(imabsdiff(im1,bgcenter{i}) >= bgthresh);
    end

    if verbose >= 1,
      estsize(i) = estsize(i) + length(idx);
    end
    
    % store different pixels
    framesloc(i,j) = sbfmf_write_frame(fp(i),stamp,idx,im1(idx));
  end

  % increment the stamp
  stamp = stamp + 1/fps;
end

% write the indices
for i = 1:nrois,
  sbfmf_write_index(fp(i),framesloc(i,:),indexlocpr(i));
end

% close the files
for i = 1:nrois,
  fclose(fp(i));
end

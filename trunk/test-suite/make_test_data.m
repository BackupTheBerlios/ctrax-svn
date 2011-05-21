function make_test_data
% make_test_data
%
% makes Python-readable test data for Ctrax test suite
%
% JAB 5/12/11

% find directory for Python code -- same as own directory!
base_name = which( 'make_test_data' );
base_dir = base_name(1:strfind( base_name, 'make_test_data.m' ) - 1);

% use Python to make lists of filenames
cmd = ['!python ' base_dir 'dump_mat_names.py'];
fprintf( 1, [cmd '\n'] )
eval( cmd )

% read in list of original filenames
fp = fopen( 'mat_name_dump.tmp', 'r' );
tracked_filenames = {};
line = fgetl( fp );
while line ~= -1
   tracked_filenames{end+1} = line;
   line = fgetl( fp );
end
fclose( fp );
delete mat_name_dump.tmp

% resave files
for li = 1:length( tracked_filenames )
   a_ = {}; b_ = {}; x_ = {}; y_ = {}; theta_ = {};
   % load data
   load( tracked_filenames{li} )
   for an = 1:length( trx )
      a_{an} = trx(an).a;
      b_{an} = trx(an).b;
      theta_{an} = trx(an).theta;
      x_{an} = trx(an).x;
      y_{an} = trx(an).y;
   end
   
   % resave data in a Python-friendly MAT-file
   newname = [tracked_filenames{li}(1:strfind( tracked_filenames{li}, '.mat' ) - 1) '_forpython.mat'];
   fprintf( 1, ['writing ' newname '\n'] )
   save( newname, 'x_', 'y_', 'theta_', 'a_', 'b_' )
end

% read in list of newly tracked filenames
fp = fopen( 'mat_new_name_dump.tmp', 'r' );
new_filenames = {}; % need load_tracks() run on them
line = fgetl( fp );
while line ~= -1
   new_filenames{end+1} = line;
   line = fgetl( fp );
end
fclose( fp );
delete mat_new_name_dump.tmp

% resave files
for li = 1:length( new_filenames )
   a_ = {}; b_ = {}; x_ = {}; y_ = {}; theta_ = {};
   [trx, matname, succeeded] = load_tracks( new_filenames{li} );
   newname = [new_filenames{li}(1:strfind( new_filenames{li}, '.mat' ) - 1) '_forpython.mat'];
   if succeeded
      for an = 1:length( trx )
         a_{an} = trx(an).a;
         b_{an} = trx(an).b;
         theta_{an} = trx(an).theta;
         x_{an} = trx(an).x;
         y_{an} = trx(an).y;
      end
      
      % resave
      fprintf( 1, ['writing ' newname '\n'] )
      save( newname, 'x_', 'y_', 'theta_', 'a_', 'b_' )
   else
      fprintf( 1, '**failed loading tracks in %s\n', newname )
      eval( ['!rm ' newname] )
   end
end


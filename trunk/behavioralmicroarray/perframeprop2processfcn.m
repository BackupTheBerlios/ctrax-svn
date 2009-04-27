% [processfcns,stillmissing] = perframeprop2processfcn(perframeprops,[issocial=false])
function [processfcns,stillmissing] = perframeprop2processfcn(perframeprops,issocial)

if ~exist('issocial','var')
  issocial = false;
end

processdata_fns = {'dtheta','du_ctr','dv_ctr','corisonfly','corfrac_maj','corfrac_min',...
  'abscorfrac_min','du_cor','dv_cor','velmag_ctr','velmag','accmag','signdtheta',...
  'absdv_cor','flipdv_cor','absdtheta','d2theta','absd2theta','smooththeta','smoothdtheta',...
  'abssmoothdtheta','smoothd2theta','abssmoothd2theta','phi','yaw','absyaw'};

processdatacrabwalks_fns = {'du_tail','dv_tail','absdu_tail','absdv_tail','dtheta_tail',...
  'absdtheta_tail','phisideways','absphisideways'};

processdataarena_fns = {'dist2wall','theta2wall','ddist2wall','absdangle2wall'};

processdataclosest_fns = {'absdangle2wall','dcenter','closestfly_center','dnose2ell',...
  'closestfly_nose2ell','dell2nose','closestfly_ell2nose','anglesub','closestfly_anglesub',...
  'magveldiff_center','magveldiff_nose2ell','magveldiff_ell2nose','magveldiff_anglesub',...
  'veltoward_center','veltoward_nose2ell','veltoward_ell2nose','veltoward_anglesub',...
  'absthetadiff_center','absthetadiff_nose2ell','absthetadiff_ell2nose','absthetadiff_anglesub',...
  'absphidiff_center','absphidiff_nose2ell','absphidiff_ell2nose','absphidiff_anglesub',...
  'absanglefrom1to2_center','absanglefrom1to2_nose2ell','absanglefrom1to2_ell2nose',...
  'absanglefrom1to2_anglesub','ddcenter','ddnose2ell','ddell2nose','danglesub'};

processdatapair_fns = {'fly1','distnose2ell','magveldiff','dcenter','veltoward','thetadiff',...
  'absthetadiff','phidiff','absphidiff','minvelmag','maxvelmag','anglefrom1to2','absanglefrom1to2',...
  'minabsanglefrom1to2','anglesub','danglefrom1to2','danglesub'};

% if this is a social behavior, then remove the pair fieldnames from the
% closest fieldnames
if issocial,  
  processdataclosest_fns = setdiff(processdataclosest_fns,processdatapair_fns);
end

processfcns = {};
idx = ismember(perframeprops,processdata_fns);
if any(idx),
  processfcns{end+1} = 'process_data';
  perframeprops(idx) = [];
end
idx = ismember(perframeprops,processdatacrabwalks_fns);
if any(idx),
  processfcns{end+1} = 'process_data_crabwalks';
  perframeprops(idx) = [];
end
idx = ismember(perframeprops,processdataarena_fns);
if any(idx),
  processfcns{end+1} = 'process_data_arena';
  perframeprops(idx) = [];
end
idx = ismember(perframeprops,processdataclosest_fns);
if any(idx),
  processfcns{end+1} = 'process_data_closestfly';
  perframeprops(idx) = [];
end
idx = ismember(perframeprops,processdatapair_fns);
if any(idx),
  processfcns{end+1} = 'process_data_pairs';
  perframeprops(idx) = [];
end
stillmissing = perframeprops;


% out = anglesubtended(x1,y1,a1,b1,theta1,x2,y2,a2,b2,theta2,fov)
% inputs:
% (x1,y1,a1,b1,theta1) are vectors of length n1 describing the ellipses fit
% to the flies in set 1
% (x2,y2,a2,b2,theta2) are vectors of length n2 describing the ellipses fit
% to the flies in set 2
% fov is the field of view of the fly in radians
% output:
% out is a n1 x n2 matrix. out(i,j) is the angle of thet vision of fly i in
% set 1 subtended by fly j in set 2

function [angle_subtended,out1,out2] = tangentangles(x1,y1,a1,b1,theta1,x2,y2,a2,b2,theta2,debug)

% time_interval = t+[-t_length:t_length];
% x1=fly.x(fly1_i,time_interval);
% y1=fly.y(fly1_i,time_interval);
% theta1 = fly.theta(fly1_i,time_interval);
% a1=fly.maj_axis(fly1_i,time_interval);
% b1=fly.min_axis(fly1_i,time_interval);
% x2=fly.x(fly2_i,time_interval);
% y2=fly.y(fly2_i,time_interval);
% theta2 = fly.theta(fly2_i,time_interval);
% a2=fly.maj_axis(fly2_i,time_interval);
% b2=fly.min_axis(fly2_i,time_interval);


  if ~exist('debug','var'),
    debug = false;
  end;

  EPS = .00001;

  INTERVAL_LEN = length(x1);
  
  out1 = zeros(1,INTERVAL_LEN);
  out2 = zeros(1,INTERVAL_LEN);
  

  for frame_i = 1:INTERVAL_LEN,

      
      % corner case
      if max(a2(frame_i),b2(frame_i)) < EPS,
        xeye1 = x1(frame_i) + 2*a1(frame_i)*cos(theta1(frame_i));
        yeye1 = y1(frame_i) + 2*a1(frame_i)*sin(theta1(frame_i));
        out1(1,frame_i) = ataframe_i(x2(frame_i)-xeye1,y2(frame_i)-yeye1);
        out2(1,frame_i) = out1(1,frame_i);
        continue;
      end;
      
      % eye of fly 1
      [c,d] = eyeoffly1givenfly2([x1(frame_i),x2(frame_i)],[y1(frame_i),y2(frame_i)],...
                                 a1(frame_i),[theta1(frame_i),theta2(frame_i)]);

      % corner case: eye of fly 1 in fly 2
      [isonborder,isinborder] = checkinborder(c,d,a2(frame_i),b2(frame_i));
      if isonborder,
        % if eye is exactly on the border, then angle subtended is
        % at most pi
        out1(1,frame_i) = -pi/2;
        out2(1,frame_i) = pi/2;
        continue;
      elseif isinborder,
        % if inside border, then entire fov is subtended
        out1(1,frame_i) = -pi;
        out2(1,frame_i) = pi;
        continue;
      end;
      
      [psi1,psi2] = computetangentpoints(c,d,a2(frame_i),b2(frame_i),[theta1(frame_i),theta2(frame_i)]);

      out1(1,frame_i) = psi1;
      out2(1,frame_i) = psi2;
      
      if debug,
        
        c1 = x1(frame_i) + a1(frame_i)*cos(theta1(frame_i));
        d1 = y1(frame_i) + a1(frame_i)*sin(theta1(frame_i));
        s = 1.5*sqrt((x1(frame_i)-x2(frame_i))^2+(y1(frame_i)-y2(frame_i))^2);
        
        figure(333); clf;
        colors = 'rg';
        hold on;
        ellipsedraw(a1(frame_i),b1(frame_i),x1(frame_i),y1(frame_i),theta1(frame_i),colors(1));
        ellipsedraw(a2(frame_i),b2(frame_i),x2(frame_i),y2(frame_i),theta2(frame_i),colors(2));
        plot(c1,d1,[colors(1),'o']);
        e1 = c1+cos(psi1+theta1(frame_i))*s;
        f1 = d1+sin(psi1+theta1(frame_i))*s;
        hpsi1 = plot([c1,e1],[d1,f1],'m')
        text(e1+1,f1,'psi1');
        e1 = c1+cos(psi2+theta1(frame_i))*s;
        f1 = d1+sin(psi2+theta1(frame_i))*s;
        hpsi2 = plot([c1,e1],[d1,f1],'m')
        text(e1+1,f1,'psi2');
        axis equal
        axis ij

        input(sprintf('psi1 = %f, psi2 = %f ',psi1,psi2));
      end;
      
    end;
    
    out1 = -out1; %% Because angles are reversed here
    out2 = -out2;
    angle_subtended = abs(modrange(out1-out2,-pi,pi));
  
  
function dpsi = limitbyfov(psi1,psi2,fov)

  % bounds on field of view
  fov1 = -fov/2;
  % make +fov/2 be after fov1
  fov2 = fov1+fov;
  
  % make psi1, psi2 be after fov1
  d = mod(psi1-fov1,2*pi);
  psi1 = fov1 + d;
  d = mod(psi2-fov1,2*pi);
  psi2 = fov1 + d;

  if (fov2 <= psi1) && (psi1 <= psi2),
    dpsi = 0;
  elseif (fov2 <= psi2) && (psi2 <= psi1),
    dpsi = mod(fov1 - psi1,2*pi) + mod(fov2 - psi2,2*pi);
  elseif (psi1 <= fov2) && (fov2 <= psi2),
    dpsi = fov2 - psi1;
  elseif (psi1 <= psi2) && (psi2 <= fov2),
    dpsi = psi2 - psi1;
  elseif (psi2 <= fov2) && (fov2 <= psi1),
    dpsi = psi2 - fov1;
  else, % psi2 <= psi1 <= fov2
    dpsi = (psi2 - fov1) + (fov2 - psi1);
  end;
    
function [psi1,psi2] = computetangentpoints(c,d,a,b,theta)
                     
  EPS = .00001;

  % solve for phi, the angles of ellipse 2 at which the tangent line
  % passes through (c,d)
  
  possiblephi = zeros(6,1);
  cost = zeros(6,1);
  
  % solve quadratic equation for cos(phi), where phi is the angle
  % of ellipse 2 at the intersection point
  % equation of the ellipse:
  % x = a * cos(phi)
  % y = b * sin(phi)
  % slope of tangent line at phi:
  % dy/dx = (dy/dphi) / (dx/dphi)
  % = (b*cos(phi))/(-a*sin(phi))
  % slope can also be determined from vector from (c,d) to
  % intersection point:
  % dy/dx = (b*sin(phi) - d) / (a*cos(phi) - c)
  % solve for cos(phi) in equation:
  % (b*cos(phi))/(-a*sin(phi)) = (b*sin(phi)-d)/(a*cos(phi)-c)
  % yields A = a^2*d^2 + b^2*c^2
  % B = - 2*a*b^2*c
  % C = a^2*(b^2-d^2)
  % A*(cos(phi))^2 + B*cos(phi) + C = 0
  A = b^2*c^2 + a^2*d^2;
  B = -2*a*b^2*c;
  C = a^2*(b^2 - d^2);
  D = B^2 - 4*A*C;
  
  % if D < 0, then
  % a^2*b^4*c^2 < (b^2*c^2 + a^2*d^2)*a^2*(b^2-d^2)
  % b^4*c^2 < b^4*c^2 + a^2*b^2*d^2 - b^2*c^2*d^2 - a^2d^4
  % 0 < a^2*b^2*d^2 - b^2*c^2*d^2 - a^2d^4
  % 0 < a^2*b^2 - b^2*c^2 - a^2d^2
  % 1 > c^2/a^2 + d^2/b^2
  % in which case (c,d) is inside the ellipse
  % this case has already been checked for
  D = sqrt(D);
  
  % if sin(phi) = 0
  % then the slope is undefined, so we will consider these as special
  % cases
  possiblephi(1) = 0;
  possiblephi(2) = pi;
  % phi = 0 only if c = a
  if abs(c - a) < EPS,
    cost(1) = 0;
  else,
    cost(1) = inf;
  end;
  % phi = pi only if c = -a
  if abs(c + a) < EPS,
    cost(2) = 0;
  else,
    cost(2) = inf;
  end;
  
  % four possible choices for phi
  cosphi = (-B + D)/(2*A);
  possiblephi(3) = acos(cosphi);
  possiblephi(4) = -possiblephi(3);
  
  cosphi = (-B - D)/(2*A);
  possiblephi(5) = acos(cosphi);
  possiblephi(6) = -possiblephi(5);
  
  sinphi = sin(possiblephi(3:end));
  cosphi = cos(possiblephi(3:end));
  cost(3:end) = abs((b*sinphi - d).*(-a*sinphi) - ...
    (a*cosphi - c).*(b*cosphi));
  
  [sortedcost,order] = sort(cost);
  phi1 = possiblephi(order(1));
  phi2 = possiblephi(order(2));
  
  x1 = a*cos(phi1);
  y1 = b*sin(phi1);
  x2 = a*cos(phi2);
  y2 = b*sin(phi2);
  
  % compute angles from (c,d) to (x1,y1), to (x2,y2)
  psi1 = atan2(y1 - d, x1 - c);
  psi2 = atan2(y2 - d, x2 - c);
    
  % compute angle from (c,d) to (0,0)
  psi0 = atan2(-d,-c);
  
  % make relative to gaze direction of fly 1
  psi0 = psi0 + theta(2) - theta(1);
  psi1 = psi1 + theta(2) - theta(1);
  psi2 = psi2 + theta(2) - theta(1);
  
  % put psi1 in [-pi,pi)
  psi1 = mod(psi1+pi,2*pi)-pi;

  % make psi0 and psi2 after psi1
  dpsi01 = mod(psi0-psi1,2*pi);
  psi0 = psi1 + dpsi01;
  dpsi21 = mod(psi2-psi1,2*pi);
  psi2 = psi1 + dpsi21;
  
  % if psi2 comes before psi0, swap
  if psi2 < psi0,
    
    tmp = psi1;
    psi1 = psi2;
    psi2 = tmp;

    % put psi1 in [-pi,pi)
    psi1 = mod(psi1+pi,2*pi)-pi;
    % make psi2 after psi1
    dpsi21 = mod(psi2-psi1,2*pi);
    psi2 = psi1 + dpsi21;
    
  end;
  
function [c,d] = eyeoffly1givenfly2(x,y,a,theta);

  % coordinate of eye of fly 1
  c1 = x(1) + a*cos(theta(1));
  d1 = y(1) + a*sin(theta(1));
  
  % translate and rotate coordinate system to be centered on
  % fly 2, aligned with fly 2
  c1 = c1 - x(2);
  d1 = d1 - y(2);
  c = c1*cos(theta(2)) + d1*sin(theta(2));
  d = d1*cos(theta(2)) - c1*sin(theta(2));
  
function [isonborder,isinborder] = checkinborder(c,d,a,b);
  
% check to see if the eye of fly i1 is inside ellipse i2
% since (c,d) is aligned with fly 2, just need to check to see 
% if (c,d) is inside the ellipse centered at (0,0), aligned
% with axes, and with x-axis length a[i2] and y-axis length
% b[i2]. 
% the equation of this ellipse is x^2/a^2 + y^2/b^2 <= 1.

  EPS = .00001;

  A = c^2/a^2 + d^2/b^2;
  if abs(A - 1) < EPS,
    % if eye is exactly on the border, then angle subtended is
    % at most pi
    isonborder = true;
    isinborder = false;
  elseif A < 1,
    isonborder = false;
    isinborder = true;
  else,
    isonborder = false;
    isinborder = false;
  end;

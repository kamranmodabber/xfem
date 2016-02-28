% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************
function [f,dfdx,dfdy] = branch(r,theta,alpha)
% Compute the branch functions spanning the near tip field for LEFM
% Inputs: 
%   (r,theta) : polar coordinates of points where the branch
%               functions are to be evaluated
%   alpha     : inclination of the crack tip segment w.r.t x axis

if( r ~=0 )
    r2 = sqrt(r);
else
    r2    = 0.1d-4;
    theta = 0.0d0 ;
end

fac  = 0.5/r2 ;
st2  = sin(theta/2.);
ct2  = cos(theta/2.);
st   = sin(theta);
ct   = cos(theta);

drdx = cos(alpha);
drdy = sin(alpha);
dtdx = -sin(alpha);
dtdy = cos(alpha);


% Functions 

f(1) = r2 * st2 ;
f(2) = r2 * ct2;
f(3) = r2 * st2 * ct;
f(4) = r2 * ct2 * ct;

% Derivatives

% first function
dPhidr  = -fac * st2;
dPhidt  =  fac * ct2;
dfdx(1) = dPhidr * drdx + dPhidt * dtdx  ;
dfdy(1) = dPhidr * drdy + dPhidt * dtdy  ;

% second function
dPhidr  = fac * ct2;
dPhidt  = fac * st2;
dfdx(2) = dPhidr * drdx + dPhidt * dtdx  ;
dfdy(2) = dPhidr * drdy + dPhidt * dtdy  ;

% third function
dPhidr  = fac * st2 * (2*st*st-ct) ;
dPhidt  = fac * ct* cos(3*theta/2.);
dfdx(3) = dPhidr * drdx + dPhidt * dtdx  ;
dfdy(3) = dPhidr * drdy + dPhidt * dtdy  ;

% fourth function
dPhidr  =  fac * ct2 * (ct + 2*st*st) ;
dPhidt  = -fac * ct * sin(3*theta/2.);
dfdx(4) = dPhidr * drdx + dPhidt * dtdx  ;
dfdy(4) = dPhidr * drdy + dPhidt * dtdy  ;
        

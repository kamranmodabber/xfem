% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************

function [ux,uy] = exact_Griffith(x,E,nu,stressState,sigmato,xTip,adv,cracklength)
% Compute the exact displacement for the infinite plate with centered crack
% Inputs:
%     - x(1,2) : coordinates of point where exact displacements are to be
%                evaluated
%     - E,nu, stressState: material properties
%     - sigmato : the loading
%     - xTip,adv,cracklength: crack data 

% inclination of local coord
alfa = atan2(adv(2),adv(1));
% transpose of the rotation matrix
QT = [cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
% local coordinates
xp = QT*(x-xTip)';         
% local polar coordinates
[theta,r] = cart2pol(xp(1),xp(2));

sigma = sigmato;
KI    = sigmato*sqrt(pi*cracklength);   % exact KI

mu = E/(2*(1+nu));

if ( strcmp(stressState,'PLANE_STRAIN') )
    ux = 2*(1+nu)*KI/E*sqrt(r/2/pi)*cos(theta/2)*(2-2*nu-cos(theta/2)^2);
    uy = 2*(1+nu)*KI/E*sqrt(r/2/pi)*sin(theta/2)*(2-2*nu-cos(theta/2)^2);
end





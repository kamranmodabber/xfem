% **********************************************************
%       TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                 Nguyen Vinh Phu
%            LTDS, ENISE, Juillet 2006
% *********************************************************

function d = signed_distance1(Qq,pt)

 x   = pt(1,1);
y   = pt(1,2);
for j=1:size(Qq,1)-1

    if x>=Qq(j,1) & x<=Qq(j+1,1)
       shomare_node=j;
    end

    if x<Qq(1,1)
       shomare_node=1;
    end
       if x>Qq(end,1)
       shomare_node=size(Qq,1)-1;
       end
   
end


   x0=Qq(shomare_node,1);    y0=Qq(shomare_node,2);
   x1=Qq(shomare_node+1,1);    y1=Qq(shomare_node+1,2);


l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0) ;
d   = phi/l;            


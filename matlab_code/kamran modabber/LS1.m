function [LS] = LS1(Qq,element,numelem,TYPE,node)
LS=zeros(4,numelem);
%shomare_node=zeros(numelem,4);
for i=1:numelem
a=TYPE(i,3);
if a~=0
for j=1:4

   x0=Qq(a,1);    y0=Qq(a,2);
   x1=Qq(a+1,1);    y1=Qq(a+1,2);

    x = node(element(i,j),1);
    y = node(element(i,j),2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
 LS(j,i) = phi/l; 
end
end
end

















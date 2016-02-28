function [TYPE,noke_tarak,noke_gere,noke_tarak_elemani]=MASIR(vorodi1,vorodi2,tedad,Q)
TYPE=zeros(tedad,3);
tip=2;
for j=1:size(Q,1)-1
    
xxtip1=Q(j,1);
yytip1=Q(j,2);
xxtip2=Q(j+1,1);
yytip2=Q(j+1,2);

a=(yytip2-yytip1);
b=-(xxtip2-xxtip1);
c=-xxtip1*a+yytip1*(-b);

for i=1:tedad
%for i=313:313
    TYPE(i,1)=i;
    x1=vorodi1(vorodi2(i,1),1);
    y1=vorodi1(vorodi2(i,1),2);
    x2=vorodi1(vorodi2(i,2),1);
    y2=vorodi1(vorodi2(i,2),2);
    x3=vorodi1(vorodi2(i,3),1);
    y3=vorodi1(vorodi2(i,3),2);
    x4=vorodi1(vorodi2(i,4),1);
    y4=vorodi1(vorodi2(i,4),2);   
    
    % yaftan noghat ghate tarak ba labe eleman
    a1=(y2-y1);
    b1=-(x2-x1);                 
    c1=-x1*a1+y1*(-b1);
    A11=[a1 b1;a b];
    A12=-[c1;c];
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    a2=(y4-y1);
    b2=-(x4-x1);
    c2=-x1*a2+y1*(-b2);
    A21=[a2 b2;a b];
    A22=-[c2;c];
    %%%%%%%%%%%%%%%%%%%%%
    a3=(y3-y2);
    b3=-(x3-x2);
    c3=-x2*a3+y2*(-b3);
    A31=[a3 b3;a b];
    A32=-[c3;c];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   a4=(y4-y3);
    b4=-(x4-x3);
    c4=-x2*a4+y2*(-b4);
    A41=[a4 b4;a b];
    A42=-[c4;c];
%%%%%%%%%%%%%%%%%%%%%


    
    
    % yaftan mahal ghate
    if det(A11)~=0
    j1=((A11)^-1)*A12;
    if j1(1,1)>=min(Q(j,1),Q(j+1,1)) & j1(1,1)<=max(Q(j,1),Q(j+1,1)) & j1(2,1)>=min(Q(j,2),Q(j+1,2)) & j1(2,1)<=max(Q(j,2),Q(j+1,2))
  
    if j1(1,1)+10^(-10)>=min(x1,x2) & j1(1,1)<=max(x1,x2)+10^(-10) & j1(2,1)+10^(-10)>=min(y1,y2) & j1(2,1)<=max(y1,y2)+10^(-10)
 
        TYPE(i,2)=TYPE(i,2)+1;
        if TYPE(i,3)==j-1
          if ((x1+x2+x3+x4)/4)<Q(j,1)
              TYPE(i,3)=j-1;
           
        else
        TYPE(i,3)=j;
          end
        else
        TYPE(i,3)=j;    
     end
    end
    end
    end
    
    
    
    
  if det(A21)~=0
    j2=((A21)^-1)*A22;
     if j2(1,1)>=min(Q(j,1),Q(j+1,1)) & j2(1,1)<=max(Q(j,1),Q(j+1,1)) & j2(2,1)>=min(Q(j,2),Q(j+1,2)) & j2(2,1)<=max(Q(j,2),Q(j+1,2))
    if j2(1,1)+10^(-10)>=min(x1,x4) & j2(1,1)<=max(x1,x4)+10^(-10) & j2(2,1)+10^(-10)>=min(y1,y4) & j2(2,1)<=max(y1,y4)+10^(-10)
     
     TYPE(i,2)=TYPE(i,2)+1;
            if TYPE(i,3)==j-1
          if ((x1+x2+x3+x4)/4)<Q(j,1)
              TYPE(i,3)=j-1;
           
        else
        TYPE(i,3)=j;
          end
        else
        TYPE(i,3)=j;    
     end
           
           
    end
    end 
       end
    
    if det(A31)~=0
    j3=((A31)^-1)*A32;
     if j3(1,1)>=min(Q(j,1),Q(j+1,1)) & j3(1,1)<=max(Q(j,1),Q(j+1,1)) & j3(2,1)>=min(Q(j,2),Q(j+1,2)) & j3(2,1)<=max(Q(j,2),Q(j+1,2))
    if j3(1,1)+10^(-10)>=min(x2,x3) & j3(1,1)<=max(x2,x3)+10^(-10) & j3(2,1)+10^(-10)>=min(y2,y3) & j3(2,1)<=max(y2,y3)+10^(-10)

     TYPE(i,2)=TYPE(i,2)+1;
            if TYPE(i,3)==j-1
          if ((x1+x2+x3+x4)/4)<Q(j,1)
              TYPE(i,3)=j-1;
           
        else
        TYPE(i,3)=j;
          end
        else
        TYPE(i,3)=j;    
     end
    end
     end
    end
    
        if det(A41)~=0
    j4=((A41)^-1)*A42;
     if j4(1,1)>=min(Q(j,1),Q(j+1,1)) & j4(1,1)<=max(Q(j,1),Q(j+1,1)) & j4(2,1)>=min(Q(j,2),Q(j+1,2)) & j4(2,1)<=max(Q(j,2),Q(j+1,2))
    if j4(1,1)+10^(-10)>=min(x3,x4) & j4(1,1)<=max(x3,x4)+10^(-10) & j4(2,1)+10^(-10)>=min(y3,y4) & j4(2,1)<=max(y3,y4)+10^(-10)
      
     TYPE(i,2)=TYPE(i,2)+1;
            if TYPE(i,3)==j-1
          if ((x1+x2+x3+x4)/4)<Q(j,1)
              TYPE(i,3)=j-1;
           
        else
        TYPE(i,3)=j;
          end
        else
        TYPE(i,3)=j;    
     end
    end
     end
     end 
    
    
end

end

for i=1:tedad
if TYPE(i,2)~=0
TYPE(i,2)=2;
end
end


xtip1=Q(1,1);
ytip1=Q(1,2);
xtip2=Q(end,1);
ytip2=Q(end,2);

noke_tarak=zeros(2,3);
noke_tarak_elemani=zeros(tedad,3);
noke_tarak_elemani(1:tedad,1)=1:tedad;

for i=1:tedad
    x1=vorodi1(vorodi2(i,1),1);
    y1=vorodi1(vorodi2(i,1),2);
    x2=vorodi1(vorodi2(i,2),1);
    y2=vorodi1(vorodi2(i,2),2);
    x3=vorodi1(vorodi2(i,3),1);
    y3=vorodi1(vorodi2(i,3),2);
    x4=vorodi1(vorodi2(i,4),1);
    y4=vorodi1(vorodi2(i,4),2); 

    
    masaht_kol=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))+abs(x1*(y3-y4)+x3*(y4-y1)+x4*(y1-y3));
        masaht_joz2=abs(xtip2*(y2-y3)+x2*(y3-ytip2)+x3*(ytip2-y2))+abs(x1*(ytip2-y4)+xtip2*(y4-y1)+x4*(y1-ytip2))+abs(x1*(y2-ytip2)+x2*(ytip2-y1)+xtip2*(y1-y2))+abs(x3*(y4-ytip2)+x4*(ytip2-y3)+xtip2*(y3-y4));
    if masaht_kol-masaht_joz2>=-10^(-13) & masaht_kol-masaht_joz2<=10^(-13)
    TYPE(i,2)=1;
    TYPE(i,3)=size(Q,1)-1;
    noke_tarak(2,1)=i;
    noke_tarak(2,2)=xtip2;
    noke_tarak(2,3)=ytip2;
    noke_tarak_elemani(i,2)=xtip2;
    noke_tarak_elemani(i,3)=ytip2;
    
    
    end
end
if tip==2
for i=1:tedad
    x1=vorodi1(vorodi2(i,1),1);
    y1=vorodi1(vorodi2(i,1),2);
    x2=vorodi1(vorodi2(i,2),1);
    y2=vorodi1(vorodi2(i,2),2);
    x3=vorodi1(vorodi2(i,3),1);
    y3=vorodi1(vorodi2(i,3),2);
    x4=vorodi1(vorodi2(i,4),1);
    y4=vorodi1(vorodi2(i,4),2); 

    
    masaht_kol=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))+abs(x1*(y3-y4)+x3*(y4-y1)+x4*(y1-y3));
        masaht_joz1=abs(xtip1*(y2-y3)+x2*(y3-ytip1)+x3*(ytip1-y2))+abs(x1*(ytip1-y4)+xtip1*(y4-y1)+x4*(y1-ytip1))+abs(x1*(y2-ytip1)+x2*(ytip1-y1)+xtip1*(y1-y2))+abs(x3*(y4-ytip1)+x4*(ytip1-y3)+xtip1*(y3-y4));
   if masaht_kol-masaht_joz1>=-10^(-13) & masaht_kol-masaht_joz1<=10^(-13)
    TYPE(i,2)=1;
    TYPE(i,3)=1;
    noke_tarak(1,1)=i;
    noke_tarak(1,2)=xtip1;
    noke_tarak(1,3)=ytip1;
    noke_tarak_elemani(i,2)=xtip1;
    noke_tarak_elemani(i,3)=ytip1;
    end
end

end
noke_gere(1:size(vorodi1,1),1)=1:size(vorodi1,1);
for i=1:size(vorodi1,1)

x=vorodi1(i,1);
y=vorodi1(i,2);

d1=(x-xtip1)^2+(y-ytip1)^2;
d2=(x-xtip2)^2+(y-ytip2)^2;
if d1<=d2
  noke_gere(i,2)=1;
else
     noke_gere(i,2)=2;
end
end
z=TYPE;






























function[darsad_khata]=adaptive2(numnode,numelem,vorodi1,vorodi2,stress)

TEDGERE=zeros(numnode,2);
TEDGERE(1:numnode,1)=1:numnode;
for i=1:numnode
for j=1:numelem
    for k=1:4
if vorodi2(j,k)==i
    TEDGERE(i,2)=TEDGERE(i,2)+1;
end
    end
end
end
 AS=zeros(numnode,3);
for i=1:numelem
  for j=1:4
AS(vorodi2(i,j),1)=AS(vorodi2(i,j),1)+stress(i,j,1);
AS(vorodi2(i,j),2)=AS(vorodi2(i,j),2)+stress(i,j,2);
AS(vorodi2(i,j),3)=AS(vorodi2(i,j),3)+stress(i,j,3);

  end
end
for i=1:numnode
AS(i,1)=AS(i,1)/TEDGERE(i,2);
AS(i,2)=AS(i,2)/TEDGERE(i,2);
AS(i,3)=AS(i,3)/TEDGERE(i,2);
end
Q=[-0.5777 -0.5777;0.5777 -0.5777;0.5777 0.5777;-0.5777 0.5777];
%Q=[-1 -1;1 -1;1 1;-1 1];


for i=1:numelem
    for j=1:4
 gpnt=Q(j,:);
 [N,dNdxi]=lagrange_basis('Q4',gpnt);
    tanesh_behbood(i,j,1)=N(1)*AS(vorodi2(i,1),1)+N(2)*AS(vorodi2(i,2),1)+N(3)*AS(vorodi2(i,3),1)+N(4)*AS(vorodi2(i,4),1);
    tanesh_behbood(i,j,2)=N(1)*AS(vorodi2(i,1),2)+N(2)*AS(vorodi2(i,2),2)+N(3)*AS(vorodi2(i,3),2)+N(4)*AS(vorodi2(i,4),2);
    tanesh_behbood(i,j,3)=N(1)*AS(vorodi2(i,1),3)+N(2)*AS(vorodi2(i,2),3)+N(3)*AS(vorodi2(i,3),3)+N(4)*AS(vorodi2(i,4),3);
    end
end
khata=zeros(numelem,4);

for i=1:numelem
for j=1:4
khata(i,j)=((tanesh_behbood(i,j,1)-stress(i,j,1))^2+(tanesh_behbood(i,j,2)-stress(i,j,2))^2+(tanesh_behbood(i,j,3)-stress(i,j,3))^2)^0.5;
end
end
for i=1:numelem
    for j=1:4
        
      if  khata(i,j)==0
        
      end
    end
    
end

khata_gere=zeros(numnode,1);
for i=1:numelem
  for j=1:4
khata_gere(vorodi2(i,j),1)=khata_gere(vorodi2(i,j),1)+khata(i,j);
  end
end
for i=1:numnode
khata_gere(i,1)=khata_gere(i,1)/TEDGERE(i,2);
end
help=khata_gere(:,1);
help1=find(help == 0);
    help(help1(1:end,1),:)=[];
MIN=min(help);
for i=1:size(help1,1)
       khata_gere(help1(i,1),1)=MIN; 
end

 NORMTANESH=zeros(numelem,4);
TANESHKOL=0;
KHATAKOL=0;
for i=1:numelem
    x1=vorodi1(vorodi2(i,1),1);
    y1=vorodi1(vorodi2(i,1),2);
    x2=vorodi1(vorodi2(i,2),1);
    y2=vorodi1(vorodi2(i,2),2);
    x3=vorodi1(vorodi2(i,3),1);
    y3=vorodi1(vorodi2(i,3),2);
    x4=vorodi1(vorodi2(i,4),1);
    y4=vorodi1(vorodi2(i,4),2); 
    masaht_kol=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))+abs(x1*(y3-y4)+x3*(y4-y1)+x4*(y1-y3));
    for j=1:4
 NORMTANESH(i,j)=(stress(i,j,1)^2+stress(i,j,2)^2+stress(i,j,3)^2)^0.5;
 TANESHKOL=TANESHKOL+NORMTANESH(i,j)*masaht_kol*0.5;
 KHATAKOL=KHATAKOL+ khata(i,j)*masaht_kol*0.5;
 
    end
end

%darsad_khata=(KHATAKOL/TANESHKOL)*100;
chegaliold=zeros(numnode,1);
 shomarande=zeros(numnode,1);
 
 for i=1:numelem
     
     a=vorodi2(i,1);
    b=vorodi2(i,2);
    c=vorodi2(i,3);
    d=vorodi2(i,4);
    chegaliold(a,1)=chegaliold(a,1)+((vorodi1(a,1)-vorodi1(b,1))^2+(vorodi1(a,2)-vorodi1(b,2))^2)^0.5;
    shomarande(a,1)=shomarande(a,1)+1;
     chegaliold(a,1)=chegaliold(a,1)+((vorodi1(a,1)-vorodi1(d,1))^2+(vorodi1(a,2)-vorodi1(d,2))^2)^0.5;
    shomarande(a,1)=shomarande(a,1)+1;
    
    chegaliold(b,1)=chegaliold(b,1)+((vorodi1(b,1)-vorodi1(a,1))^2+(vorodi1(b,2)-vorodi1(a,2))^2)^0.5;
    shomarande(b,1)=shomarande(b,1)+1;
    chegaliold(b,1)=chegaliold(b,1)+((vorodi1(b,1)-vorodi1(c,1))^2+(vorodi1(b,2)-vorodi1(c,2))^2)^0.5;
    shomarande(b,1)=shomarande(b,1)+1;
    
     chegaliold(c,1)=chegaliold(c,1)+((vorodi1(c,1)-vorodi1(b,1))^2+(vorodi1(c,2)-vorodi1(b,2))^2)^0.5;
    shomarande(c,1)=shomarande(c,1)+1;
    chegaliold(c,1)=chegaliold(c,1)+((vorodi1(c,1)-vorodi1(d,1))^2+(vorodi1(c,2)-vorodi1(d,2))^2)^0.5;
    shomarande(c,1)=shomarande(c,1)+1;
    
     chegaliold(d,1)=chegaliold(d,1)+((vorodi1(d,1)-vorodi1(c,1))^2+(vorodi1(d,2)-vorodi1(c,2))^2)^0.5;
    shomarande(d,1)=shomarande(d,1)+1;
    chegaliold(d,1)=chegaliold(d,1)+((vorodi1(d,1)-vorodi1(a,1))^2+(vorodi1(d,2)-vorodi1(a,2))^2)^0.5;
    shomarande(d,1)=shomarande(d,1)+1;
    
    
     end

for i=1:numnode
chegaliold1(i,1)=chegaliold(i,1)/shomarande(i,1);
end
chegaliold1;


normmiyangin=0;
for i=1:numnode
    normtanesh(i,1)=((AS(i,1))^2+(AS(i,2))^2+(AS(i,3))^2)^0.5;
normmiyangin=normmiyangin+normtanesh(i,1);
end
normmiyangin=normmiyangin/numnode;

%chegaliold1=chegali_gid;

for i=1:numnode
    
    eror(i,1)=((khata_gere(i,1))/normmiyangin)*100;
    KHATAHADAF(i,1)=0.5*normmiyangin;
    chegalinew(i,1)=((KHATAHADAF(i,1))/(khata_gere(i,1)))*chegaliold1(i,1)*1.78;
    %chegalinew(i,1)=((KHATAHADAF(i,1))/(khata_gere(i,1)))*1.5;
end

%had1=min(chegaliold1);
had2=max(chegaliold1);
had1=0.1;
%had2=2;

for i=1:numnode
   if chegalinew(i,1)<=had1
chegalinew(i,1)=had1;
    end
    if chegalinew(i,1)>=had2;  
   chegalinew(i,1)=had2;
    end
end
chegalinew ;

mod=fopen('gid.bgm', 'wt');
fprintf(mod, 'BackgroundMesh V 1.0\n');
fprintf(mod, 'MESH    Dimension  2 ElemType Quadrilateral  Nnode 4\n');
fprintf(mod, 'Coordinates\n');
for i=1:numnode
    fprintf(mod, '%6.0f %12.8f %12.8f\n',i,vorodi1(i,1),vorodi1(i,2));
end
fprintf(mod, 'End Coordinates\n');
fprintf(mod, 'Elements\n');
for i=1:numelem
    fprintf(mod, '%6.0f %6.0f %6.0f %6.0f %6.0f\n',i,vorodi2(i,1),vorodi2(i,2),vorodi2(i,3),vorodi2(i,4));
end
fprintf(mod, 'end elements\n');
fprintf(mod, 'DesiredSize Nodes\n');
for i=1:numnode
   %fprintf(mod, '%6.0f %12.8f\n',i,chegalinew(i,1)); 
   %fprintf(mod, '%6.0f %12.8f\n',i,((KHATAHADAF(i,1))/(khata_gere(i,1)))); 
   %fprintf(mod, '%6.0f %12.8f\n',i,2);
  fprintf(mod, '%6.0f %12.8f\n',i,chegalinew(i,1));
end
fprintf(mod, 'End DesiredSize\n');

fid = fopen('xfem.dat', 'wt');
fprintf(fid, 'TITLE="HAMID MOSLEMI"\n');
fprintf(fid, 'VARIABLES = "X","Y","SX","SY","SXY","EROR"\n');
fprintf(fid, 'ZONE N= %6.0f ,E= %6.0f ,F=FEPOINT ,ET=Quadrilateral\n',numnode,numelem);
for i=1:numnode
    fprintf(fid, '%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n', vorodi1(i,1),vorodi1(i,2),AS(i,1),AS(i,2),AS(i,3), chegalinew(i,1));
end
for i=1:numelem
    fprintf(fid, '%6.0f %6.0f %6.0f %6.0f\n',vorodi2(i,1),vorodi2(i,2),vorodi2(i,3),vorodi2(i,4));
end

mod=fopen('chegali_gid.txt', 'wt');
for i=1:numnode
fprintf(mod, '%6.0f %12.8f\n',i,chegalinew(i,1));
end




darsad_khata=(KHATAKOL/TANESHKOL)*100;
















































function [z]=tarsim(vorodi1,vorodi2,tedad)
 RR=vorodi2;
%tarsim sheklha
R=vorodi1(:,1:3);
X1=R(:,2);
Y1=R(:,3);
for i=1:tedad
    XW1=[X1(RR(i,2)) X1(RR(i,3))];
 YW1=[Y1(RR(i,2)) Y1(RR(i,3))];
plot(XW1,YW1);
hold on
end
for i=1:tedad
    XW2=[X1(RR(i,2)) X1(RR(i,5))];
 YW2=[Y1(RR(i,2)) Y1(RR(i,5))];
plot(XW2,YW2)
hold on
end
for i=1:tedad
    XW3=[X1(RR(i,3)) X1(RR(i,4))];
 YW3=[Y1(RR(i,3)) Y1(RR(i,4))];
plot(XW3,YW3)
hold on
end
for i=1:tedad
    XW4=[X1(RR(i,4)) X1(RR(i,5))];
 YW4=[Y1(RR(i,4)) Y1(RR(i,5))];
plot(XW4,YW4)
hold on
end







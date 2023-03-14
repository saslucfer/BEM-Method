clc
clear all
r=xlsread('1-r.xlsx')
C=xlsread('1-chord_distribution')
B=3
R=1.8
landa=6
for i=1:15
landa_r(i)=landa*(r(i)/R)
Beta(i)=asind((C(i)*3.*landa_r(i).*B)./(8*pi*r(i)))
end
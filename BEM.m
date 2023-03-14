clc;
clear all;
close all;
R=1.8;
R_hub=0.4;
A=pi*(R^2);
B=3;
N=15;
Re=81712;
ro=1.2;
landa=6;
e=10e-9;
V0=8;
omega=(landa*V0)/R;
U=[1:1:30];
landa=omega*R./U;
r=xlsread('1-r.xlsx');
C=xlsread('1-chord_distribution.xlsx');
beta=xlsread('1-twist_distribution.xlsx');
S=xlsread('1-E387-R81712.xlsx');
a=zeros(15,300,30);
a_prime=zeros(15,300,30);

%%
for k=1:30
      
   for i=1:15
       sigma(i)=(B*C(i))/(2*pi*r(i));
   
       for j=1:300
              v_wind(i,j,k)=(1-a(i,j,k))*U(k);
              vrot(i,j,k)=r(i)*omega*(1+a_prime(i,j,k));
              vrel(i,j,k)=sqrt(((v_wind(i,j,k)).^2)+((vrot(i,j,k)).^2));
              landa_r(i,j,k)=vrot(i,j,k)/v_wind(i,j,k);
%             landa_r(i,j,k)=(landa(k)*r(i))/R;      %2nd equ
%             phi(i,j,k)=atand((v_wind(i,j,k))/(vrot(i,j,k))); %2nd equ
              phi(i,j,k)=atand((1-a(i,j,k))/((1+a_prime(i,j,k))*landa_r(i,j,k)));
              
              AoA=S(:,1);
              lift=S(:,2);
              drag=S(:,3);
              alpha(i,j,k)=phi(i,j,k)-beta(i);
              alpha1=floor(alpha(i,j,k));
              [rw1,col1]=find(alpha1==AoA);
               CL(i,j,k)=lift(rw1);
               Cd(i,j,k)=drag(rw1);
               Cn(i,j,k)=CL(i,j,k)*sind(phi(i,j,k))+Cd(i,j,k)*cosd(phi(i,j,k));
               Ct(i,j,k)=CL(i,j,k)*sind(phi(i,j,k))-Cd(i,j,k)*cosd(phi(i,j,k));
               
%              f(i,j,k)=((B/2)*(R-r(i))/(r(i)*sind(phi(i,j,k))));
%              F(i,j,k)=((2/pi)*acos(exp(-f(i,j,k))));            %2nd equ
               
               f_tip(i,j,k)=(B/2)*(R-r(i))/(r(i)*sind(phi(i,j,k)));
               F_tip(i,j,k)=(2/pi)*acosd(exp(-f_tip(i,j,k)));
               f_hub(i,j,k)=(B/2)*(r(i)-R_hub)/(r(i)*sind(phi(i,j,k)));
               F_hub(i,j,k)=(2/pi)*acosd(exp(-f_hub(i,j,k)));
               F(i,j,k)=F_hub(i,j,k)*F_tip(i,j,k);
               
               
%                CTr(i,j,k)=(sigma(i)*(1-a(i,j,k))^2*(CL(i,j,k)*cosd(phi(i,j,k))+Cd(i,j,k)*sind(phi(i,j,k))))/(sind(phi(i,j,k))^2);
%                
%                if  CTr<0.96 
%                          a(i,j+1,k)=1/(1+((4*F(i,j,k)*sind(phi(i,j,k)))/(sigma(i)*CL(i,j,k)*cosd(phi(i,j,k)))));
%                else
%                          a(i,j+1,k)=(1/F(i,j,k))*(0.143+sqrt(0.0203-0.6427*(0.889-CTr(i,j,k))));
%                end
%                
%                          a_prime(i,j+1,k)=1/(((4*F(i,j,k)*cosd(phi(i,j,k)))/(sigma(i)*CL(i,j,k)))-1);


             O=((4*F(i,j,k)*(sind(phi(i,j,k))).^2)/(sigma(i)*Cn(i,j,k)));
             
             if a(i,j,k)<0.33
                    a(i,j+1,k)=1/(((4*F(i,j,k)*(sind(phi(i,j,k))).^2)/(sigma(i)*Cn(i,j,k)))+1);
             else
                    a(i,j+1,k)=(1/2)*(2+O*(1-(0.66))-sqrt(((O*(1-(0.66))+2).^2)+(4*(O*((0.33).^2)-1))));
             end
                    a_prime(i,j+1,k)=((1/2)*(sqrt(1+((4*a(i,j+1,k)*(1-a(i,j+1,k)))/(landa_r(i,j,k).^2)))-1));
               
              e1=(a(i,j+1,k)-a(i,j,k));
              e2=(a_prime(i,j+1,k)-a_prime(i,j,k));
              
              if (e1&&e2)<e
                  break
              end
       end
       f_a(i,k)=a(i,j+1,k);
       f_a_prime(i,k)=a_prime(i,j+1,k);
       f_phi(i,k)=phi(i,j,k);
       f_CL(i,k)=CL(i,j,k);
       f_Cd(i,k)=Cd(i,j,k);
       f_Ct(i,k)=Ct(i,j,k);
       f_Cn(i,k)=Cn(i,j,k);
       f_v_wind(i,k)=v_wind(i,j,k);
       f_vrot(i,k)=vrot(i,j,k);
       f_F(i,k)=F(i,j,k);
       f_landa_r(i,k)=landa_r(i,j,k);
       

       
   end
end


%%
Q=zeros(30,1);
P=zeros(30,1);
Cp=zeros(30,1);

for k=1:30
    for i=1:15
    Q(k)=Q(k)+(R/15)*((0.5*1.2*B)*f_v_wind(i,k)*(f_vrot(i,k))*C(i)*(f_Ct(i,k))*r(i))/(sind(f_phi(i,k))*cosd(f_phi(i,k)));
    end
    P(k)=Q(k)*omega;
    Cp(k)=P(k)/(1.2*pi*(R.^2)*(U(k).^3));
end  

%%
% Cp=zeros(30,1)
% for k=1:30
%       for i=1:N
% Cp(k)=Cp(k)+(8*(landa(k)/N)/(landa(k)^2))*f_F(i,k)*(sind(f_phi(i,k))^2)*(cosd(f_phi(i,k))-f_landa_r(i,k)*sind(phi(i,k)))*(sind(phi(i,k))+f_landa_r(i,k)*cosd(phi(i,k)))*(1-(Cd(i,k)/CL(i,k))*cotd(phi(i,k)))*f_landa_r(i,k)
%       end
% end

%%

figure
poly1=polyfit(landa.',Cp,5);
landapol=linspace(0,10,30);
Cp1=polyval(poly1,landapol);
plot(landapol,Cp1,'-ok','Markerfacecolor','B');
hold on;
grid on;

Cp_article=xlsread('2-Cp_exp.xlsx');
landa_article=xlsread('2-landa_exp.xlsx');
poly2=polyfit(landa_article,Cp_article,3);
Cpex2=polyval(poly2,landa_article);
plot(landa_article,Cpex2,'r')

xlabel('Tip speed ratio')
ylabel('power coefficient')
title('power curve at Re = 81,712 for E780')
legend('my BEM','articles BEM')


%%
figure
plot(U,P,'-ok','Markerfacecolor','K')
xlabel('wind velocity (m/s)')
ylabel('power(W)')
title('P vs U curve')
grid on


%%
alphaq=xlsread('2-alpha.xlsx');
CLq=xlsread('2-CL.xlsx');
Cdq=xlsread('2-Cd.xlsx');

figure
plot(alphaq,CLq,'-ok','Markerfacecolor','G')
xlabel('Angle of attack')
ylabel('Lift coefficiet')
grid on

figure
plot(alphaq,Cdq,'-ok','Markerfacecolor','C')
xlabel('Angle of attack')
ylabel('Drag coeffiecient')
grid on

%%
ratio=zeros(77,1);
for n=1:77
    ratio(n)=CLq(n)/Cdq(n);
end
figure
plot(alphaq,ratio,'-ok','Markerfacecolor','Y')
title('CL/Cd vs alpha')
xlabel('Angle of attack')
ylabel('CL/Cd')
grid on



clc;
clear all;
close all;
format long
%%Coupled wave approach
lambda=10^-6;
n1=1.55;
n2=1.5;
epso=8.854*10^-12;
er1=n1^2;
er2=n2^2;
deler=er1-er2;
eo=n2^2;
L1=lambda/(4*n1);
L2=lambda/(4*n2);
L_unit=L1+L2;
half=L_unit/2;
N=100;
L=N*L_unit;
syms x;
if(L1==half)
    kappa=0;
elseif (L1>half)
    er_minus1=(deler/L_unit)*(-1j*L_unit/(2*pi))*(-1+exp(1j*pi));
    kappa=er_minus1*pi/(n2*lambda)
else
    er_minus1=(deler/L_unit)*(-1j*L_unit/(2*pi))*(-1+exp(1j*2*pi*L1/L_unit));
    kappa=er_minus1*pi/(n2*lambda)
end

kx1=2*pi*n1/lambda;
kx2=2*pi*n2/lambda;

A_zero=1;
kx=2*pi*n2/lambda;
delk= (2*kx)-(2*pi/L_unit)
s=sqrt( abs(kappa)^2 - (delk/2)^2  );
s=abs(s);
A= exp(1j.*delk.*x/2).*[s.*cosh(s.*(L-x)) + (1j.*(delk./2).*sinh(s.*(L-x)))]./[s.*cosh(s.*(L)) + (1j.*(delk./2).*sinh(s.*(L)))];
B= exp(-1j.*delk.*x/2) .*(-1j.* conj(kappa).*sinh(s.*(L-x)))./ [s.*cosh(s.*(L)) + ((1j.*delk./2).*sinh(s.*(L)))];

figure('units','normalized','outerposition',[0 0 1 1]);
fplot(abs(A),[0 L]);
xlabel('x(m)'); 
hold on;
fplot(abs(B),[0 L]);
xlim([0 3.279569892473118e-05]);
ylim([-0.05 1.05]);
axis on
legend('A(x)','B(x)');
title(['Field amplitudes using coupled wave approach with n_1 = ' num2str(n1) ' and n_2 = ' num2str(n2)]);
hold off
 set(gcf, 'PaperPositionMode', 'auto');
 saveas(gcf,"coupled.png");

RCA =abs((-1j* conj(kappa)*sinh(s*(L)))/ [s*cosh(s*(L)) + ((1j*delk/2)*sinh(s*(L)))])^2
TCA =abs((exp(1j*delk*L/2)*s)/[s*cosh(s*L) + (1j*(delk/2)*sinh(s*L))])^2

%%Matrix approach
lambda=10^-6;

n2=1.5;
N=200;
L1= lambda/(4*n1);
L2=lambda/(4*n2);
LL=100*(L1+L2);
k1x=2*pi*n1/lambda;
k2x=2*pi*n2/lambda;
kix=2*pi*n2/lambda;
ktx=2*pi*n2/lambda;
A=[1 1; k2x -1*k2x];
A=inv(A);
B=[1 1; k1x -1*k1x];
C=[exp(-1j*k1x*L1) 0; 0 exp(1j*k1x*L1)];
ABC=A*B;
ABC=ABC*C;
D=[1 1; k1x -1*k1x];
D=inv(D);
E=[1 1; k2x -1*k2x];
F=[exp(-1j*k2x*L2) 0; 0 exp(1j*k2x*L2)];
DEF=D*E;
DEF=DEF*F;
X=[1 1; ktx -1*ktx];
X=inv(X);
Y=[1 1; k2x -1*k2x];
Z=[exp(-1j*k2x*L2) 0; 0 exp(1j*k2x*L2)];
XYZ=X*Y;
XYZ=XYZ*Z;
R=[1 1; k1x -1*k1x];
R=inv(R);
S=[1 1; kix -1*kix];
RS=R*S;
T=XYZ;
for i =1:N-1
    if(mod(i,2))
    T=T*DEF;
    else
         T=T*ABC;
    end
end
T=T*RS;
T=real(T);
Tin=inv(T);
Ai=real(Tin(1,1))/real(Tin(1,1));
Bi=real(Tin(2,1))/real(Tin(1,1));
j=1;

 [CC1]=RS*[Ai;Bi];
       a(j)=CC1(1,1);
       b(j)=CC1(2,1);
       j=j+1;
 

for i =2:N

    if(mod(i,2))
       [CC1]=ABC*[a(i-1);b(i-1)];
       a(j)=CC1(1,1);
       b(j)=CC1(2,1);
       j=j+1;
    else
       [CC1]=DEF*[a(i-1);b(i-1)];
       a(j)=CC1(1,1);
       b(j)=CC1(2,1);
       j=j+1;
   end
end
[CC1]=XYZ*[a(i-1);b(i-1)];
       a(j)=CC1(1,1);
       b(j)=CC1(2,1);
       
xa=linspace(0,LL,201);
figure('units','normalized','outerposition',[0 0 1 1]);
plot(xa,abs(a))
hold on
plot(xa,abs(b))
hold off
xlim([0 3.279569892473118e-05]);
ylim([-0.05 1.05]);
xlabel('x(m)');
legend('A(x)','B(x)');
title(['Field amplitudes using matrix approach with n_1 = ' num2str(n1) ' and n_2 = ' num2str(n2)]);
set(gcf, 'PaperPositionMode', 'auto');
 saveas(gcf,"matrix.png");

RMA= abs(Bi/Ai)^2
TMA= abs(1/ Tin(1,1))^2


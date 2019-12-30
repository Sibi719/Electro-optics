clc;
clear all;
close all;
n1=1.55;
n2=1.5;
lambda=10^-6;
L1=lambda/(4*n1);
L2=lambda/(4*n2);
AA=L1+L2;
deler=n1^2-n2^2;
ao=deler*L1/(L1+L2);
N=1000;%%%Change this to vary the number of fourier co-efficients
erx=zeros(1,2000);
x=linspace(0,AA,2000);
for m=-N:N
    if(m==0)
        delerm=ao;
    else
       delerm=(-deler/(1j*m*2*pi))*[-1+exp(-1j*2*m*pi*n2/(n1+n2))];
    end 
  temp= delerm.*exp((1j*m*2*pi.*x)/AA);
  erx=erx+temp;
end
figure
plot(x,real(erx));
title(['\Delta\epsilon_r(x) for m = -' num2str(N) ' to ' abs(num2str(N))]);
xlabel('x(m)');
xlim([0 AA])

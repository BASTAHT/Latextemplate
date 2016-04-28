
clear all
clc
tic
% V=0;
% n=1;
hbar=1.05457e-34;
m=0.3*9.10938e-31;
e=1.602e-19;
L=10e-9;
ns=10^27;
T=10^-4;
% En=pi^2*hbar^2*n^2/(2*m*L^2)

V=0:0.5e-3:1;
for i=1:1:length(V)
nmin(i)=0;
nmax(i)=floor(sqrt((3*pi^2*ns)^(2/3)*L^2/pi^2+m*L^2/ (pi^2*hbar^2)*e*V(i))); %blabalbalbtestest om te kijken of hij afbreekt
nergie=(3*pi^2*ns)^(2/3)*L^2/pi^2;
voltenergie=m*L^2/(pi^2*hbar^2)*e*V(i);
ding(i)=nergie-voltenergie;
if ding(i)>=0
nmin(i)=ceil(sqrt((3*pi^2*ns)^(2/3)*L^2/pi^2-m*L^2/(pi^2*hbar^2)*e*V(i)));
else
    nmin(i)=1;
end
for n=nmin(i):1:nmax(i)
    E(n,i)=pi^2*hbar^2*n^2/(2*m*L^2)-0.5*e*V(i);
end
end

E(E<0)=0;
E(E>0)=1;

for j=1:max(nmax)
    En=pi^2*hbar^2*j^2/(2*m*L^2);
    Er(j,:)=E(j,:)*En;
    gamma=hbar/(2*L)*sqrt(2/m)*T*sqrt(Er);
    I=e/hbar*gamma;
end

% gamma=hbar/(2*L)*sqrt(2/m)*T*sqrt(Er)


c = sum(I);
d=sum(E);
%
% figure(1)
% plot(V,c)
figure(2)
plot(V,d)
xlabel('Voltage (V)')
ylabel('Number of available states')
title('Available states')
Process('filename','Available states')

I10=e/(2*L)*sqrt(2*pi^2*hbar^2*10^2/(2*m*L^2)/m)*T;

L=2*sqrt(1e-12/pi)/2*sqrt(pi);
ak=4*pi^2/L^2;
kf=(3*pi^2*ns)^(1/3);
af=pi*kf^2;
transstates=af/ak;
totalcurrent=2.*c.*transstates;
figure(3)
plot(V,totalcurrent)
xlabel('Voltage (V)')
ylabel('Total current (A)')
title('Total current')
Process('filename','Total current')
toc
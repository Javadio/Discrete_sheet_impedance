clc 
close all
clear

x = [-16.3786741065241	-132.448018643830	-84.6850846659483	-142.977204168630	-164.664661862392	-80.1099147534326	-73.3764427503033	-74.7878847300909	-86.2855489375774];



%% basic parameters
N=10; % 2N+1 harmonics in total
NumEl = 9;

LossTan = 0;

eps_d=2.2*(1-1i*LossTan);
d=1.57*1e-3;

f0=1e9;
freq=8*f0;

eps_0=8.854e-12; % permittivity of freespace
mue_0=4*pi*10^-7; % permeability of free space
c_0=1/sqrt(eps_0*mue_0);
th_d=60;

D_m=c_0/freq/sind(th_d); % periodicity of modulation
beta_m=2*pi/D_m;

z1=1i*x(1);
z2=1i*x(2);
z3=1i*x(3);
z4=1i*x(4);
z5=1i*x(5);
z6=1i*x(6);
z7=1i*x(7);
z8=1i*x(8);
z9=1i*x(9);

M=500;

ys=zeros(1,2*M+1);

W1=D_m/NumEl;
W2=D_m/NumEl;
W3=D_m/NumEl;
W4=D_m/NumEl;
W5=D_m/NumEl;
W6=D_m/NumEl;
W7=D_m/NumEl;
W8=D_m/NumEl;
W9=D_m/NumEl;


N_1=round(W1/(D_m/(2*M+1)));
N_2=round((W1+W2)/(D_m/(2*M+1)));
N_3=round((W1+W2+W3)/(D_m/(2*M+1)));
N_4=round((W1+W2+W3+W4)/(D_m/(2*M+1)));
N_5=round((W1+W2+W3+W4+W5)/(D_m/(2*M+1)));
N_6=round((W1+W2+W3+W4+W5+W6)/(D_m/(2*M+1)));
N_7=round((W1+W2+W3+W4+W5+W6+W7)/(D_m/(2*M+1)));
N_8=round((W1+W2+W3+W4+W5+W6+W7+W8)/(D_m/(2*M+1)));
N_9=round((W1+W2+W3+W4+W5+W6+W7+W8+W9)/(D_m/(2*M+1)));


ys(1:N_1)=1/z1; 
ys(N_1+1:N_2)=1/z2; 
ys(N_2+1:N_3)=1/z3; 
ys(N_3+1:N_4)=1/z4; 
ys(N_4+1:N_5)=1/z5;
ys(N_5+1:N_6)=1/z6;
ys(N_6+1:N_7)=1/z7;
ys(N_7+1:N_8)=1/z8;
ys(N_8+1:N_9)=1/z9;
% ys=adm_array(y1,y2);
yp=fliplr(fftshift(fft(ys)/(2*M+1)));     % notice that fft is defined as exp(1i*m*2*pi/D), but all of our formulas is defined as exp(-1i*m*2*pi/D), therefore we should flip yp!!!!!!!!

Ys=zeros(2*N+1,2*N+1);

for rr=1:2*N+1
    
    for cc=1:2*N+1

    Ys(rr,cc)=yp(rr-cc+M+1)*exp(-1i*(rr-cc)*pi);  
    
    end
end

w_0=2*pi*freq; % angular frequency of incident wave  % 1.0004 have very strong evanecent waves

%% build matrix

%  define wave admittance of free space

k_0=w_0*sqrt(eps_0*mue_0); % free space wavevector
    
theta_i=0;

k_0z=k_0*sind(theta_i); % transverse component of wavevector
% d=2*pi/(k_0*sqrt(eps_d-(cosd(theta_i))^2))/4.00000000001;

Z0=zeros(2*N+1,2*N+1); % define a matrix that is have 2*N+1 rows and 2*N+1 columns

for rr=1:2*N+1
    for cc=1:2*N+1
        k_0zn=k_0z+(cc-N-1)*beta_m; % transverse wavevector of n harmonic at cc column
        k_0n=w_0*sqrt(mue_0*eps_0); % wavevector of n harmonics at cc column
%         k_0xn=sqrt(k_0n^2-k_0zn^2); % normal wavevector of n harmonics at cc column
      
     if abs(k_0zn)<=abs(k_0n) %propatating wave condition
      k_0xn=sqrt(k_0n^2-k_0zn^2); % normal wavevector of n harmonics at cc column
     else
      k_0xn=-sqrt(k_0n^2-k_0zn^2); % normal wavevector of n harmonics at cc column   
     end
      
        if cc==rr
%           Z0(rr,cc)=k_0xn/eps_0/w_n;  % for TM mode
          Z0(rr,cc)=mue_0*w_0/k_0xn;  % for TE mode
        end 
        
    end
    
end 

Y0=pinv(Z0);

% define wave admittance of grounded substrate

Yd=zeros(2*N+1,2*N+1);

for rr=1:2*N+1
    for cc=1:2*N+1

        k_dzn=k_0z+(cc-N-1)*beta_m; % transverse wavevector of n harmonic at cc column
        k_dn=w_0*sqrt(mue_0*eps_0*eps_d); % wavevector of n harmonics at cc column
%         k_dxn=sqrt(k_dn^2-k_dzn^2); % normal wavevector of n harmonics at cc column

     if abs(k_dzn)<=abs(k_dn) %propatating wave condition
      k_dxn=sqrt(k_dn^2-k_dzn^2); % normal wavevector of n harmonics at cc column
     else
      k_dxn=-sqrt(k_dn^2-k_dzn^2); % here, k_dxn minus and plus have no difference because Zd and gamma have the same sign. 
     end
     
        gamma_n=1i*k_dxn;
%         Zd_n=k_dxn/eps_d/eps_0/w_n;  % characteristic impedance substrate for TM mode
        Zd_n=mue_0*w_0/k_dxn;  % characteristic impedance substrate for TE mode
    if cc==rr
      Yd(rr,cc)=1/(Zd_n*tanh(gamma_n*d));  % ??????
    end

    end
end

Yt=Ys+Yd;

% I=eye(2*N+1);

% Gamma=pinv(I+Yt*Z0)*(Yt*Z0-I);  % TM polarization
Gamma=pinv(Yt+Y0)*(Y0-Yt); % for TE polarization
R=abs(Gamma(:,N+1));

figure
stem(-N:N, R)
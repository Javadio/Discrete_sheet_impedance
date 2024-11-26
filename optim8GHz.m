
close all
clear

rng default % For reproducibility
ms=MultiStart('UseParallel',true,'Display','iter');
tic % start time

%% initial point
% this part is to specify initial point. Matlab will serarch for
% solutions near the initial point. Normally, for unsmooth objective function,
% one intial point is not enough. Therefore we should use multiple intial
% points. MultiStart function can provide multiple user-defined initial
% points

% specify the initial points, this determines the accuracy and
% time of optimization. More intial points will take more optimization time
% but will provide more accurate results.

Num_xx = 2; % specify the numbers of each x

Num_x1=Num_xx; % specify the numbers of x1 overal number should be not more than 1000
Num_x2=Num_xx; % specify the numbers of x2
Num_x3=Num_xx; % specify the numbers of x3
Num_x4=Num_xx; % specify the numbers of x4
Num_x5=Num_xx; % specify the numbers of x4
Num_x6=Num_xx; % specify the numbers of x4
Num_x7=Num_xx; % specify the numbers of x4
Num_x8=Num_xx; % specify the numbers of x4
Num_x9=Num_xx; % specify the numbers of x4

x1=linspace(-1000, -30, Num_x1);
x2=linspace(-1000, -30, Num_x2);
x3=linspace(-1000, -30, Num_x3);
x4=linspace(-1000, -30, Num_x4);
x5=linspace(-1000, -30, Num_x5);
x6=linspace(-1000, -30, Num_x6);
x7=linspace(-1000, -30, Num_x7);
x8=linspace(-1000, -30, Num_x8);
x9=linspace(-1000, -30, Num_x9);

% specify the first initial point. this is just a guess but it is necessary
% in the program.
x0=[-400 -519 -200 -480 -500 -600 -350 -200 -200];

% store initial points as a matrix
ptmatrix =(combvec(x1, x2, x3, x4, x5, x6, x7, x8, x9))';
tpoints = CustomStartPointSet(ptmatrix);

%% define objective function

objective_function=@fun;

%% variable bounds
lb = [];
ub = [];

%% linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

%% nonlinear constraints
nonlincon = @nlcons;

%% define problem

problem = createOptimProblem('fmincon','x0',x0,...
    'objective',objective_function,'lb',[],'ub',[],'nonlcon',nonlincon);
% [x,fval,ef,output,lambda] = fmincon(objective_function,x0,A,b,Aeq,beq,lb,ub,nonlincon);

% options = optimoptions('fmincon','FiniteDifferenceStepSize',1e-16);

[x,fval,exitflag,output,solutions] = run(ms,problem,tpoints);

[c,ceq]=nlcons(x);

% show final optimized results
disp(['Final Objective: ' num2str(objective_function(x))])
% print solution
disp('Solution')
disp(['x1 = ' num2str(x(1))])
disp(['x2 = ' num2str(x(2))])
disp(['x3 = ' num2str(x(3))])
disp(['x4 = ' num2str(x(4))])
disp(['x5 = ' num2str(x(5))])
disp(['x6 = ' num2str(x(6))])
disp(['x7 = ' num2str(x(7))])
disp(['x8 = ' num2str(x(8))])
disp(['x9 = ' num2str(x(9))])

%% basic parameters

f0=1e9;
freq=8*f0;
MunEl = 9; % specify the number of elements in the supercell

eps_0=8.854e-12; % permittivity of freespace
mue_0=4*pi*10^-7; % permeability of free space
c_0=1/sqrt(eps_0*mue_0);
th_d=60;

D_m=c_0/freq/sind(th_d); % periodicity of modulation
beta_m=2*pi/D_m;

% Num_D_m=100;

M=500; % Here 2N<M  % numbers of Fourier coefficients in total is 2*M+1, considering -M to M

W1=D_m/MunEl;
W2=D_m/MunEl;
W3=D_m/MunEl;
W4=D_m/MunEl;
W5=D_m/MunEl;
W6=D_m/MunEl;
W7=D_m/MunEl;
W8=D_m/MunEl;
W9=D_m/MunEl;

%%
z1=1i*x(1);
z2=1i*x(2);
z3=1i*x(3);
z4=1i*x(4);
z5=1i*x(5);
z6=1i*x(6);
z7=1i*x(7);
z8=1i*x(8);
z9=1i*x(9);
%%construct ys

N_1=round(W1/(D_m/(2*M+1)));
N_2=round((W1+W2)/(D_m/(2*M+1)));
N_3=round((W1+W2+W3)/(D_m/(2*M+1)));
N_4=round((W1+W2+W3+W4)/(D_m/(2*M+1)));
N_5=round((W1+W2+W3+W4+W5)/(D_m/(2*M+1)));
N_6=round((W1+W2+W3+W4+W5+W6)/(D_m/(2*M+1)));
N_7=round((W1+W2+W3+W4+W5+W6+W7)/(D_m/(2*M+1)));
N_8=round((W1+W2+W3+W4+W5+W6+W7+W8)/(D_m/(2*M+1)));
N_9=round((W1+W2+W3+W4+W5+W6+W7+W8+W9)/(D_m/(2*M+1)));

ys=zeros(1,2*M+1);
ys(1:N_1)=1/z1;
ys(N_1+1:N_2)=1/z2;
ys(N_2+1:N_3)=1/z3;
ys(N_3+1:N_4)=1/z4;
ys(N_4+1:N_5)=1/z5;
ys(N_5+1:N_6)=1/z6;
ys(N_6+1:N_7)=1/z7;
ys(N_7+1:N_8)=1/z8;
ys(N_8+1:N_9)=1/z9;


zs(1:N_1)=z1;
zs(N_1+1:N_2)=z2;
zs(N_2+1:N_3)=z3;
zs(N_3+1:N_4)=z4;
zs(N_4+1:N_5)=z5;
zs(N_5+1:N_6)=z6;
zs(N_6+1:N_7)=z7;
zs(N_7+1:N_8)=z8;
zs(N_8+1:N_9)=z9;
% yp=fliplr(fftshift(fft(ys)/(2*M+1)));     % notice that fft is defined as exp(1i*m*2*pi/D), but all of our formulas is defined as exp(-1i*m*2*pi/D), therefore we should flip yp!!!!!!!!
zp=fliplr(fftshift(fft(zs)/(2*M+1)));     % notice that fft is defined as exp(1i*m*2*pi/D), but all of our formulas is defined as exp(-1i*m*2*pi/D), therefore we should flip yp!!!!!!!!

%% reconstruct the admittance profile
Num_R=200; % reduced numbers of Fourier coefficients in total is 2*Num_R+1, considering -Num_R to Num_R

zp_red=zp(M+1-Num_R:M+1+Num_R); % reduced array

znn=linspace(0, D_m-D_m/(2*M+1), 10001); % make the points enough to ensure accuracy
Adm_sum=zeros(1, 10001);

for kk=1:2*Num_R+1
    
    Adm=zp_red(kk)*exp(-1i*(kk-Num_R-1)*beta_m.*znn);
    Adm_sum=Adm+Adm_sum;
    
end

figure
plot(znn/D_m,real(Adm_sum),znn/D_m,imag(Adm_sum))

toc % end time

function z=fun(x)

LossTan = 0;

eps_d=2.2*(1-1i*LossTan);
d=1.57*1e-3; %thickness

%% basic parameters
f0=1e9;
freq=8*f0;
MunEl = 9;

eps_0=8.854e-12; % permittivity of freespace
mue_0=4*pi*10^-7; % permeability of free space
c=1/sqrt(eps_0*mue_0);
th_d=60;

D_m=c/freq/sind(th_d); % periodicity of modulation
beta_m=2*pi/D_m;

D_m=2*pi/beta_m; % periodicity of modulation

% Num_D_m=100;

M=500; % Here 2N<M  % numbers of Fourier coefficients in total is 2*M+1, considering -M to M

W1=D_m/MunEl;
W2=D_m/MunEl;
W3=D_m/MunEl;
W4=D_m/MunEl;
W5=D_m/MunEl;
W6=D_m/MunEl;
W7=D_m/MunEl;
W8=D_m/MunEl;
W9=D_m/MunEl;
%% basic parameters
N=20; % 2N+1 harmonics in total

z1=1i*x(1);
z2=1i*x(2);
z3=1i*x(3);
z4=1i*x(4);
z5=1i*x(5);
z6=1i*x(6);
z7=1i*x(7);
z8=1i*x(8);
z9=1i*x(9);
ys=zeros(1,2*M+1);

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

% z=abs(Gamma(N+1,N+1))+abs(abs(Gamma(N,N+1))-0.5);
z=abs(R(N+1)-0)+abs(R(N+2)-1/sqrt(cosd(th_d)))+abs(R(N)-0)+0.3*abs(sum(R(N+5:2*N+1))-0);

end

%% define nonlinear constraints

function [c,ceq] = nlcons(x)
c1 = -x(1) -800;
c2 = x(1) -1000;
c3 = -x(2) -800;
c4 = x(2) -1000;
c5 = -x(3) -800;
c6 = x(3) -1000;
c7 = -x(4) -800;
c8 = x(4) -1000;
c9 = -x(5) -800;
c10 = x(5) -1000;
c11 = -x(6) -800;
c12 = x(6) -1000;
c13 = -x(7) -800;
c14 = x(7) -1000;
c15 = -x(8) -800;
c16 = x(8) -1000;
c17 = -x(9) -800;
c18 = x(9) -1000;

c=[ c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18];
ceq = [] ;
end

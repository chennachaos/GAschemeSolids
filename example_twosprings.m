%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script for computing the numerical solutions of two-spring model
%
% Author    : Dr Chennakesava Kadapa
% Date      : 26-Apr-2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
more off;
format long;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters

m1 = 0.0; m2 = 1.0; m3 = 1.0;

k1 = 10^4; k2 = 1.0;

M = [m2 0.0;0.0 m3];
C = [0.0 0.0;0.0 0.0];
K = [k1+k2 -k2; -k2 k2];

u0 = [0.0 0.0]';
v0 = [0.0 0.0]';
a0 = M\(-C*v0 - K*u0);

wp = 1.2;
T = 2*pi/wp;

% time step size
dt = 0.1;

% final time
tf = 20.0;

% constants for the time integration scheme
rhoInf= 0.5; % spectral radius for generalised-alpha schemes
gamma = 0.5; beta  = 0.25;

% gamma and beta are used only for the Newmark-beta scheme
% Otherwise, they are reset
td = timeSteppingParameters_Solid(5, rhoInf, dt, gamma, beta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrays to store the solution

tt = 0.0:dt:tf;

Nt=max(size(tt));

uExact=zeros(Nt,2);
vExact=zeros(Nt,2);
aExact=zeros(Nt,2);

uNum=zeros(Nt,2);
vNum=zeros(Nt,2);
aNum=zeros(Nt,2);

uNum(1,:) = u0;
vNum(1,:) = v0;
aNum(1,:) = a0;

R1=zeros(Nt,1);


dispPrev = u0;
veloPrev = v0;
accePrev = a0;
dispDotPrev = v0;

fsCur = [0.0 0.0]';

% For non-linear problems the following two steps should be done inside the loops
Ktemp = td(5)*M + td(6)*C + td(7)*K; % effective stiffness matrix

invKtemp = inv(Ktemp);


% time loop
for  timeStep=2:Nt
    tnp1 = tt(timeStep); % current time instal
    tn   = tt(timeStep-1); % previous time instant
    tCur = td(2)*tnp1 + (1.0-td(2))*tn;

    % predictor (initial guess) for the Newton-Raphson scheme
    % d_{n+1}
    disp = dispPrev;

    % Newton-Raphson iterations loop
    for iter = 1:10
        % compute v_{n+1}, a_{n+1}, ...
        velo     = td(10)*disp + td(11)*dispPrev + td(12)*veloPrev + td(13)*accePrev + td(14)*dispDotPrev;
        acce     = td(15)*disp + td(16)*dispPrev + td(17)*veloPrev + td(18)*accePrev + td(19)*dispDotPrev;
        dispDot  = td(20)*disp + td(21)*dispPrev + td(22)*veloPrev + td(23)*accePrev + td(24)*dispDotPrev;

        % compute d_{n+af}, v_{n+af}, a_{n+am}, ...
        dispCur = td(2)*disp + (1.0-td(2))*dispPrev;
        veloCur = td(2)*velo + (1.0-td(2))*veloPrev;
        acceCur = td(1)*acce + (1.0-td(1))*accePrev;
        dispDotCur = td(1)*dispDot + (1.0-td(1))*dispDotPrev;

        % external force
        fsCur(1) = k1*sin(wp*tCur);

        resi = fsCur - M*acceCur - C*veloCur - K*dispCur;

        rNorm = norm(resi,2);

        fprintf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);

        if( rNorm < 1.0e-6 ) 
            break;
        end

        disp = disp + invKtemp*resi; % incremental solution is computed and added on the fly
    end
    

    % reaction force
    %R1(timeStep) = 1.0e7*(sin(1.2*tCur)-disp(1));

    uNum(timeStep,:) = disp;
    vNum(timeStep,:) = velo;
    aNum(timeStep,:) = acce;

    % copy variables ()_{n} <-- ()_{n+1}
    dispPrev    = disp;
    dispDotPrev = dispDot;
    veloPrev    = velo;
    accePrev    = acce;
end

plot(tt, uNum(:,1), 'b')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script for the pendulum example by Crisfield
% Nonlinear truss element is used
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

NEL = 1; % number of elements
nen = 2; % number of nodes per element
ndof = 2; % number of DOFs per node



elmDat = zeros(10,1);

elmDat(1) = 1 ;
elmDat(2) = 6.57; % rho0
elmDat(3) = 1.0;  % A0
elmDat(4) = 10^4;    % E


nlocal = nen;
nsize = nlocal*ndof;
NN = NEL + 1;

XX = [0.0 -3.0443;0.0 0.0];

IEN = zeros(NEL, nen);
LM  = zeros(NEL, nsize);

count=1;
for e=1:NEL
    for jj=1:nen
        IEN(e,jj) = count;
        count = count + 1;
    end
    count = count - 1;
end

for e=1:NEL
    count = 1;
    for jj=1:nen
        ind = ndof*(IEN(e,jj)-1)+1;
        for kk=1:ndof
            LM(e,count) = ind;
            ind = ind + 1;
            count = count + 1;
        end
    end
    count = count - 1;
end


% time step size
dt = 0.2;

% final time
tf = 20.0;

% constants for the time integration scheme
rhoInf= 0.5; % spectral radius for generalised-alpha schemes
gamma = 0.5; beta  = 0.25;

% gamma and beta are used only for the Newmark-beta scheme
% Otherwise, they are reset
% See the function for different  schemes available
td = timeSteppingParameters_Solid(4, rhoInf, dt, gamma, beta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrays to store the solution

tt = 0.0:dt:tf;
Nt = max(size(tt));

neq = NN*ndof;
soln=zeros(neq,1);

disp = soln;
velo = soln;
acce = soln;

dispPrev2 = soln;
dispPrev3 = soln;
dispPrev4 = soln;

dispPrev = soln;
veloPrev = soln;
accePrev = soln;
dispDotPrev = soln;

dispCur = soln;
veloCur = soln;
acceCur = soln;


v0 = 7.72; % initial velocity
a0 = 0.0;  % initial acceleration


veloPrev(1) = v0;
accePrev(2) = a0;
dispDotPrev(1) = v0;

assy4r=1:2;

dx = zeros(1,Nt);
dy = zeros(1,Nt);
rot = zeros(1,Nt);

vx = zeros(1,Nt);
vy = zeros(1,Nt);
ax = zeros(1,Nt);
ay = zeros(1,Nt);

bf=[0.0 0.0];

vx(1) = v0;
ay(1) = a0;

% global stiffness matrix and force vector
K_global = zeros(neq,neq);
F_global = zeros(neq,1);

% time loop
for  timeStep=2:Nt
    tnp1 = tt(timeStep); % current time instal
    tn   = tt(timeStep-1); % previous time instant
    tCur = td(2)*tnp1 + (1.0-td(2))*tn;

    % predictor (initial guess) for the Newton-Raphson scheme
    % d_{n+1}
    disp = dispPrev;

    % Newton-Raphson iterations loop
    for iter = 1:20
        % compute v_{n+1}, a_{n+1}, ...
        velo     = td(10)*disp + td(11)*dispPrev + td(12)*veloPrev + td(13)*accePrev + td(14)*dispDotPrev;
        acce     = td(15)*disp + td(16)*dispPrev + td(17)*veloPrev + td(18)*accePrev + td(19)*dispDotPrev;
        dispDot  = td(20)*disp + td(21)*dispPrev + td(22)*veloPrev + td(23)*accePrev + td(24)*dispDotPrev;

        % compute d_{n+af}, v_{n+af}, a_{n+am}, ...
        dispCur = td(2)*disp + (1.0-td(2))*dispPrev;
        veloCur = td(2)*velo + (1.0-td(2))*veloPrev;
        acceCur = td(1)*acce + (1.0-td(1))*accePrev;

        K_global(1:end,1:end) = 0.0;
        F_global(1:end) = 0.0;

        % Loop over elements
        for e = 1:NEL
            [Klocal, Flocal] = Elem_Truss_2D(elmDat, IEN, e, XX, td, disp, dispPrev, veloCur, acceCur, bf);

            K_global = Assembly_Matrix(K_global,Klocal,LM,e,nsize);
            F_global = Assembly_Vector(F_global,Flocal,LM,e,nsize);
        end
        
        %%% Applying Boundary Conditions
        K1 = K_global(assy4r,assy4r);
        F1 = F_global(assy4r);

        rNorm = norm(F1,2);

        printf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
        if(rNorm < 1.0e-8)
            break;
        end
        
        %% solve for the solution and update the displacement
        disp(assy4r) = disp(assy4r) + K1\F1;
    end

    dispPrev4 = dispPrev3;
    dispPrev3 = dispPrev2;
    dispPrev2 = dispPrev;

    dispPrev = disp;
    veloPrev = velo;
    accePrev = acce;
    dispDotPrev = dispDot;

    % data for plotting
    dx(timeStep) = disp(1);
    dy(timeStep) = disp(2);

    vx(timeStep) = velo(1);
    vy(timeStep) = velo(2);

    ax(timeStep) = acce(1);
    ay(timeStep) = acce(2);
end


plot(tt, dy,'k-')


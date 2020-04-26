%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Element subroutine for nonlinear truss element in 2D
%
% Author    : Dr Chennakesava Kadapa
% Date      : 26-Apr-2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Klocal, Flocal] = Elem_Truss_2D(elmDat, IEN, e, XX, td, soln, solnPrev, veloCur, acceCur, bf)

ndof = 2;

af = td(2);
am = td(1);
d1 = td(5);


finite = (elmDat(1) == 1) ; % linear/nonlinear formulation
rho0   = elmDat(2); % density
A0     = elmDat(3); % area
E      = elmDat(4); % Young's modulus


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp     = zeros(4,1);
dispPrev = zeros(4,1);
accC     = zeros(4,1);
velC     = zeros(4,1);
B0       = zeros(4,1);
H        = zeros(4,4);


% rotate nodal displacements and compute nodal positions on element axis   

x0(1) = XX(IEN(e,1),1);
y0(1) = XX(IEN(e,1),2);
x0(2) = XX(IEN(e,2),1);
y0(2) = XX(IEN(e,2),2);

disp(1) = soln(ndof*(IEN(e,1)-1)+1);
disp(2) = soln(ndof*(IEN(e,1)-1)+2);
disp(3) = soln(ndof*(IEN(e,2)-1)+1);
disp(4) = soln(ndof*(IEN(e,2)-1)+2);

dispPrev(1) = solnPrev(ndof*(IEN(e,1)-1)+1);
dispPrev(2) = solnPrev(ndof*(IEN(e,1)-1)+2);
dispPrev(3) = solnPrev(ndof*(IEN(e,2)-1)+1);
dispPrev(4) = solnPrev(ndof*(IEN(e,2)-1)+2);


dispC = af*disp + (1.0-af)*dispPrev;


% compute the orientation of the element
dx = x0(2) - x0(1);
dy = y0(2) - y0(1);
    
h0 = sqrt(dx*dx+dy*dy);

c0x = dx/h0;
c0y = dy/h0;

B0(1) = -c0x;   B0(2) = -c0y;   B0(3) =  c0x;   B0(4) =  c0y;
B0 = B0/h0;

H(1,1) =  1.0;    H(1,3) = -1.0;
H(2,2) =  1.0;    H(2,4) = -1.0;
H(3,1) = -1.0;    H(3,3) =  1.0;
H(4,2) = -1.0;    H(4,4) =  1.0;
H = H/(h0*h0);


if(finite == 0)
   fprintf('small strain algorithm \n');
end


Klocal = zeros(4,4);
Mlocal = zeros(4,4);
Flocal = zeros(4,1);

B = B0;
strain = B0'*dispC ;

if(finite == 1)
   B = B0 + H*dispC;
   strain = strain + (0.5*dispC'*(H*dispC));
end

stress = E * strain;
F = A0*stress;

% residual

Flocal = Flocal - (F*h0)*B;

% material stiffness
    
fact = E*A0*h0;
Klocal = Klocal + ((fact*B)*B') ;

% geometric stiffness
if(finite == 1)
   Klocal = Klocal + ((F*h0)*H);
end

Klocal = Klocal*af;

% inertia
% lumped mass matrix
    
%     fact = rho0*A0*h0/2.0;
% 
%     Mlocal(1,1) = fact;
%     Mlocal(2,2) = fact;
%     Mlocal(3,3) = fact;
%     Mlocal(4,4) = fact;

% consistent mass matrix

fact  = rho0*A0*h0/6.0;
fact2 = 2.0*fact;

Mlocal(1,1) = fact2;
Mlocal(2,2) = fact2;
Mlocal(3,3) = fact2;
Mlocal(4,4) = fact2;
    
Mlocal(1,3) = fact;
Mlocal(2,4) = fact;
Mlocal(3,1) = fact;
Mlocal(4,2) = fact;
 

accC(1) = acceCur(ndof*(IEN(e,1)-1)+1);
accC(2) = acceCur(ndof*(IEN(e,1)-1)+2);
accC(3) = acceCur(ndof*(IEN(e,2)-1)+1);
accC(4) = acceCur(ndof*(IEN(e,2)-1)+2);

Klocal = Klocal + d1*Mlocal;
Flocal = Flocal - Mlocal*accC;


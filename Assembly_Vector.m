%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subroutine for assemblling the element vector into global vector
%
% Author    : Dr Chennakesava Kadapa
% Date      : 26-Apr-2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function global_vector=Assembly_Vector(global_vector,local_vector,LM,e,nen)

for aa=1:nen
    mm=LM(e,aa);
    if mm~=0
        global_vector(mm,1)=global_vector(mm,1)+local_vector(aa,1);
    end
end

endfunction
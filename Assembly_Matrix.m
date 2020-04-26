%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subroutine for assemblling the element matrix into global matrix
%
% Author    : Dr Chennakesava Kadapa
% Date      : 26-Apr-2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function global_matrix=Assembly_Matrix(global_matrix,local_matrix,LM,e,nen)

for aa=1:nen
    mm=LM(e,aa);
    if mm~=0
        for bb=1:nen
            nn=LM(e,bb);
            if nn~=0
                global_matrix(mm,nn)=global_matrix(mm,nn)+local_matrix(aa,bb);
            end
        end
    end
end

endfunction
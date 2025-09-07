function isposdef = check_posdef(M)

test1 = issymmetric(M);

if test1 == 1
    
    
    if all(eig(M)>0) == 1
        disp('Positive definite')
    else
        disp('Not positive definite')
    end
    
else
    disp('Matrix is not symmetric')
end


end

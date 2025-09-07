function constraint = constraintmaker(N, E_constr, e_input)
    
    E_constr2 = kron(eye(N),E_constr);
    e_input2 = repmat(e_input,N,1);
    
    constraint.E = E_constr2;
    constraint.e = e_input2;
end
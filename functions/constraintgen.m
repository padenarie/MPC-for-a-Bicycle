function constraints=constraintgen(predmod,dim)

constraints.Ae=predmod.S(end-dim.nx+1:end,:);
constraints.be=predmod.T(end-dim.nx+1:end,:);

end
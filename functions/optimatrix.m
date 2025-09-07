function [H, F] = optimatrix(predmod,dim,weight)

Q = weight.Q;
R = weight.R;
P = weight.P;
W = predmod.W;

Qbar = kron(Q,eye(dim.N-1));

Qbarbar = [Qbar zeros(size(Qbar,1),size(Q,2));...
           zeros(size(Q,1),size(Qbar,2)) P];
       
Rbar = kron(R,eye(dim.N));

H = predmod.S(dim.nx+1:end,:)'*Qbarbar*predmod.S(dim.nx+1:end,:) + Rbar;
F = (predmod.T(dim.nx+1:end,:))'*Qbarbar*predmod.S(dim.nx+1:end,:);

end
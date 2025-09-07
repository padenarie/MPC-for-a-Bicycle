function [H, F, omega, delta] = optimatrix(sys,N,Q,R,P)
A = sys.A;
B = sys.B;
C = sys.c;

for i = 1:N
    for j = 1:N
        if j-i < 0
            omega(j*size(B,1)-(size(B,1)-1):j*size(B,1),i*(size(B,2))-(size(B,2)-1):i*size(B,2)) = zeros(size(B));
        else
            omega(j*size(B,1)-(size(B,1)-1):j*size(B,1),i*(size(B,2))-(size(B,2)-1):i*size(B,2)) =  A^(j-i)*B;
        end
    end
   delta(i*size(A,1)-(size(A,1)-1):i*size(A,1),:) = A^i;
end

Qbar = kron(Q,eye(N-1));

Qbarbar = [Qbar zeros(size(Qbar,1),size(Q,2));...
           zeros(size(Q,1),size(Qbar,2)) P];
       
Rbar = kron(R,eye(N));

H = omega'*Qbarbar*omega + Rbar;
F = delta'*Qbarbar*omega;

end
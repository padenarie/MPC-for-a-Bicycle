function feasibleset = Feasibleset(LTI,dim,predmod,constraint)

% State constraints
A = constraint.F;
b = constraint.es;

% Input constraints
C = constraint.E;
d = constraint.ei;

% Select the first component, for each temporal step
G_x = predmod.T(1:dim.nx:end,:);
H_x = predmod.S(1:dim.nx:end,:);

Psi_x = b*ones(dim.N+1,1);

% Define input constraints
H_u = kron(eye(dim.N),C);
Psi_u = repmat(d,dim.N,1);

G = [G_x; zeros(size(H_u,1),dim.nx)];
H = [H_x; H_u];
Psi = -[Psi_x; Psi_u];

%%
% Algorithm implementation

G_i = [G H(:,1:end-1)];
H_i = H(:,end);
Psi_i = Psi;

for i = dim.N-1:-1:0

    [P_i, gamma_i] = single_input(G_i,H_i,Psi_i);
    
    G_i = P_i(:,1:end-1);
    H_i = P_i(:,end);
    Psi_i = gamma_i;
    
end

predmod.T = P_i;
gamma = gamma_i;

% Eliminate repeated rows
norep = unique([predmod.T gamma], 'rows');

predmod.T = norep(:,1:end-1);
gamma = norep(:,end);

% Plot region
x_1 = -3:.01:3;
x_2 = -3:.01:3;

[X_1, X_2] = meshgrid(x_1, x_2);

for i = 1:size(predmod,1)
   
    map(:,:,i) = (predmod(i,2).*X_2 <= -predmod(i,1).*X_1 - gamma(i));
    
end

figure,
contourf(X_1, X_2, min(double(map), [], 3), [1 1]);
xlabel('$x_1$','Interpreter','latex'), ylabel('$x_2$','Interpreter','latex'), grid on;
xlim([-2 2]), ylim([-2 2]);


end























% Fixed point solver
% super simple

% -------- Problem set up ---------
Re = 10;
R = 50;
M = 150; N = 150;

xiMax = log(R)/pi;
dXi = xiMax/(N-1.5);
dEta = 2/M;
eta = [-1:dEta:1-dEta];
% xi = [0:dXi:xiMax+dXi/2];
xi = [xiMax+dXi/2:-dXi:0];

% Form operator and RHS
psOp = formOps(M,N,R);
rhs = formRHS(M,N,R,Re);

% Solve
psiOm = psOp\rhs;
psi = reshape(psiOm(1:N*M), [N, M]);
omega = reshape(psiOm(1+N*M:2*M*N), [N, M]);

% Plot
figure(1);
subplot(1,2,1);
contourf(eta, xi, psi, 50);
subplot(1,2,2);
contourf(eta, xi, omega, 50);

% Hans' plotter here
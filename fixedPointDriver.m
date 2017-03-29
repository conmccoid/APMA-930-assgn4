% Fixed point solver
% super simple

% -------- Problem set up ---------
Re = 1;
R = 10;
M = 200; N = 200;
OmegaPsi = zeros(2*M*N,1);

xiMax = log(R)/pi;
dXi = xiMax/(N-1.5);
dEta = 2/M;
eta = -1:dEta:1-dEta;
% xi = [0:dXi:xiMax+dXi/2];
xi = xiMax+dXi/2:-dXi:0;

% Form operator and RHS
psOp = formOps(M,N,R);
[LL,UU,PP,QQ,RR] = lu(psOp);
% rhs = formRHS(OmegaPsi, M,N,R,Re);

% Solve
IterMax = 1000;
Iter = 0;
test = 1;
tol = 1e-9;
% for k = 1 : IterMax
while Iter<IterMax && test>tol && isnan(test)==0

    Iter = Iter + 1;

    rhs = formRHS(OmegaPsi, M, N, Re, xi, eta, dXi, dEta);
%     rhs = formRHS_2(OmegaPsi, M, N, Re, R); % takes longer
    c = PP * ( RR \ rhs );
    OmegaPsi_new = QQ * ( UU \ ( LL \ c ) );
    
    test = norm(OmegaPsi_new - OmegaPsi)/norm(OmegaPsi_new);
    OmegaPsi = OmegaPsi_new;
    
    disp(['Iteration: ',num2str(Iter),' Test: ',num2str(test)])
end
psiOm = OmegaPsi;
psi = reshape(psiOm(1:N*M), [N, M]);
omega = reshape(psiOm(1+N*M:2*M*N), [N, M]);

%%

% Free stream
psiFree = exp(pi*xi)'*sin(pi*eta);

% Plot
figure(1);
subplot(1,2,1);
contourf(eta, xi, psi+psiFree, 50);
subplot(1,2,2);
contourf(eta, xi, omega, 50);
% Hans' plotter here

rtPsi = WAvg_EXtransform(psi+psiFree,R);
rtOmega = WAvg_EXtransform(omega,R);
figure(2);
subplot(1,2,1)
contourf(rtPsi,100);

subplot(1,2,2)
contourf(rtOmega,100);

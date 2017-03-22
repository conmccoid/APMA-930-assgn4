function rhs = formRHS(Omega, Psi, M,N,R,Re)
xiMax = log(R)/pi;
dXi = xiMax/(N-1.5);
dEta = 2/M;
unkOrd = reshape([1:M*N]', N, M);
rhs = zeros(2*M*N,1);
% xi = 0,xiMax 
for i=1:M
    % Dirichlet on bottom (psi)
    eta = (i-1)*dEta-1;
    rhs(unkOrd(N,i)) = -sin(pi*eta);
end

E = xiMax-dXi/2:dXi:0; E = exp(pi*E);
S = 0:M-1; C = cos(-pi + pi*dEta*S);
S = sin(-pi + pi*dEta*S);

for i=2:M-1         % Eta variable
    for j=2:N-1     % Xi variable
        % interior omega
        idx = unkOrd(j,i)+M*N;
        
        Etaplus = unkOrd(j,i+1);
        Etaminus = unkOrd(j,i-1);
        Xiplus = unkOrd(j+1,i);
        Ximinus = unkOrd(j-1,i);
        
        dOmEta = (Omega(Etaplus+M*N) - Omega(Etaminus+M*N))/(2*dEta);
        dOmXi = (Omega(Xiplus+M*N) - Omega(Ximinus+M*N))/(2*dXi);
        dPsiEta = (Psi(Etaplus) - Psi(Etaminus))/(2*dEta);
        dPsiXi = (Psi(Xiplus) - Psi(Ximinus))/(2*dXi);
        
        rhs(idx) = Re/2*(pi*E(j)*(C(i)*dOmXi - S(i)*dOmEta) - dPsiXi*dOmEta + dPsiEta*dOmXi);
    end
end
for j=2:N-1         % Xi variable
    % omega periodic walls
    idx = unkOrd(j,M)+M*N;
    
    Etaplus = unkOrd(j,1);
    Etaminus = unkOrd(j,M-1);
    Xiplus = unkOrd(j+1,M);
    Ximinus = unkOrd(j-1,M);
    
        dOmEta = (Omega(Etaplus+M*N) - Omega(Etaminus+M*N))/(2*dEta);
        dOmXi = (Omega(Xiplus+M*N) - Omega(Ximinus+M*N))/(2*dXi);
        dPsiEta = (Psi(Etaplus) - Psi(Etaminus))/(2*dEta);
        dPsiXi = (Psi(Xiplus) - Psi(Ximinus))/(2*dXi);
        
        rhs(idx) = Re/2*(pi*E(j)*(C(M)*dOmXi - S(M)*dOmEta) - dPsiXi*dOmEta + dPsiEta*dOmXi);
        
    idx = unkOrd(j,1)+M*N;
    
    Etaplus = unkOrd(j,2);
    Etaminus = unkOrd(j,M);
    Xiplus = unkOrd(j+1,1);
    Ximinus = unkOrd(j-1,1);
    
        dOmEta = (Omega(Etaplus+M*N) - Omega(Etaminus+M*N))/(2*dEta);
        dOmXi = (Omega(Xiplus+M*N) - Omega(Ximinus+M*N))/(2*dXi);
        dPsiEta = (Psi(Etaplus) - Psi(Etaminus))/(2*dEta);
        dPsiXi = (Psi(Xiplus) - Psi(Ximinus))/(2*dXi);
        
        rhs(idx) = Re/2*(pi*E(j)*(C(1)*dOmXi - S(1)*dOmEta) - dPsiXi*dOmEta + dPsiEta*dOmXi);

end
% xi = 0,xiMax 
for i=1:M
    % Dirichlet on bottom (omega)
    eta = (i-1)*dEta-1;
    rhs(unkOrd(N,i)+M*N) = 2*sin(pi*eta)*(dXi*pi-1)/(dXi^2);
end
end

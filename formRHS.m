function rhs = formRHS(OmegaPsi, M, N, Re, xi, eta, dXi, dEta)
% xiMax = log(R)/pi;
% dXi = xiMax/(N-1.5);
% dEta = 2/M;
unkOrd = reshape([1:M*N]', N, M);
rhs = zeros(2*M*N,1);
% xi = 0,xiMax 

% E = xiMax+dXi/2:-dXi:0; 
E = exp(pi*xi);
% S = 0:M-1; 
C = cos(pi*eta);
S = sin(pi*eta);

for i=2:M-1         % Eta variable
    for j=2:N-1     % Xi variable
        % interior omega
        idx = unkOrd(j,i)+M*N;
        
        Etaplus = unkOrd(j,i+1);
        Etaminus = unkOrd(j,i-1);
        Xiplus = unkOrd(j+1,i);
        Ximinus = unkOrd(j-1,i);
%         Ximinus = unkOrd(j+1,i);
%         Xiplus = unkOrd(j-1,i);
        
        dOmEta = (OmegaPsi(Etaplus+M*N) - OmegaPsi(Etaminus+M*N))/(2*dEta);
        dOmXi = (OmegaPsi(Xiplus+M*N) - OmegaPsi(Ximinus+M*N))/(2*dXi);
        dPsiEta = (OmegaPsi(Etaplus) - OmegaPsi(Etaminus))/(2*dEta);
        dPsiXi = (OmegaPsi(Xiplus) - OmegaPsi(Ximinus))/(2*dXi);
        
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
%     Ximinus = unkOrd(j+1,M);
%     Xiplus = unkOrd(j-1,M);
    
        dOmEta = (OmegaPsi(Etaplus+M*N) - OmegaPsi(Etaminus+M*N))/(2*dEta);
        dOmXi = (OmegaPsi(Xiplus+M*N) - OmegaPsi(Ximinus+M*N))/(2*dXi);
        dPsiEta = (OmegaPsi(Etaplus) - OmegaPsi(Etaminus))/(2*dEta);
        dPsiXi = (OmegaPsi(Xiplus) - OmegaPsi(Ximinus))/(2*dXi);
        
        rhs(idx) = Re/2*(pi*E(j)*(C(M)*dOmXi - S(M)*dOmEta) - dPsiXi*dOmEta + dPsiEta*dOmXi);
        
    idx = unkOrd(j,1)+M*N;
    
    Etaplus = unkOrd(j,2);
    Etaminus = unkOrd(j,M);
    Xiplus = unkOrd(j+1,1);
    Ximinus = unkOrd(j-1,1);
%     Ximinus = unkOrd(j+1,1);
%     Xiplus = unkOrd(j-1,1);
    
        dOmEta = (OmegaPsi(Etaplus+M*N) - OmegaPsi(Etaminus+M*N));
        dOmXi = (OmegaPsi(Xiplus+M*N) - OmegaPsi(Ximinus+M*N));
        dPsiEta = (OmegaPsi(Etaplus) - OmegaPsi(Etaminus));
        dPsiXi = (OmegaPsi(Xiplus) - OmegaPsi(Ximinus));
        
        rhs(idx) = Re/2*(pi*E(j)*(C(1)*dOmXi/(2*dXi) - S(1)*dOmEta/(2*dEta)) - ( dPsiXi*dOmEta - dPsiEta*dOmXi )/(4*dXi*dEta));

end
% xi = 0,xiMax 
for i=1:M
    % Dirichlet on bottom
    rhs(unkOrd(N,i)) = -S(i);
    rhs(unkOrd(N,i)+M*N) = 2*pi*S(i)/dXi;
%     rhs(unkOrd(N,i)+M*N) = (- 2/(dXi^2) - pi^2)*S(i);
end
end

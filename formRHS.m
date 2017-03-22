function rhs = formRHS(M,N,R,Re)
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

for i=2:M-1
    for j=2:N-1
        % interior omega
        idx = unkOrd(j,i)+M*N;
        eta = (i-1)*dEta-1;
        xi = (N-j)*dXi;
        r = exp(pi*xi);
        theta = pi*xi;
        % Connor's function here, see formOps for method of referencing
        % other indices/neighboring points
    end
end
for j=2:N-1
    % omega periodic walls
    idx = unkOrd(j,i)+M*N;
    eta = (i-1)*dEta-1;
    xi = (N-j)*dXi;
    r = exp(pi*xi);
    theta = pi*xi;
%     Connor's function here, see formOps for method of referencing
%     other indices/neighboring points
end
% xi = 0,xiMax 
for i=1:M
    % Dirichlet on bottom (omega)
    eta = (i-1)*dEta-1;
    rhs(unkOrd(N,i)+M*N) = 2*sin(pi*eta)*(dXi*pi-1)/(dXi^2);
end
end
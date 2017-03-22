function psOp = formOps(M,N,R)
xiMax = log(R)/pi;
dXi = xiMax/(N-1.5);
dEta = 2/M;
unkOrd = reshape([1:M*N]', N, M);

% -------- Form the laplacian --------
ctr = 1;
% Interior
for i=2:M-1
    for j=2:N-1
        % Self
        iIdx(ctr,1) = unkOrd(j,i);
        jIdx(ctr,1) = unkOrd(j,i);
        val(ctr,1) = -2/(dEta^2)-2/(dXi^2);
        ctr = ctr + 1;
        % Psi omega connection
        iIdx(ctr,1) = unkOrd(j,i);
        jIdx(ctr,1) = unkOrd(j,i)+M*N;
        xi = (N-j)*dXi;
        val(ctr,1) = -pi^2*exp(2*pi*xi);
        ctr = ctr + 1;
        % West
        iIdx(ctr,1) = unkOrd(j,i);
        jIdx(ctr,1) = unkOrd(j,i-1);
        val(ctr,1) = 1/(dEta^2);
        ctr = ctr + 1;
        % East
        iIdx(ctr,1) = unkOrd(j,i);
        jIdx(ctr,1) = unkOrd(j,i+1);
        val(ctr,1) = 1/(dEta^2);
        ctr = ctr + 1;
        % South
        iIdx(ctr,1) = unkOrd(j,i);
        jIdx(ctr,1) = unkOrd(j+1,i);
        val(ctr,1) = 1/(dXi^2);
        ctr = ctr + 1;
        % North
        iIdx(ctr,1) = unkOrd(j,i);
        jIdx(ctr,1) = unkOrd(j-1,i);
        val(ctr,1) = 1/(dXi^2);
        ctr = ctr + 1;
    end
end
% eta = -1,1 periodic
for j=2:N-1
    % Self
    iIdx(ctr,1) = unkOrd(j,1);
    jIdx(ctr,1) = unkOrd(j,1);
    val(ctr,1) = -2/(dEta^2)-2/(dXi^2);
    ctr = ctr + 1;
    % Self
    iIdx(ctr,1) = unkOrd(j,M);
    jIdx(ctr,1) = unkOrd(j,M);
    val(ctr,1) = -2/(dEta^2)-2/(dXi^2);
    ctr = ctr + 1;
    % West
    iIdx(ctr,1) = unkOrd(j,1);
    jIdx(ctr,1) = unkOrd(j,M);
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % West
    iIdx(ctr,1) = unkOrd(j,M);
    jIdx(ctr,1) = unkOrd(j,M-1);
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % East
    iIdx(ctr,1) = unkOrd(j,1);
    jIdx(ctr,1) = unkOrd(j,2);
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % East
    iIdx(ctr,1) = unkOrd(j,M);
    jIdx(ctr,1) = unkOrd(j,1);
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % South
    iIdx(ctr,1) = unkOrd(j,1);
    jIdx(ctr,1) = unkOrd(j+1,1);
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
    % South
    iIdx(ctr,1) = unkOrd(j,M);
    jIdx(ctr,1) = unkOrd(j+1,M);
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
    % North
    iIdx(ctr,1) = unkOrd(j,1);
    jIdx(ctr,1) = unkOrd(j-1,1);
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
    % North
    iIdx(ctr,1) = unkOrd(j,M);
    jIdx(ctr,1) = unkOrd(j-1,M);
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
end
% xi = 0,xiMax 
for i=1:M
    % Homogeneous Neumann on top
    iIdx(ctr,1) = unkOrd(1,i);
    jIdx(ctr,1) = unkOrd(2,i);
    val(ctr,1) = 1;
    ctr = ctr + 1;
    % Dirichlet on bottom
    iIdx(ctr,1) = unkOrd(N,i);
    jIdx(ctr,1) = unkOrd(N,i);
    val(ctr,1) = 1;
    ctr = ctr + 1;
end

for i=2:M-1
    for j=2:N-1
        % Self
        iIdx(ctr,1) = unkOrd(j,i)+M*N;
        jIdx(ctr,1) = unkOrd(j,i)+M*N;
        val(ctr,1) = -2/(dEta^2)-2/(dXi^2);
        ctr = ctr + 1;
        % West
        iIdx(ctr,1) = unkOrd(j,i)+M*N;
        jIdx(ctr,1) = unkOrd(j,i-1)+M*N;
        val(ctr,1) = 1/(dEta^2);
        ctr = ctr + 1;
        % East
        iIdx(ctr,1) = unkOrd(j,i)+M*N;
        jIdx(ctr,1) = unkOrd(j,i+1)+M*N;
        val(ctr,1) = 1/(dEta^2);
        ctr = ctr + 1;
        % South
        iIdx(ctr,1) = unkOrd(j,i)+M*N;
        jIdx(ctr,1) = unkOrd(j+1,i)+M*N;
        val(ctr,1) = 1/(dXi^2);
        ctr = ctr + 1;
        % North
        iIdx(ctr,1) = unkOrd(j,i)+M*N;
        jIdx(ctr,1) = unkOrd(j-1,i)+M*N;
        val(ctr,1) = 1/(dXi^2);
        ctr = ctr + 1;
    end
end
% eta = -1,1 periodic
for j=2:N-1
    % Self
    iIdx(ctr,1) = unkOrd(j,1)+M*N;
    jIdx(ctr,1) = unkOrd(j,1)+M*N;
    val(ctr,1) = -2/(dEta^2)-2/(dXi^2);
    ctr = ctr + 1;
    % Self
    iIdx(ctr,1) = unkOrd(j,M)+M*N;
    jIdx(ctr,1) = unkOrd(j,M)+M*N;
    val(ctr,1) = -2/(dEta^2)-2/(dXi^2);
    ctr = ctr + 1;
    % West
    iIdx(ctr,1) = unkOrd(j,1)+M*N;
    jIdx(ctr,1) = unkOrd(j,M)+M*N;
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % West
    iIdx(ctr,1) = unkOrd(j,M)+M*N;
    jIdx(ctr,1) = unkOrd(j,M-1)+M*N;
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % East
    iIdx(ctr,1) = unkOrd(j,1)+M*N;
    jIdx(ctr,1) = unkOrd(j,2)+M*N;
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % East
    iIdx(ctr,1) = unkOrd(j,M)+M*N;
    jIdx(ctr,1) = unkOrd(j,1)+M*N;
    val(ctr,1) = 1/(dEta^2);
    ctr = ctr + 1;
    % South
    iIdx(ctr,1) = unkOrd(j,1)+M*N;
    jIdx(ctr,1) = unkOrd(j+1,1)+M*N;
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
    % South
    iIdx(ctr,1) = unkOrd(j,M)+M*N;
    jIdx(ctr,1) = unkOrd(j+1,M)+M*N;
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
    % North
    iIdx(ctr,1) = unkOrd(j,1)+M*N;
    jIdx(ctr,1) = unkOrd(j-1,1)+M*N;
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
    % North
    iIdx(ctr,1) = unkOrd(j,M)+M*N;
    jIdx(ctr,1) = unkOrd(j-1,M)+M*N;
    val(ctr,1) = 1/(dXi^2);
    ctr = ctr + 1;
end
% xi = 0,xiMax 
for i=1:M
    % Homogeneous Neumann on top
    iIdx(ctr,1) = unkOrd(1,i)+M*N;
    jIdx(ctr,1) = unkOrd(2,i)+M*N;
    val(ctr,1) = 1;
    ctr = ctr + 1;
    % Dirichlet on bottom
    iIdx(ctr,1) = unkOrd(N,i)+M*N;
    jIdx(ctr,1) = unkOrd(N,i)+M*N;
    val(ctr,1) = 1;
    ctr = ctr + 1;
    iIdx(ctr,1) = unkOrd(N,i)+M*N;
    jIdx(ctr,1) = unkOrd(N-1,i);
    val(ctr,1) = -2/(dXi^2);
    ctr = ctr + 1;
end
psOp = sparse(iIdx,jIdx,val);
end
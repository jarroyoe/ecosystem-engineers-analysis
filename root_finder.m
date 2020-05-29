syms N1 E1 N2 E2
%define parameter values
gamma=0.2;
delta=0.02;
etaR=0.01;
epsilon=0.11;
eta=etaR*epsilon;
muR=1;
mu=eta*muR;

%use w as a scaling factor
w=[(0.001:0.005:2)];

n=length(w);
uniqueSolutions=[];
numSolutions=zeros(n,1);

parfor i=1:n
    %solve equations of N1 and N2
    allSol=vpasolve(0==(-1+E1-(gamma+delta*E1)*N1)*N1+w(i)*mu*(N2-N1),...
        0==(-1+E2-(gamma+delta*E2)*N2)*N2+w(i)*mu*(N1-N2),...
        0==epsilon*(N1-E1)+w(i)*eta*N2,...
        0==epsilon*(N2-E2)+w(i)*eta*N1);

    %make matrix of solutions with numerical values
    sol=[repmat(w(i)*mu,length(allSol.N1),1),double(allSol.N1),...
        double(allSol.E1),double(allSol.N2),...
        double(allSol.E2)];

    %remove complex solutions
    sol(~imag(sol(:,2))==0 | ~imag(sol(:,3))==0 ...
        |~imag(sol(:,4))==0 | ~imag(sol(:,5))==0 ,:)=[];
   
    %remove negative solutions
    sol(sol(:,2)<0 | sol(:,3)<0 | sol(:,4)<0 | sol(:,5)<0,:)=[];
    
    %add new solutions and number of solutions
    uniqueSolutions=vertcat(uniqueSolutions,sol);
    numSolutions(i)=length(sol);
end

syms N1 E1 N2 E2
%calculate eigenvalues for all solutions and the dimension of the stable manifold
vars=[N1 E1 N2 E2];
m=length(uniqueSolutions);
isStable=zeros(m,1);
stableManifoldDimension=zeros(m,1);
parfor i=1:m
    dN1=(-1+E1-(gamma+delta*E1)*N1)*N1+uniqueSolutions(i,1)*mu*(N2-N1);
    dE1=epsilon*(N1-E1)+uniqueSolutions(i,1)*eta*N2;
    dN2=(-1+E2-(gamma+delta*E2)*N2)*N2+uniqueSolutions(i,1)*mu*(N1-N2);
    dE2=epsilon*(N2-E2)+uniqueSolutions(i,1)*eta*N1;
    
    J=[double(subs(gradient(dN1,vars),vars,uniqueSolutions(i,2:5))),...
        double(subs(gradient(dE1,vars),vars,uniqueSolutions(i,2:5))),...
        double(subs(gradient(dN2,vars),vars,uniqueSolutions(i,2:5))),...
        double(subs(gradient(dE2,vars),vars,uniqueSolutions(i,2:5)))];
    realPartEigenvals=real(eig(J));
    isStable(i)=max(realPartEigenvals)<0;
    stableManifoldDimension(i)=length(realPartEigenvals(realPartEigenvals<0));
end

%rescale solutions
uniqueSolutions(:,1)=uniqueSolutions(:,1)*mu;

%save solutions
uniqueSolutions=[uniqueSolutions, isStable, stableManifoldDimension];
save("uniqueSolutions.mat");
writematrix(uniqueSolutions,'uniqueSolutions.csv');

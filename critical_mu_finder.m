syms N1 E1 N2 E2
%define parameter values
gamma=0.2;
delta=0.02;
%etaR=0.01;
%epsilon=(0.01:0.1:2);
epsilon=0.06;
eta=0.1;
muR=1;
mu=eta*muR;

%define resolution of search for critical mu
kmax=50;
n=length((0.0001:0.001:1));
m=length(epsilon);
saddle_bifurcation=zeros(m,1);

%for each value of epsilon, find the smallest mu that reduces the number of solutions
for j=1:m
    k=0;
    numSolutions=9;
    ep=epsilon(j);
    while k<kmax && min(numSolutions)==9
        kn=k*n;
        k=k+1;
        w=(0.0001:0.001:k);
        parfor i=1:n
            %solve equations of N1 and N2
            allSol=vpasolve([0==(-1+E1-(gamma+delta*E1)*N1)*N1+w(kn+i)*mu*(N2-N1);...
                 0==(-1+E2-(gamma+delta*E2)*N2)*N2+w(kn+i)*mu*(N1-N2);...
                 0==ep*(N1-E1)+w(kn+i)*eta*N2;...
                 0==ep*(N2-E2)+w(kn+i)*eta*N1]);
            
            %make matrix of solutions with numerical values
            sol=[double(allSol.N1),...
                double(allSol.E1),double(allSol.N2),...
                double(allSol.E2)];
            
            %remove complex solutions
            sol(~imag(sol(:,1))==0 | ~imag(sol(:,2))==0 ...
                |~imag(sol(:,3))==0 | ~imag(sol(:,4))==0 ,:)=[];
            
            %remove negative solutions
            sol(sol(:,1)<0 | sol(:,2)<0 | sol(:,3)<0 | sol(:,4)<0,:)=[];
            
            %add new number of solutions
            numSolutions(kn+i)=length(sol);
            allSol=[];
            sol=[];
        end
    end
    numSolutions(1)=[];
    if k<kmax
        saddle_bifurcation(j)=mu*w(find(numSolutions<9,1,'first'));
    end
    w=[];
end

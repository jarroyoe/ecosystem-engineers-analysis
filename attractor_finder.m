%define parameters
gamma=0.2;
delta=0.02;
epsilon=(0.01:0.05:1);
mu=(1e-4:2e-5:2e-3);
muR=1;
eta=mu*muR;

n=length(epsilon);
m=length(mu);

%define search space for critical density
kmax=10000;
Tmax=100;
min_invasion_size=zeros(n,m);

%solve the ODE for increasing initial condition until the second patch has a density higher than 10
for i=1:n
    for j=1:m
        if j==1
           k=0;
            ep=epsilon(i);
            u=mu(m-j+1);
            et=eta(m-j+1);
            ODEs=@(t,y) [(-1+y(2)-(gamma+delta*y(2))*y(1))*y(1)+u*(y(3)-y(1)),...
                ep*(y(1)-y(2))+et*y(3),...
                (-1+y(4)-(gamma+delta*y(4))*y(3))*y(3)+u*(y(1)-y(3)),...
                ep*(y(3)-y(4))+et*y(1)]';
            thresholdFlag=false;
            solFlag=[];
            while k<kmax && ~thresholdFlag
                kn=k*100;
                k=k+1;
                N0=(1:1:k*100);
                sol=zeros(1,100);
                parfor l=1:100
                    odeSol=ode45(ODEs,[0 Tmax],[N0(kn+l),38.7,0,0]);
                    timeSteps=length(odeSol.y);
                    sol(l)=odeSol.y(3,timeSteps);
                end
                solFlag=[solFlag,sol];
                thresholdFlag=max(solFlag)>10;
            end
            if(thresholdFlag)
                min_invasion_size(i,m-j+1)=N0(find(solFlag>10,1,'first'));
            else
                min_invasion_size(i,m-j+1)=Inf;
                disp(Inf)
            end 
        elseif min_invasion_size(i,m-j+2)== Inf
            min_invasion_size(i,m-j+1)=Inf;
        else
            k=0;
            ep=epsilon(i);
            u=mu(m-j+1);
            et=eta(m-j+1);
            ODEs=@(t,y) [(-1+y(2)-(gamma+delta*y(2))*y(1))*y(1)+u*(y(3)-y(1)),...
                ep*(y(1)-y(2))+et*y(3),...
                (-1+y(4)-(gamma+delta*y(4))*y(3))*y(3)+u*(y(1)-y(3)),...
                ep*(y(3)-y(4))+et*y(1)]';
            thresholdFlag=false;
            solFlag=[];
            while k<kmax && ~thresholdFlag
                kn=k*100;
                k=k+1;
                N0=(1:1:k*100);
                sol=zeros(1,100);
                parfor l=1:100
                    odeSol=ode45(ODEs,[0 Tmax],[N0(kn+l),38.7,0,0]);
                    timeSteps=length(odeSol.y);
                    sol(l)=odeSol.y(3,timeSteps);
                end
                solFlag=[solFlag,sol];
                thresholdFlag=max(solFlag)>10;
            end
            if(thresholdFlag)
                min_invasion_size(i,m-j+1)=N0(find(solFlag>10,1,'first'));
            else
                min_invasion_size(i,m-j+1)=Inf;
                disp(Inf)
            end
        end
    end
end

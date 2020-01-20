function BinarySwitch(alpha,beta,R,M,r)
if nargin<5
    r=1;
end

S=@(X) r*(1-alpha)*(1-X).^6;

for i=2:6
    if (i-1)<=M
        S=@(X) S(X) +r* X.^(i-1).*(1-X).^(7-i).*...
            (nchoosek(5,i-1)-alpha*nchoosek(6,i-1));
    else
        S=@(X) S(X) + R*X.^(i-1).*(1-X).^(7-i).*...
            (nchoosek(5,i-1)-beta*nchoosek(6,i-1));
    end
    
end

S=@(X) S(X)-R*beta*X.^6;
X=linspace(0,1,1000);

figure(1011)
plot(X,X.*S(X),'b',[0 1],[0 0],'k','LineWidth',2)
box on
pause(0.1)

run = input('Run Discrete Simulation? Yes=1, No=0  ');

if run
    C0=input('Specify initial population density C0 in the interval [0,1]  ');
    
    if C0<0 || C0>1
        error('ERROR: C0 must be in the interval [0,1]')
    end
    
    RP=R*ones(1,7);
    RP(1:M+1)=r;
    RD=R*beta*ones(1,7);
    RD(1:M+1)=r*alpha;
    RM=ones(1,7)*max(RP)*100;
    
    if r>0
        TEND=10*max(r,R);
    else
        TEND=10;
    end
    time=linspace(0,TEND,1e3);
    
    [QQ,X]=gillBinary(RM,RP,RD,C0,TEND);
    
    options=odeset('RelTol',1e-8,'AbsTol',1e-8);
    sol = ode45(@(t,x) x.*S(x), time,C0,options);
    
    
    
    C=deval(sol,time)';
    
    figure(100)
    hold off
    plot(X,QQ,'r--','LineWidth',2);
    hold on
    plot(time,C,'k','LineWidth',2);
    axis([0 TEND 0 1])
end

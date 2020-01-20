function Err=BinarySwitchData(file,M,r,R,alpha,beta,C0)
% Inputs: file: an xls dataset (first column is time points, second column 
% is rescaled population density datapoints)
% M,r,R,alpha,beta: Binary Switch Model parameters
% C0: Binary Switch Model initial population density
%
% Output: Err: combined least-squares error for all time points.

D =importdata(strcat(file,'.xls'));

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

figure(4)
plot(X,X.*S(X),'b',[0 1],[0 0],'k','LineWidth',2)
box on
pause(0.1)

TEND=D(end,1);
time=linspace(0,TEND,1e3);
    
    options=odeset('RelTol',1e-8,'AbsTol',1e-8);
    sol = ode45(@(t,x) x.*S(x), time,C0,options);
    Cdata=deval(sol,D(:,1))';
    
    
    C=deval(sol,time)';


figure(5)
hold off
plot(D(:,1),D(:,2),'ko','MarkerSize',5)
hold on
plot(time,C,'r','LineWidth',2)


Err = norm(abs(Cdata-D(:,2))).^2;



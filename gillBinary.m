function [QQ,X]=gillBinary(RM,RP,RD,C0,TEND)
%modified from 2A, but uses individual rates rather than r_I <M and G
%otherwise.

Y=10;
m=C0; % initial seeding
n=100; %lattice nodes in each direction
n2=115;

Tend=TEND;
NN=length(RM);



%% Construct nearest neighbour index structure (S)
XX=zeros(n,n2);
YY=zeros(n,n2);
for i=1:n
    for j=1:n2
        if mod(j,2)==0
            XX(i,j)=i;
            YY(i,j)=j*sqrt(3)/2;
        else
            XX(i,j)=i+1/2;
            YY(i,j)=j*sqrt(3)/2;
        end
    end
end


r=(-1+sqrt(1+4*(NN-1)/3))/2;
S=(n*n2+1)*ones(n*n2,NN-1);

for i=1:n
    for j=1:n2
        Z=(XX-XX(i,j)).^2+(YY-YY(i,j)).^2-r^2;
        V=find(Z<=1e-1 & Z>-r^2);
        S(i+(j-1)*n,1:length(V))=V;

        
    end
end
S(S==0)=n*n2+1;





%% Gillespie
parfor i=1:Y
    C=double(rand(n,n2)<m);
    Q0=sum(sum(C));
    while Q0==0
        C=double(rand(n,n2)<m);
        Q0=sum(sum(C));
    end
    Q=Q0;
    C(1,n2+1)=0; %this is the pointer for all 'null values' of boundary sites
    
    
    T=0;
    j=1;
    tau=0;
    
    
    
   
    Qend=Q0;
    
    Jx=randi([1,n],1,10000);
    Jy=randi([1,n2],1,10000);
    W=rand(1,10000);
    W2=rand(1,10000);
    W3=rand(1,10000);
    y=1;
    while T<Tend && Qend<n*n2 && Qend>0
        if C(1,n2+1)~=0
            error('FATAL ERROR: Null pointer has been altered.')
        end
        %find a random occupied site
        while C(Jx(y),Jy(y))==0
            y=y+1;
            if y==10001
                Jx=randi([1,n],1,10000);
                Jy=randi([1,n2],1,10000);
                W=rand(1,10000);
                W2=rand(1,10000);
                W3=rand(1,10000);
                y=1;
            end
        end
        JX=Jx(y);
        JY=Jy(y);
        JJ=JX+(JY-1)*n;
        Near=sum(C(S(JJ,:)));
        
        P1=RM(Near+1);
        P2=RP(Near+1);
        P3=RD(Near+1);
        tau(j+1)=tau(j)+log(1/W3(y))/((P1+P2+P3)*Q(j));
        R=(P1+P2+P3)*Q(j)*W2(y);
        Q(j+1)=Q(j);
        
        if R<P1*Q(j)
            % cell movement
            I=ceil(W(y)*6);
            
            if I==1
                Ind = [JX+1 JY];
            elseif I==2
                Ind = [JX-1 JY];
            elseif I==3
                Ind = [JX JY+1];
            elseif I==4
                Ind = [JX JY-1];
            elseif I==5
                if mod(JX,2)==0
                    Ind = [JX-1 JY+1];
                else
                    Ind = [JX-1 JY-1];
                end
            else
                if mod(JX,2)==0
                    Ind = [JX+1 JY+1];
                else
                    Ind = [JX+1 JY-1];
                end
            end
            
            if Ind(1)<=n && Ind(2)<=n2 && Ind(1)>0 && Ind(2)>0 && ...
                    C(Ind(1),Ind(2))==0
                C(Ind(1),Ind(2))=1;
                C(JX,JY)=0;
            else
                %movement does not happen
            end
            
            
        elseif R<(P1+P2)*Q(j)
            % cell proliferation
            
            I=ceil(W(y)*6);
            if I==1
                Ind = [JX+1 JY];
            elseif I==2
                Ind = [JX-1 JY];
            elseif I==3
                Ind = [JX JY+1];
            elseif I==4
                Ind = [JX JY-1];
            elseif I==5
                if mod(JX,2)==0
                    Ind = [JX-1 JY+1];
                else
                    Ind = [JX-1 JY-1];
                end
            else
                if mod(JX,2)==0
                    Ind = [JX+1 JY+1];
                else
                    Ind = [JX+1 JY-1];
                end
            end
         
            if Ind(1)<=n && Ind(2)<=n2 && Ind(1)>0 && Ind(2)>0 && ...
                    C(Ind(1),Ind(2))==0
                C(Ind(1),Ind(2))=1;
                Q(j+1)=Q(j+1)+1;
            else
                %prolif does not happen
            end
            
        else
            %cell death
            C(JX,JY)=0;
            Q(j+1)=Q(j+1)-1;
        end
        
        
        T=tau(j+1);
        Qend=Q(j+1);
        j=j+1;
         y=y+1;

       
        
        if y==10001
            Jx=randi([1,n],1,10000);
            Jy=randi([1,n2],1,10000);
            W=rand(1,10000);
            W2=rand(1,10000);
            W3=rand(1,10000);
            y=1;
        end
    end
    X=linspace(0,Tend,1000);
    [tau,U]=unique(tau,'first');
    Q=Q(U);
    G2=griddedInterpolant(tau,Q/(n*n2),'nearest');
    Qi(i,:)=G2(X);
end
X=linspace(0,Tend,1000);

QQ=mean(Qi,1);


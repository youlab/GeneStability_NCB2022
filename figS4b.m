clear;
clc;
global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu;
PoolNumSpecies=100;
PoolNumPlasmid=20;
NumPlasmid=PoolNumPlasmid;
timespan=0:500;
Master=1;
Repeat=40;
plasmid=0*ones(Repeat,PoolNumPlasmid);
species=0*ones(Repeat,PoolNumSpecies);
for i=1:Master
gammaPool=0.1*rand(PoolNumSpecies,PoolNumSpecies)-0.05*ones(PoolNumSpecies,PoolNumSpecies);    
sigmaPool=0.2*rand(PoolNumSpecies,1);

kappamax=0.05;
kappaPool=kappamax*rand(PoolNumSpecies,PoolNumPlasmid);
lambdaPool=0.05-0.1*rand(PoolNumSpecies,PoolNumPlasmid);
muPool=0.2+0.4*rand(PoolNumSpecies,1);
D=0.2;

etaPool=0*ones(PoolNumPlasmid,PoolNumSpecies,PoolNumSpecies);
HostOrNot=1*ones(PoolNumSpecies,PoolNumPlasmid);
    for i1=1:PoolNumPlasmid
        etamax=0.01*2*rand;% Here the mean value of HGT rate is defined
        for i2=1:PoolNumSpecies
            for i3=1:PoolNumSpecies
                etaPool(i1,i2,i3)=etamax;
            end
        end
    end

    for j=1:Repeat
        j+(i-1)*Repeat
        for fg=1:1000
            NumSpecies=poissrnd(0.5*PoolNumSpecies,1);
            if NumSpecies>=1&&NumSpecies<=PoolNumSpecies
                break;
            end
        end
        rr=randperm(PoolNumSpecies);
        MapToPool=rr(1:NumSpecies);
        
        kappa=kappaPool(MapToPool,:);
        
        lambda=lambdaPool(MapToPool,:);
        
        mu=muPool(MapToPool,1);
        gamma=gammaPool(MapToPool, MapToPool);    
        for k=1:NumSpecies
            gamma(k,k)=-mu(k)/rand;
        end
        sigma=sigmaPool(MapToPool);
        eta=etaPool(:,MapToPool, MapToPool);
         initial=0;
        for k=1:NumSpecies
            initial(k)=0.01;
                for kk=1:NumPlasmid
                        initial(NumSpecies+(k-1)*NumPlasmid+kk)=0.1*initial(k)*HostOrNot(MapToPool(k),kk);
                end
        end
        warning('');
        opts = odeset('RelTol',1e-3,'NonNegative',1);
        [t,y]=ode45(@multi_plasmid,timespan,initial);
        [warnMsg, warnId] = lastwarn;
        if isempty(warnMsg)
            for k=1:NumPlasmid
                plasmid(j,k)=sum(y(end,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))./sum(y(end,1:NumSpecies));
            end
            for k=1:NumSpecies
                species(j,MapToPool(k))=y(end,k)./sum(y(end,1:NumSpecies));
            end
        end
    end
end
subplot(1,2,1);
bar(1:Repeat,species,1,'stack','edgecolor','none');
colormap lines;
axis([0.5 Repeat+0.5 0 1]);
set(gca,'fontsize',16);
xlabel('community index','fontsize',20);
ylabel('species abundance','fontsize',20);

subplot(1,2,2);
bar(1:Repeat,plasmid,1,'stack','edgecolor','none');
colormap lines;
axis([0.5 Repeat+0.5 0 max(sum(plasmid,2))*1.1]);
set(gca,'fontsize',16);
xlabel('community index','fontsize',20);
ylabel('plasmid abundance','fontsize',20);

set(gcf,'position',[100 100 400 300]);

saveas(gcf,'fig4b.fig');
saveas(gcf,'fig4b.png');

function dydt=multi_plasmid(t,y)
        global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu;
        dydt(NumSpecies*(1+NumPlasmid),1)=0;
        for i=1:NumSpecies
            sum=0;
            for j=1:NumPlasmid
            sum=sum+lambda(i,j)*y(NumSpecies+NumPlasmid*(i-1)+j);
            end
            Neg=0;
            Pos=0;
            for j=1:NumSpecies
                if gamma(j,i)<0
                Neg=Neg-gamma(j,i)*y(j);
                end
                if gamma(j,i)>0
                Pos=Pos+gamma(j,i)*y(j);
                end
            end
            
            mui=mu(i)-Neg+sigma(i)*Pos/(Pos+1);
            dydt(i,1)=y(i)/(sum+y(i))*mui*y(i)-D*y(i);

            for j=1:NumPlasmid
            betaij=(1+lambda(i,j))/(1+lambda(i,j)+(sum-lambda(i,j)*y(NumSpecies+NumPlasmid*(i-1)+j))/y(i));
            muij=mui/(1+lambda(i,j));
                summ=0;
                for k=1:NumSpecies
                    summ=summ+eta(j,k,i)*y(NumSpecies+(k-1)*NumPlasmid+j);
                end
                dydt(NumSpecies+(i-1)*NumPlasmid+j,1)=betaij*muij*y(NumSpecies+(i-1)*NumPlasmid+j)+(y(i)-y(NumSpecies+(i-1)*NumPlasmid+j))*summ-(kappa(i,j)+D)*y(NumSpecies+(i-1)*NumPlasmid+j);
            end
        end      
end

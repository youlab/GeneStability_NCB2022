clear;
clc;
global NumSpecies NumPlasmid eta kappa D lambda delta gamma mu;
PoolNumSpecies=100; %define the number of species in each species pool;
PoolNumPlasmid=20; %define the number of plasmids in each species pool;
NumPlasmid=PoolNumPlasmid;% define the number of plasmids transferred in each local community
timespan=0:200;% define the timespan of numerical simulations
Master=5;%define the number of species pools in this simulation
Repeat=40;%define the number of local communities created for each species pool

for i=1:Master %for each species pool
gammaPool=0.1*rand(PoolNumSpecies,PoolNumSpecies)-0.05*ones(PoolNumSpecies,PoolNumSpecies); %randomize the gamam values, which describe the species interactions  
deltaPool=0.2*rand(PoolNumSpecies,1);% randomize delta, which describes the maximum growth benefit provided by positive species interactions

kappamax=0.05;
kappaPool=kappamax*rand(PoolNumSpecies,PoolNumPlasmid);% randomize kappa, which describe the plasmid loss rate
lambdaPool=0.05-0.1*rand(PoolNumSpecies,PoolNumPlasmid);%radomize lambda, which describe the plasmid fitness effects
muPool=0.2+0.4*rand(PoolNumSpecies,1);% randomize mu, the maximum growth rate of each species
D=0.2*rand;%randomize dilution rate

etaPool=0*ones(PoolNumPlasmid,PoolNumSpecies,PoolNumSpecies);

% Next, we randomize the host range of each plasmid
HostOrNot=0*ones(PoolNumSpecies,PoolNumPlasmid);
for i1=1:PoolNumPlasmid
    temp1=min(fix(rand*PoolNumSpecies)+1,PoolNumSpecies);
    temp2=randperm(PoolNumSpecies);
    HostOrNot(temp2(1:temp1),i1)=1;
end

%Based on the plasmid host range, we randomize the plasmid transfer rates
    for i1=1:PoolNumPlasmid
        etamax=0.2*rand;
        for i2=1:PoolNumSpecies
            for i3=1:PoolNumSpecies
                etaPool(i1,i2,i3)=etamax*rand*HostOrNot(i2,i1)*HostOrNot(i3,i1);
            end
        end
    end

    for j=1:Repeat %assembly local communities for each species pool.
        j+(i-1)*Repeat
        for fg=1:1000
            NumSpecies=poissrnd(0.5*PoolNumSpecies,1);% randomize the number of species in each local population
            if NumSpecies>=1&&NumSpecies<=PoolNumSpecies
                break;
            end
        end
        rr=randperm(PoolNumSpecies);%randomly sample species from the pool to local communities
        MapToPool=rr(1:NumSpecies);
        
        %Next, determine all the kinetic parameters in the local
        %communities
        kappa=kappaPool(MapToPool,:);
        lambda=lambdaPool(MapToPool,:);
        mu=muPool(MapToPool,1);
        gamma=gammaPool(MapToPool, MapToPool);    
        for k=1:NumSpecies
            gamma(k,k)=-mu(k)/rand;
        end
        delta=deltaPool(MapToPool);
        eta=etaPool(:,MapToPool, MapToPool);
        
        %Next,determine the initial densities of all the species in the
        %local communities
         initial=0;
        for k=1:NumSpecies
            initial(k)=0.01;
                for kk=1:NumPlasmid
                        initial(NumSpecies+(k-1)*NumPlasmid+kk)=0.1*initial(k)*HostOrNot(MapToPool(k),kk);%Also randomize the initial densities of plasmid-carrying cells
                end
        end
        
        [t,y]=ode45(@multi_plasmid,timespan,initial);%run the ODE simulations

        for k=1:NumPlasmid
            plasmid(i,j,k)=sum(y(end,NumSpecies+k:NumPlasmid:NumSpecies+(NumSpecies-1)*NumPlasmid+k))./sum(y(end,1:NumSpecies));%caculate the relative abundance of plasmid-carrying cells
        end
        
    end
end

threshold=10.^(-2.5);% set the threshold for plasmid being 'lost'; if plasmid abundance is smaller than this threhold, the plasmid is lost in this community and plasmid abundance is treated as 0 
plasmid_cleared=plasmid.*(plasmid>threshold);

fp=mean(plasmid_cleared,2);%calculate the mean abundance of each plasmid
stability=mean(plasmid_cleared,2)./std(plasmid_cleared,0,2);%calculate the functional stability
plot(fp(:),stability(:),'k.','markersize',20);
hold on;

set(gca,'fontsize',16);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('f_p','fontsize',20);
ylabel('functional stability','fontsize',20);
axis([threshold 1 10^(-0.5) 10^2]);
set(gcf,'position',[100 100 350 350]);

function dydt=multi_plasmid(t,y)% the ODEs of the complex communities
        global NumSpecies NumPlasmid eta kappa D lambda delta gamma mu;
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
            
            mui=mu(i)-Neg+delta(i)*Pos/(Pos+1);
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

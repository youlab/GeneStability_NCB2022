clear;
clc;
close all;
global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu;
load('MyColormap.mat');
PoolNumSpecies=100;
PoolNumPlasmid=20;
NumPlasmid=PoolNumPlasmid;
timespan=0:500;
Master=200;
Repeat=40;
plasmid=0*ones(Master,Repeat,PoolNumSpecies,PoolNumPlasmid);
transfer=0*ones(Master,PoolNumPlasmid);
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
        etamax=0.1*rand;
        for i2=1:PoolNumSpecies
            for i3=1:PoolNumSpecies
                etaPool(i1,i2,i3)=etamax*rand;
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

        [t,y]=ode45(@multi_plasmid,timespan,initial);

        for jkl=1:NumSpecies
            for k=1:NumPlasmid
                plasmid(i,j,jkl,k)=y(end,NumSpecies+(jkl-1)*NumPlasmid+k)./sum(y(end,1:NumSpecies));
                temp=etaPool(k,:,:);
                transfer(i,k)=mean(temp(:));
            end
        end

    end
end

subplot(1,4,1);
PlasmidExpress=0*ones(Master,Repeat,PoolNumPlasmid);
ExpressPool=1+0*rand(Master,PoolNumSpecies,PoolNumPlasmid);
for i=1:Master
    for j=1:Repeat
        for k=1:PoolNumPlasmid
            for jkl=1:PoolNumSpecies
                PlasmidExpress(i,j,k)=PlasmidExpress(i,j,k)+plasmid(i,j,jkl,k)*ExpressPool(i,jkl,k);
            end
        end
    end
end
C=linspecer(5);
x=mean(PlasmidExpress,2);
yy=mean(PlasmidExpress,2)./std(PlasmidExpress,0,2);

X=transfer(:);
Y=yy(:);
XXX=[X,Y];
xmi=0;
xma=0.05;
ymi=0;
yma=100;
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{xmi:0.0005:xma ymi:0.5:yma});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(xmi,xma,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(ymi,yma,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
caxis([0 6]);
c=mymap;
colormap(c);
set(gca,'fontsize',16);
xlabel('transfer rate','fontsize',20);
ylabel('functional stability','fontsize',20);
set(gcf,'position',[100 100 400 400])
hold on;
plot([0 0],[ymi yma],'k-');hold on;
plot([xmi xma],[0 0],'k-');hold on;
plot([xmi xma],[yma yma],'k-');hold on;
plot([xma xma],[ymi yma],'k-');hold on;

X_input=X;
Y_input=Y;
rep=size(X_input);
group_meds=xmi:0.005:xma;
group_wid=0.0025;
med=0;
err=0;
for i=1:length(group_meds)
    pin=1;
    zan=0;
    group_med=group_meds(i);
    for j=1:rep
        if X_input(j)>=group_med-group_wid&&X_input(j)<group_med+group_wid
            zan(pin)=Y_input(j);
            pin=pin+1;
        end
    end
    med(i)=mean(zan);
    err(i)=std(zan);
end
h_err=errorbar(group_meds,med,err,'o-','MarkerSize',6);
h_err.CapSize=1;
h_err.Color='k';%[0.5 0.5 0.5];
h_err.LineWidth=2;
hold on;

subplot(1,4,2);
PlasmidExpress=0*ones(Master,Repeat,PoolNumPlasmid);
ExpressPool=0.8+0.2*rand(Master,PoolNumSpecies,PoolNumPlasmid);
for i=1:Master
    for j=1:Repeat
        for k=1:PoolNumPlasmid
            for jkl=1:PoolNumSpecies
                PlasmidExpress(i,j,k)=PlasmidExpress(i,j,k)+plasmid(i,j,jkl,k)*ExpressPool(i,jkl,k);
            end
        end
    end
end
C=linspecer(5);
x=mean(PlasmidExpress,2);
yy=mean(PlasmidExpress,2)./std(PlasmidExpress,0,2);

X=transfer(:);
Y=yy(:);
XXX=[X,Y];
xmi=0;
xma=0.05;
ymi=0;
yma=100;
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{xmi:0.0005:xma ymi:0.5:yma});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(xmi,xma,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(ymi,yma,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
caxis([0 6]);
c=mymap;
colormap(c);
set(gca,'fontsize',16);
xlabel('transfer rate','fontsize',20);
ylabel('functional stability','fontsize',20);
set(gcf,'position',[100 100 400 400])
hold on;
plot([0 0],[ymi yma],'k-');hold on;
plot([xmi xma],[0 0],'k-');hold on;
plot([xmi xma],[yma yma],'k-');hold on;
plot([xma xma],[ymi yma],'k-');hold on;

X_input=X;
Y_input=Y;
rep=size(X_input);
group_meds=xmi:0.005:xma;
group_wid=0.0025;
med=0;
err=0;
for i=1:length(group_meds)
    pin=1;
    zan=0;
    group_med=group_meds(i);
    for j=1:rep
        if X_input(j)>=group_med-group_wid&&X_input(j)<group_med+group_wid
            zan(pin)=Y_input(j);
            pin=pin+1;
        end
    end
    med(i)=mean(zan);
    err(i)=std(zan);
end
h_err=errorbar(group_meds,med,err,'o-','MarkerSize',6);
h_err.CapSize=1;
h_err.Color='k';%[0.5 0.5 0.5];
h_err.LineWidth=2;
hold on;

subplot(1,4,3);
PlasmidExpress=0*ones(Master,Repeat,PoolNumPlasmid);
ExpressPool=0.5+0.5*rand(Master,PoolNumSpecies,PoolNumPlasmid);
for i=1:Master
    for j=1:Repeat
        for k=1:PoolNumPlasmid
            for jkl=1:PoolNumSpecies
                PlasmidExpress(i,j,k)=PlasmidExpress(i,j,k)+plasmid(i,j,jkl,k)*ExpressPool(i,jkl,k);
            end
        end
    end
end
C=linspecer(5);
x=mean(PlasmidExpress,2);
yy=mean(PlasmidExpress,2)./std(PlasmidExpress,0,2);

X=transfer(:);
Y=yy(:);
XXX=[X,Y];
xmi=0;
xma=0.05;
ymi=0;
yma=50;
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{xmi:0.0005:xma ymi:0.25:yma});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(xmi,xma,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(ymi,yma,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
caxis([0 6]);
c=mymap;
colormap(c);
set(gca,'fontsize',16);
xlabel('transfer rate','fontsize',20);
ylabel('functional stability','fontsize',20);
set(gcf,'position',[100 100 400 400])
hold on;
plot([0 0],[ymi yma],'k-');hold on;
plot([xmi xma],[0 0],'k-');hold on;
plot([xmi xma],[yma yma],'k-');hold on;
plot([xma xma],[ymi yma],'k-');hold on;

X_input=X;
Y_input=Y;
rep=size(X_input);
group_meds=xmi:0.005:xma;
group_wid=0.0025;
med=0;
err=0;
for i=1:length(group_meds)
    pin=1;
    zan=0;
    group_med=group_meds(i);
    for j=1:rep
        if X_input(j)>=group_med-group_wid&&X_input(j)<group_med+group_wid
            zan(pin)=Y_input(j);
            pin=pin+1;
        end
    end
    med(i)=mean(zan);
    err(i)=std(zan);
end
h_err=errorbar(group_meds,med,err,'o-','MarkerSize',6);
h_err.CapSize=1;
h_err.Color='k';%[0.5 0.5 0.5];
h_err.LineWidth=2;
hold on;

subplot(1,4,4);
PlasmidExpress=0*ones(Master,Repeat,PoolNumPlasmid);
ExpressPool=1*rand(Master,PoolNumSpecies,PoolNumPlasmid);
for i=1:Master
    for j=1:Repeat
        for k=1:PoolNumPlasmid
            for jkl=1:PoolNumSpecies
                PlasmidExpress(i,j,k)=PlasmidExpress(i,j,k)+plasmid(i,j,jkl,k)*ExpressPool(i,jkl,k);
            end
        end
    end
end
C=linspecer(5);
x=mean(PlasmidExpress,2);
yy=mean(PlasmidExpress,2)./std(PlasmidExpress,0,2);

X=transfer(:);
Y=yy(:);
XXX=[X,Y];
xmi=0;
xma=0.05;
ymi=0;
yma=20;
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{xmi:0.0005:xma ymi:0.1:yma});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(xmi,xma,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(ymi,yma,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
caxis([0 6]);
c=mymap;
colormap(c);
set(gca,'fontsize',16);
xlabel('transfer rate','fontsize',20);
ylabel('functional stability','fontsize',20);
set(gcf,'position',[100 100 400 400])
hold on;
plot([0 0],[ymi yma],'k-');hold on;
plot([xmi xma],[0 0],'k-');hold on;
plot([xmi xma],[yma yma],'k-');hold on;
plot([xma xma],[ymi yma],'k-');hold on;

X_input=X;
Y_input=Y;
rep=size(X_input);
group_meds=xmi:0.005:xma;
group_wid=0.0025;
med=0;
err=0;
for i=1:length(group_meds)
    pin=1;
    zan=0;
    group_med=group_meds(i);
    for j=1:rep
        if X_input(j)>=group_med-group_wid&&X_input(j)<group_med+group_wid
            zan(pin)=Y_input(j);
            pin=pin+1;
        end
    end
    med(i)=mean(zan);
    err(i)=std(zan);
end
h_err=errorbar(group_meds,med,err,'o-','MarkerSize',6);
h_err.CapSize=1;
h_err.Color='k';%[0.5 0.5 0.5];
h_err.LineWidth=2;
hold on;

set(gcf,'position',[100 100 1300 300]);


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



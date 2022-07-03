clear;
clc;
close all;
global mu eta D lambda kappa;
CO=linspecer(3);
D=0;
kappa=[0.05 0.05];
muTop=[0.77 0.77 0.69 0.8 0.88];
muMG1655=[0.96 0.73 0.57 0.04 0.08];
eta0=[6.2 6;4.2 13.6]/30;
lambda=[0.07 0.06];
inhis=[1 0.5545 0.377];
dilu=[10^4 10^5];

for i=1:length(dilu)
    for j=1:length(inhis)        
        for k=1
            eta=eta0*inhis(j);
            mu=[muMG1655(k) muTop(k)];
            initial=[0.1 0 0.05 0];
            Top10=0;
            R388=0;
            for kk=1:15
            timespan=0:24;
            [t,y]=ode45(@TwoStrains,timespan,initial);
            Top10(kk)=y(end,2)./(y(end,1)+y(end,2));
            R388(kk)=(y(end,3)+y(end,4))./(y(end,1)+y(end,2));
            initial=y(end,:)/dilu(i);
            end
            
            figure(1);
            subplot(1,2,i);
            plot(0:15,[0.5 R388],'-','linewidth',3,'color',CO(j,:));hold on;
            set(gca,'YScale','log');
            axis([0 15 10^(-4) 10^0.5]);
            set(gca,'fontsize',16);
            yticks(10.^[-4 -3 -2 -1 0]);
            yticklabels('');
        end
        
    end
end
figure(1);
subplot(1,2,1);
xlabel('time (days)','fontsize',24);
ylabel('plasmid abundance','fontsize',24);
set(gcf,'position',[100 100 450 220]);
saveas(gcf,'SerialDilution_MG1655.fig');
saveas(gcf,'SerialDilution_MG1655.png');

function dydt=TwoStrains(t,y)
global mu eta D lambda kappa;
s1=y(1);
s2=y(2);
p1=y(3);
p2=y(4);
st=y(1)+y(2);
dydt=[mu(1)*s1/(s1+lambda(1)*p1)*s1*(1-st)-D*s1;
    mu(2)*s2*(1-st)-D*s2;
    mu(1)/(1+lambda(1))*p1*(1-st)+(s1-p1)*(eta(1,1)*p1+eta(2,1)*p2)-(kappa(1)+D)*p1;
    mu(2)/(1+lambda(2))*p2*(1-st)+(s2-p2)*(eta(1,2)*p1+eta(2,2)*p2)-(kappa(2)+D)*p2];
end

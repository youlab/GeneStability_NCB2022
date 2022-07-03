clear;
clc;
close all;
global mu eta D lambda kappa;
CO=linspecer(5);
D=0;
kappa=[0.05 0.05];
% muTop=[0.8 0.8 0.8 0.8 0.8];
muTop=[0.77 0.77 0.69 0.77 0.81];
muMG1655=[0.96 0.73 0.57 0.04 0.08];
% muMG1655=[1.24 0.89 0.65 0.24 0.00];

eta0=[6.2 6;4.2 13.6]/30;
lambda=[0.07 0.06];
inhis=[1 0.5545 0.377];
dilu=[10^4 10^5];

for i=1:length(dilu)
    for j=1:length(inhis)        
        for k=1:5
            eta=eta0*inhis(j);
            mu=[muMG1655(k) muTop(k)];
            initial=[0.1 0.1 0.05 0.05];
            Top10=0;
            R388=0;
            for kk=1:15
            timespan=0:24;
            [t,y]=ode45(@TwoStrains,timespan,initial);
            Top10(kk)=y(end,2)./(y(end,1)+y(end,2));
            R388(kk)=(y(end,3)+y(end,4))./(y(end,1)+y(end,2));
            initial=y(end,:)/dilu(i);
            end
            R388_5(k)=R388(5);
            R388_10(k)=R388(10);
            R388_15(k)=R388(15);
            figure(1);
            subplot(3,2,2*(j-1)+i);
            plot(0:15,[0.5 R388],'.-','linewidth',3,'color',CO(k,:));hold on;
            set(gca,'YScale','log');
            axis([0 15 10^(-4) 10^0.5]);
            set(gca,'fontsize',16);
            yticks(10.^[-4 -3 -2 -1 0]);
            yticklabels('');
            figure(2);
            subplot(3,2,2*(j-1)+i);
            plot(0:15,[0.5 Top10],'.-','linewidth',3,'color',CO(k,:));hold on;
            axis([0 15 0 1]);
            set(gca,'fontsize',16);
        end
        stab5(i,j)=mean([R388_5])/std([R388_5]);
        stab10(i,j)=mean([R388_10])/std([R388_10]);
        stab15(i,j)=mean([R388_15])/std([R388_15]);
    end
end
figure(1);
subplot(3,2,5);
xlabel('time (days)','fontsize',24);
ylabel('plasmid abundance','fontsize',24);
set(gcf,'position',[100 100 450 650]);
saveas(gcf,'SerialDilutionPlasmid.fig');
saveas(gcf,'SerialDilutionPlasmid.png');
figure(2);
subplot(3,2,5);
xlabel('time (days)','fontsize',24);
ylabel('Top10 abundance','fontsize',24);
set(gcf,'position',[100 100 450 650]);
saveas(gcf,'SerialDilutionStrain.fig');
saveas(gcf,'SerialDilutionStrain.png');

figure(3);
plot([0 5 8],fliplr(stab5(1,:)),'.-','markersize',40,'color',CO(1,:),'linewidth',3);hold on;

plot([0 5 8],fliplr(stab5(2,:)),'o-','markersize',10,'color',CO(1,:),'linewidth',3);hold on;

plot([0 5 8],fliplr(stab10(1,:)),'.--','markersize',40,'color',CO(1,:),'linewidth',3);hold on;

plot([0 5 8],fliplr(stab10(2,:)),'o--','markersize',10,'color',CO(1,:),'linewidth',3);hold on;

plot([0 5 8],fliplr(stab15(1,:)),'.:','markersize',40,'color',CO(1,:),'linewidth',3);hold on;

plot([0 5 8],fliplr(stab15(2,:)),'o:','markersize',10,'color',CO(1,:),'linewidth',3);hold on;
set(gca,'YScale','log');
axis([-0.5 8.5 0.5 10]);
set(gca,'YScale','log');
set(gca,'fontsize',16);
set(gcf,'position',[100 100 300 300]);
xticks([0 5 8]);
xticklabels({'8','3','0'});
xlabel('LAC concnetration (mM)','fontsize',24);
ylabel('decoupling coefficient','fontsize',24);
saveas(gcf,'SerialDilutionStability.fig');
saveas(gcf,'SerialDilutionStability.png');

function dydt=TwoStrains(t,y)
global mu eta D lambda kappa;
s1=y(1);
s2=y(2);
p1=y(3);
p2=y(4);
st=y(1)+y(2);
dydt=[mu(1)*s1/(s1+lambda(1)*p1)*s1*(1-st)-D*s1;
    mu(2)*s2/(s2+lambda(2)*p2)*s2*(1-st)-D*s2;
    mu(1)/(1+lambda(1))*p1*(1-st)+(s1-p1)*(eta(1,1)*p1+eta(2,1)*p2)-(kappa(1)+D)*p1;
    mu(2)/(1+lambda(2))*p2*(1-st)+(s2-p2)*(eta(1,2)*p1+eta(2,2)*p2)-(kappa(2)+D)*p2];
end

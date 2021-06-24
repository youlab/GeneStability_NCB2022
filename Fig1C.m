clear;
clc;
beta=0.5; % define the beta value; 
epsilons=[2 5 20];% define the effective plasmid transfer rate, epsilon

s1=0:0.01:1;% changing s1, which represents the community composition

for i=1:length(epsilons)
    epsilon=epsilons(i);
    pt=1-1/epsilon*1./((1-beta)*s1+beta);% for each epsilon, calculating the plasmid abundance pt as a function of s1
    patch([s1 fliplr(s1)], [min(pt)*ones(1,length(pt)) fliplr(pt)],[0.9 0.9 0.9],'LineStyle','none');hold on;
    plot(s1,pt,'k-','linewidth',5);
    hold on
end
set(gca,'fontsize',16);
xlabel('s_1','fontsize',24);
ylabel('plasmid abundance p_t','fontsize',24);
set(gcf,'position',[100 100 300 300]);
box on;



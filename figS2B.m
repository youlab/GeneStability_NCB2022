clear;
clc;
beta=0.5;
theta=0.1;
lambdas=[0.9 0.4 0.2];%define the value of lamdba1, the selection coefficient of this gene in species 1
s1=0:0.01:1;% s1 changes between 0 and 1
for i=1:length(lambdas)
    lambda=lambdas(i);
    pt=1-((beta*theta-theta)*s1+theta)/beta/lambda;% calculate the gene abundance pt as a function of s1
    patch([s1 fliplr(s1)], [min(pt)*ones(1,length(pt)) fliplr(pt)],[0.9 0.9 0.9],'LineStyle','none');hold on;
    plot(s1,pt,'k-','linewidth',5);
    hold on
end
set(gca,'fontsize',16);
xlabel('s_1','fontsize',24);
ylabel('plasmid abundance p_t','fontsize',24);
set(gcf,'position',[100 100 300 300]);
box on;



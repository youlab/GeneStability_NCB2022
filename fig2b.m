clear;
clc;
beta=0.5;
eps=[2 5 20];
s1=0:0.01:1;
for i=1:length(eps)
    ep=eps(i);
    y=1-1/ep*1./((1-beta)*s1+beta);
    patch([s1 fliplr(s1)], [min(y)*ones(1,length(y)) fliplr(y)],[0.9 0.9 0.9],'LineStyle','none');hold on;
    plot(s1,y,'-','color',[255, 87, 51]/256,'linewidth',5);
    hold on
end
set(gca,'fontsize',16);
xlabel('s_1','fontsize',24);
ylabel('plasmid abundance p_t','fontsize',24);
set(gcf,'position',[100 100 300 300]);
box on;
saveas(gcf,'fig2b.fig');
saveas(gcf,'fig2b.png');


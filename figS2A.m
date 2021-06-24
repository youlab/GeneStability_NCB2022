clear;
clc;
fp=10.^[-0.000001:-0.01:-2.5];%define the range of gene abundance fp
ms=10.^[2 4 6];%define the number of species in the local community
theta=1;%define the theta value
for i=1:length(ms)
    m=ms(i);
    stability=sqrt(m*fp./(theta-fp));%calculate functional stability as a function of gene abundance
    plot(fp,stability,'k-','linewidth',5);hold on;
end

set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'fontsize',16);
xlabel('gene relative abundance','fontsize',24);
ylabel('functional stability','fontsize',24);
axis(10.^[-2.5 0 -1 6]);
set(gcf,'position',[100 100 350 350]);


close;
clear;
s1=rand(300,1);
eps=0:0.01:6;
cv=0;
fp=0;
for i=1:length(eps)
    ep=eps(i);
    y=1-2./(s1+1)/ep;
    sp(i)=std(y.*(y>=0));
    fp(i)=mean(y.*(y>=0));
end
plot(eps,fp./sp,'-','color',[255, 87, 51]/256,'linewidth',5);hold on;
set(gca,'fontsize',16);
xlabel('relative transfer rate','fontsize',24);
ylabel('functional stability','fontsize',24);
set(gcf,'position',[100 100 400 400]);
axis([1 max(eps) min(fp./sp) max(fp./sp)]);
saveas(gcf,'fig2c.fig');
saveas(gcf,'fig2c.png');
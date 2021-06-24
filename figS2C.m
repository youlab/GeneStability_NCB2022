close;
clear;
s1=rand(300,1);% creat 300 parallel communities; s1 fluctuates between 0 and 1 following uniform distribution
lambdas=0.1:0.001:0.9;% lambda1 changes from 0.1 to 0.9
beta=0.5;
theta=0.1;
fp=0;%initialize fp, the mean gene abundance
sp=0;%initialize sp, the std of gene abundance
for i=1:length(lambdas)
    lambda=lambdas(i);
    p1=(1-theta/lambda)*s1;%calculate the gene abundance in species 1
    p2=(1-theta/lambda/beta)*(1-s1);%calculate the gene abundance in species 2
    pt=p1.*(p1>=0)+p2.*(p2>=0);% combine the gene abundance in 2 species; gene abundance must be non-negative
    sp(i)=std(pt);
    fp(i)=mean(pt);
end
    
plot(fp,fp./sp,'-','color',[255, 87, 51]/256,'linewidth',5);hold on;%plot functional stability as a funciton of mean gene abundance fp
set(gca,'fontsize',16);
xlabel('fp','fontsize',24);
ylabel('functional stability','fontsize',24);
set(gcf,'position',[100 100 400 400]);
axis([min(fp) max(fp) min(fp./sp) max(fp./sp)]);

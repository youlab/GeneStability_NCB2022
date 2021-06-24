close;
clear;
clc;
s1=rand(300,1);%s1 fluctuates between 0 and 1 following uniform distribution; create 300 parallel communities
epsilons=1:0.05:100; % increasing effective plasmid transfer rate from 1 to 100;
beta=0.5;%define the beta value
fp=0;%initialize the plasmid mean abundance
sp=0;%initialize the standard deviation of plasmid abundance
for i=1:length(epsilons)
    epsilon=epsilons(i);
    pt=1-1./((1-beta)*s1+beta)/epsilon;%calculating the plasmid abundance pt as a function of s1
    sp(i)=std(pt.*(pt>=0));%calculate the standard deviation of plasmid abundance; plasmid abundance must be non-negative.
    fp(i)=mean(pt.*(pt>=0));%calculate the mean plasmid abundance
end
    
plot(fp,1./(sp./fp),'-','color',[255, 87, 51]/256,'linewidth',5);hold on;%plot stability as a function of mean plasmid abundance
set(gca,'fontsize',16);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('fp','fontsize',24);
ylabel('functional stability','fontsize',24);
set(gcf,'position',[100 100 400 400]);
axis([min(fp) max(fp) min(1./(sp./fp)) max(1./(sp./fp))]);

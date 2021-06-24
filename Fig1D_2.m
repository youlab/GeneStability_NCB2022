close;
clear;
num=300;%define the number of parallel communities
s1=rand(1,num);% s1 fluctuates between 0 and 1 following uniform distribution
temp=[1:num];

subplot(4,1,1);% visualize the compositonal variations
patch([temp fliplr(temp)], [0*temp fliplr(s1)],[0 175 80]/256,'LineStyle','none');hold on;
patch([temp fliplr(temp)], [s1 ones(1,num)],[239 148 85]/256,'LineStyle','none');hold on;
axis([1 num 0 1]);

epsilon=5;%define the effective plasmid transfer rate
beta=0.5;%define the beta value
pt=1-1./((1-beta)*s1+beta)/epsilon;%calculating the plasmid abundance pt as a function of s1
pt=pt.*(pt>=0);%plasmid abundance must be non-negative.

subplot(4,1,2);
plot(1:num,pt./mean(pt),'k-');hold on;% plot the fluctuations of plasmid abundances normalized by its mean abundance
axis([0 num 0 2]);

epsilon=2;% repeat the calculation with another value of effective plasmid transfer rate
beta=0.5;%define the beta value
pt=1-1./((1-beta)*s1+beta)/epsilon;%calculating the plasmid abundance pt as a function of s1
pt=pt.*(pt>=0);%plasmid abundance must be non-negative.

subplot(4,1,3);
plot(1:num,pt./mean(pt),'k-');hold on;% plot the fluctuations of plasmid abundances normalized by its mean abundance
axis([0 num 0 2]);

epsilon=1.2;% repeat the calculation with another value of effective plasmid transfer rate
beta=0.5;%define the beta value
pt=1-1./((1-beta)*s1+beta)/epsilon;%calculating the plasmid abundance pt as a function of s1
pt=pt.*(pt>=0);%plasmid abundance must be non-negative.

subplot(4,1,4);
plot(1:num,pt./mean(pt),'k-');hold on;% plot the fluctuations of plasmid abundances normalized by its mean abundance
axis([0 num 0 2]);

set(gcf,'position',[100 100 200 400]);

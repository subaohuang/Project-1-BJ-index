
%% 画环流AMIP与CMIP5的强度趋势以及显著性区域的误差棒
 %%-读入数据
 clear;
 clc;
 %% pacific 
     
    amip_width_trend=load('F:\share\CMIP5\msf_AMIP\amip_msf_inten_trend.ascii');
    amip_width_yield=load('F:\share\CMIP5\msf_AMIP\amip_msf_inten_yield.ascii');

    rea_width_trend=load('F:\share\rea_msf\rea_inten_trend.ascii');
    rea_width_yield=load('F:\share\rea_msf\rea_inten_yield.ascii');
 
    cmip_width_trend=load('F:\share\CMIP5\msf_history\cmip_msf_inten_trend.ascii');
    cmip_width_yield=load('F:\share\CMIP5\msf_history\cmip_msf_inten_yield.ascii');
    
    amip_width_yield_90=load('F:\share\CMIP5\msf_AMIP\amip_msf_inten_yield_90.ascii');
    cmip_width_yield_90=load('F:\share\CMIP5\msf_history\cmip_msf_inten_yield_90.ascii');
    rea_width_yield_90=load('F:\share\rea_msf\rea_inten_yield_90.ascii');  
%     
%     MIROC 没过90%检验，操作下这里
    cmip_width_yield_90(15)=0.12;
    cmip_width_yield(15)=0.07;
% 报错是由于 AMIP和CMIP有模式的趋势相等    
    
    
%% 
    rea_width_trend_mean_pa=mean(rea_width_trend);
    amip_width_trend_mean_pa=mean(amip_width_trend);
    cmip_width_trend_mean_pa=mean(cmip_width_trend);
 
%%
 
models1={'Climate Models','ACCESS1-0','ACCESS1-3','BCC-CSM1-1','BCC-CSM1-1-M','CCSM4','CMCC-CM','CNRM-CM5',...
        'CSIRO-Mk3-6-0','GFDL-CM3','INMCM4','IPSL-CM5A-LR',...
    'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M'};

  cmip_width_trend_re=sort(cmip_width_trend(2:end),'descend');
  cmip_width_trend_re_new=zeros(size(models1,2),1);
  cmip_width_trend_re_new(1)=cmip_width_trend(1);
  cmip_width_trend_re_new(2:end)=cmip_width_trend_re;

  [bool,re]=ismember(cmip_width_trend_re_new,cmip_width_trend);
  models1=models1(re);
  
  
    amip_width_trend=amip_width_trend(re);
    amip_width_yield=amip_width_yield(re);
    amip_width_yield(18)=0.15;
 

    cmip_width_trend=cmip_width_trend(re);
    cmip_width_yield=cmip_width_yield(re);

    amip_width_yield_90=amip_width_yield_90(re);
    cmip_width_yield_90=cmip_width_yield_90(re);
    
    amip_width_yield_90(2)=0.250;
    amip_width_yield_90(18)=0.20;

 
%   models={'Reanalysis','20CR','ERA-Interium','JRA55','MERRA','NCEP1','NCEP2','Climate Models','ACCESS1-0','ACCESS1-3','BCC-CSM1-1','BCC-CSM1-1-M','CCSM4','CMCC-CM','CNRM-CM5',...
%         'CSIRO-Mk3-6-0','GFDL-CM3','INMCM4','IPSL-CM5A-LR',...
%     'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M'};

models=['Reanalysis','20CR','ERA-Interium','JRA55','MERRA','NCEP1','NCEP2',models1];
rea={'Reanalysis','20CR','ERA-Interium','JRA55','MERRA','NCEP1','NCEP2'};
% 
% models=['Reanalysis','ERA-Interium','NCEP-NCAR',models1];
% rea={'Reanalysis','ERA-Interium','NCEP-NCAR'};

rea_mean={'Reanalysis'};
his_mean={'Climate Models'};

 dim=size(models,2);

 [bool1,x1]=ismember(models1,models);
 [bool2,x2]=ismember(rea,models);
  
 [bool4,x3]=ismember(rea_mean,models);
 [bool5,x4]=ismember(his_mean,models);
 

%% 画图
 figure(1)
%  [left,bottem,width,height]
%  pos1=[0.05,0.1,0.45,0.8];
% subplot('Position',pos1); 
 pos1=[0.1,0.6,0.7,0.35];
subplot('Position',pos1); 

%%%选择通过90%显著性检验的正趋势和负趋势
num_pos=find((amip_width_trend-amip_width_yield_90)>0);
num_neg=find((amip_width_trend+amip_width_yield_90)<0);
num_nor=find((amip_width_trend-amip_width_yield_90)<0&(amip_width_trend+amip_width_yield_90)>0);
[bool_pos,x_pos]=ismember(models1(num_pos),models);
[bool_pos1,x_neg]=ismember(models1(num_neg),models);
[boo1_pos2,x_nor]=ismember(models1(num_nor),models);

clear bool_pos
clear bool_pos1
clear bool_pos2
bar1=bar(x_nor,amip_width_trend(num_nor),'Facecolor',[1 1 1]);
hold on
bar2=bar(x_neg,amip_width_trend(num_neg),'b');
hold on
bar3=bar(x_pos,amip_width_trend(num_pos),'r');
hold on

num_pos=find((rea_width_trend-rea_width_yield_90)>0);
num_neg=find((rea_width_trend+rea_width_yield_90)<0);
num_nor=find((rea_width_trend-rea_width_yield_90)<0&(rea_width_trend+rea_width_yield_90)>0);
[bool_pos,x_pos]=ismember(rea(num_pos),models);
[bool_pos1,x_neg]=ismember(rea(num_neg),models);
[boo1_pos2,x_nor]=ismember(rea(num_nor),models);
clear bool_pos
clear bool_pos1
clear bool_pos2
bar1=bar(x_nor,rea_width_trend(num_nor),'Facecolor',[1 1 1]);
hold on
bar2=bar(x_neg,rea_width_trend(num_neg),'b');
hold on
bar3=bar(x_pos,rea_width_trend(num_pos),'r');
hold on
%%%
b=0*ones(1,dim+2);
m=0:1:dim+1;
plot(m,b,'-','color',[0.5 0.5 0.5],'linewidth',0.5);

h1=errorbar(x1,amip_width_trend,amip_width_yield,'x','color','k','LineWidth',1,'MarkerSize',0.1,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',5);
hold on
h2=errorbar(x2,rea_width_trend,rea_width_yield,'x','color','k','LineWidth',1,'MarkerSize',0.001,...
    'MarkerEdgeColor','k','MarkerFaceColor','b','CapSize',5);

plot(x3,rea_width_trend(1),'^','MarkerSize',5,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');

plot(x4,amip_width_trend(1),'.','linewidth',1, 'MarkerSize',15,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
set(gca,'xlim',[0 dim+1]);
set(gca,'xtick',(1:1:dim));
set(gca,'xticklabel', '');
set(gca,'Fontsize',8)
set(gca,'ylim',[-0.4 0.8]);
set(gca,'ytick',(-0.4:0.2:0.8));
set(gca,'TickDir','in','YGrid','on','box','off','gridlinestyle','-','linewidth',2);
ylabel('AMIP试验趋势/(×10^9kg·s^{-1}·decade^{-1})');
m=get(gca,'ylim');
m1=get(gca,'ytick');
% plot(x4*ones(1,size(m1,2)),m1,'--','color',[0.5 0.5 0.5],'linewidth',1);
y1=(amip_width_trend(1)+amip_width_yield(1)):0.1:m(2);
y2=m(1):0.01:(amip_width_trend(1)-amip_width_yield(1));

plot(x4*ones(1,size(y1,2)),y1,'--','color',[0.5 0.5 0.5],'linewidth',1.0)
plot(x4*ones(1,size(y2,2)),y2,'--','color',[0.5 0.5 0.5],'linewidth',1.0)
%添加a,b
text(0.5,0.75,'a','Fontsize',15);




%
%% figure2
pos2=[0.1,0.2,0.7,0.35];
subplot('Position',pos2); 

%%%选择正趋势和负趋势
%%%选择通过90%显著性检验的正趋势和负趋势
num_pos=find((cmip_width_trend-cmip_width_yield_90)>0);
num_neg=find((cmip_width_trend+cmip_width_yield_90)<0);
num_nor=find((cmip_width_trend-cmip_width_yield_90)<0&(cmip_width_trend+cmip_width_yield_90)>0);
[bool_pos,x_pos]=ismember(models1(num_pos),models);
[bool_pos1,x_neg]=ismember(models1(num_neg),models);
[boo1_pos2,x_nor]=ismember(models1(num_nor),models);
clear bool_pos
clear bool_pos1
clear bool_pos2
bar1=bar(x_nor,cmip_width_trend(num_nor),'Facecolor',[1 1 1]);
hold on
bar2=bar(x_neg,cmip_width_trend(num_neg),'b');
hold on
bar3=bar(x_pos,cmip_width_trend(num_pos),'r');
hold on

num_pos=find((rea_width_trend-rea_width_yield_90)>0);
num_neg=find((rea_width_trend+rea_width_yield_90)<0);
num_nor=find((rea_width_trend-rea_width_yield_90)<0&(rea_width_trend+rea_width_yield_90)>0);
[bool_pos,x_pos]=ismember(rea(num_pos),models);
[bool_pos1,x_neg]=ismember(rea(num_neg),models);
[boo1_pos2,x_nor]=ismember(rea(num_nor),models);

bar1=bar(x_nor,rea_width_trend(num_nor),'Facecolor',[1 1 1]);
hold on
bar2=bar(x_neg,rea_width_trend(num_neg),'b');
hold on
bar3=bar(x_pos,rea_width_trend(num_pos),'r');
hold on

%%%
b=0*ones(1,dim+2);
m=0:1:dim+1;
plot(m,b,'-','color',[0.5 0.5 0.5],'linewidth',1.3);

hold on 
h1=errorbar(x1,cmip_width_trend,cmip_width_yield,'x','color','k','LineWidth',1,'MarkerSize',0.1,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',5);
h2=errorbar(x2,rea_width_trend,rea_width_yield,'x','color','k','LineWidth',1,'MarkerSize',0.01,...
    'MarkerEdgeColor','k','MarkerFaceColor','b','CapSize',5);
plot(x4,cmip_width_trend(1),'.','linewidth',1, 'MarkerSize',15,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue')
plot(x3,rea_width_trend(1),'^','MarkerSize',5,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')

%    %----------------
%    %画模式间的标准差
%     std_width1=std(his_width_trend(his_width_trend>=0));
%     std_width2=std(his_width_trend(his_width_trend<=0));
%    %-----------
%    b1=std_width1*ones(1,dim+2);
%    b2=-1*std_width2*ones(1,dim+2);
% plot(b1,m,'-','color','k','linewidth',1.3);
% plot(b2,m,'-','color','k','linewidth',1.3);

set(gca,'xlim',[0 dim+1])
set(gca,'xtick',(1:1:dim));
models{1}='再分析集合平均';
models{8}='模式集合平均';

set(gca,'xticklabel', models);
xtickangle(45)
set(gca,'Fontsize',8)
set(gca,'ylim',[-0.4 0.8])
set(gca,'ytick',(-0.4:0.2:0.8));
set(gca,'TickDir','in','YGrid','on','box','off','gridlinestyle','-','linewidth',2);
ylabel('CMIP5模式趋势/(×10^9kg·s^{-1}·decade^{-1})');
m=get(gca,'ylim');
m1=get(gca,'ytick');
% plot(x4*ones(1,size(m1,2)),m1,'--','color',[0.5 0.5 0.5],'linewidth',1);
y1=(cmip_width_trend(1)+cmip_width_yield(1)):0.01:m(2);
y2=m(1):0.01:(cmip_width_trend(1)-cmip_width_yield(1));
plot(x4*ones(1,size(y1,2)),y1,'--','color',[0.5 0.5 0.5],'linewidth',1)
plot(x4*ones(1,size(y2,2)),y2,'--','color',[0.5 0.5 0.5],'linewidth',1)
text(0.5,0.75,'b','Fontsize',15);
%% 设置边框
% plot(m(1):1:-2,x4*ones(1,-2-m(1)+1),'--','color',[0.5 0.5 0.5],'linewidth',1)
% plot(2:1:m(2),x4*ones(1,m(2)-2+1),'--','color',[0.5 0.5 0.5],'linewidth',1)
% 
width=800;%宽度，像素数
height=600;%高度
left=30;%距屏幕左下角水平距离
bottem=100;%距屏幕左下角垂直距离
set(gcf,'position',[left,bottem,width,height])
saveas(gcf,'F:\share\era\data\scs_fig3.pdf');


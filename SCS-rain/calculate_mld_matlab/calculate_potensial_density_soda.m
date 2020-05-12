clear;
clc;
%% 读入温度数据
nc1=netcdf('F:\share\SODA\SODA3.4.2_temp_5m.nc');
lat=nc1{'lat'}(:);
lon=nc1{'lon'}(:);
temp=nc1{'temp'}(:,:,:);
FillValue=nc1{'temp'}.FillValue_(:);
temp(temp==FillValue)=nan;
close(nc1)
%% 读入盐度数据
nc2=netcdf('F:\share\SODA\soda3.4.2_salt_mn_1980-2015_chazhi_new.nc');
salt=nc2{'salt'}(1:408,:,:);
FillValue=nc2{'salt'}.FillValue_(:);
salt(salt==FillValue)=nan;
close(nc2)

%% 计算密度
 dim=size(salt);
 P=zeros(dim)+10;
 %1000hPa=100000pa = 10db 1db= 0.1bar = 10000pa
 PR=zeros(dim);%reference pressure
 T_diff=temp-0.8;
 pden=sw_pden(salt,T_diff,P,PR);
 pden_clm= squeeze(nanmean(pden,1));
 pden_clm=pden_clm-1000.0;
 %% 画图
 pden1=squeeze(pden(1,:,:));
 pden1=pden1-1000.0;
 [lonn,latt]=meshgrid(lon,lat);
m_proj('Equidistant Cylindrical','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
hold on
% 
colormap('jet')
[c,h]=m_contourf(lonn,latt,pden_clm,20:1:30);  % 选择等值线间隔和范围
clabel(c,h);    %添加等值线的label
% m_pcolor(lonn,latt,pden1);
m_coast('patch',[0 0 0],'edgecolor','k');
set(gca,'Ytick',[min(lat):5:max(lat)]);
m_grid('linewi',1,'linest','none','box','on','tickdir','out','yaxisloc','left');
 %%  写入nc文件：位势温度

d0=408;
d1=size(lon,1);
d2=size(lat,1);
time=1:1:408;

ncid = netcdf.create('F:\share\SODA\SODA3.4.2_potensial_density_mld_0.8.nc','CLOBBER');

%%%定义维度
time_dim = netcdf.defDim(ncid,'time',d0);
lat_dim = netcdf.defDim(ncid,'lat',round(d2));
lon_dim = netcdf.defDim(ncid,'lon',round(d1));

%%%定义变量
varid1= netcdf.defVar(ncid,'lat','double',lat_dim);
varid2= netcdf.defVar(ncid,'lon','double',lon_dim);
varid3= netcdf.defVar(ncid,'time','double',time_dim);
varid4 = netcdf.defVar(ncid,'phro_mld','double',[time_dim lat_dim lon_dim]);
netcdf.endDef(ncid);%定义结束
%%%赋值
netcdf.putVar(ncid,varid1,lat);
netcdf.putVar(ncid,varid2,lon);
netcdf.putVar(ncid,varid3,time);
netcdf.putVar(ncid,varid4,pden);
netcdf.reDef(ncid);%赋值结束

%%属性
netcdf.putAtt(ncid,varid1,'units','degrees_north');
netcdf.putAtt(ncid,varid2,'units','degrees_east');
netcdf.putAtt(ncid,varid3,'units','1');
netcdf.putAtt(ncid,varid1,'long_name','latitude');
netcdf.putAtt(ncid,varid2,'long_name','longitude');
netcdf.putAtt(ncid,varid3,'long_name','time');
netcdf.putAtt(ncid,varid4,'long_name','potensial_density_at_mld');
netcdf.endDef(ncid);
netcdf.close(ncid);


 
 
 
 
 
function m2y(data)   ;;;
begin
 size=dimsizes(data)
 year=size(0)/12
 print(year)
 sst_annual=new((/year,size(1),size(2)/),typeof(data))
 do i=0,year-1
 sst_annual(i,:,:)=dim_avg_n(data(i*12:i*12+11,:,:),0)
 end do
  sst_annual!0="year"
  sst_annual!1="lat"
  sst_annual!2="lon"
  sst_annual&year=ispan(1, year, 1)
  sst_annual&lat=data&lat
  sst_annual&lon=data&lon
return(sst_annual)
end 


begin
;;;读入hadisst
  diri="/mnt/f/share/GISS/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"air.2x2.1200.mon.anom.comb.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)

  timeind=ind(date.le.201212.and.date.ge.195801)
  sst=f1->air(timeind,{60:-60},:)
  sst=(/rmMonAnnCycTLL(sst)/)

  ;;;mask掉降水和风场的陆地数据
  fils2=systemfunc("ls /mnt/f/share/shixiong/landsea.nc")
  a=addfile(fils2, "r")
  lsdata = a->LSMASK
  ; lsm  = landsea_mask(lsdata,p_sm_lin&LAT,p_sm_lin&LON)   ;;得到的结果，1表示land   
  ; p_sm_lin = mask(p_sm_lin,lsm.eq.1,False)      

  lsm1 = landsea_mask(lsdata,sst&lat,sst&lon)   ;;得到的结果，1表示land   
  sst_land=sst
  sst_sea=sst

  sst_land = mask(sst,lsm1.eq.0,False)
  sst_sea = mask(sst,lsm1.eq.1,False)     
  
  printVarSummary(sst_land)
  printVarSummary(sst_sea)
  


;;;
  sst_annual_sea=m2y(sst_sea)
  printVarSummary(sst_annual_sea)
  
  sst_annual_sea&year=ispan(1,dimsizes(sst_annual_sea(:,1,1)),1)
  sst_sm_sea=runave_n_Wrap(sst_annual_sea, 3, 1, 0)

   sst_annual_land=m2y(sst_land)
  sst_annual_land&year=ispan(1,dimsizes(sst_annual_land(:,1,1)),1)
  sst_sm_land=runave_n_Wrap(sst_annual_land, 3, 1, 0)


  land_col=reshape(sst_annual_land(lat|:,lon|:,year|:), (/dimsizes(sst&lat)*dimsizes(sst&lon),\
                                     dimsizes(sst_annual_sea&year)/))
  
  sea_col=reshape(sst_annual_sea(lat|:,lon|:,year|:), (/dimsizes(sst&lat)*dimsizes(sst&lon),\
                                   dimsizes(sst_annual_sea&year)/))
  printVarSummary(land_col)
  
  
  
  land_col=where(ismissing(land_col), 999.0, land_col)
  sea_col=where(ismissing(sea_col), 999.0, sea_col)
   
   
   opt1=True
   opt1@fout=diri_pi+"/land_col.txt"

   opt2=True
   opt2@fout=diri_pi+"/sea_col.txt"
   write_matrix(land_col, "55f7.2",opt1)
   write_matrix(sea_col, "55f7.2", opt2)
  
   exit()

  dsize=dimsizes(land_col)

  ncols=dsize(0)
  nsvd=3
 
  homlft = new((/nsvd,ncols/),float)
  hetlft = new((/nsvd,ncols/),float)
  homrgt = new((/nsvd,ncols/),float)
  hetrgt = new((/nsvd,ncols/),float)

  x= svdstd(sea_col,land_col, nsvd, homlft, hetlft, homrgt, hetrgt)
  printVarSummary(x)


  end 
;;SVD分析
begin
;;;读入hadisst
  diri="/mnt/f/share/ersst/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201212.and.date.ge.195801)
  sst=lonFlip(f1->sst(timeind,{60:-60},:))
  ; sst=(/rmMonAnnCycTLL(sst)/)
  sst_annual=month_to_annual(sst, 1)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  printVarSummary(sst_annual)
  sst_sm=runave_n_Wrap(sst_annual, 3, 1, 0)
  printVarSummary(sst_sm)
  
  ;;插值到2.5乘2.5网格
   yo=fspan(-59.5,58,48)
   xo=fspan(0.5, 358, 144)
   
  sst_sm_lin=linint2_Wrap(sst_sm&longitude, sst_sm&latitude(::-1), sst_sm(:,::-1,:), \
                          True, xo, yo, 0)
  sst_sm_lin=(/dim_standardize_n_Wrap(sst_sm_lin, 0, 0)/)

  printVarSummary(sst_sm_lin)
  delete(diri)  
  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)
  delete(sst)

;;;读入降水
  diri="/mnt/f/share/era/data/"
  fils1=systemfunc("ls /mnt/f/share/NOAA/precip.mon.anom.nc") 
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201212.and.date.ge.195801)
  p=f1->precip(timeind,{60:-60},:)
  
  ; p=(/rmMonAnnCycTLL(p)/)
  p_annual=month_to_annual(p, 1)
  p_annual&year=ispan(1,dimsizes(p_annual(:,1,1)),1)
  printVarSummary(p_annual)
  p_sm=runave_n_Wrap(p_annual, 3, 1, 0)
  printVarSummary(p_sm)

  p_sm_lin=linint2_Wrap(p_sm&lon, p_sm&lat(::-1), p_sm(:,::-1,:), \
                        True, xo, yo, 0)
  p_sm_lin=(/dim_standardize_n_Wrap(p_sm_lin, 0, 0)/)
  printVarSummary(p_sm_lin)
  
  delete(diri)  
  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)



  sst_col=reshape(sst_sm_lin(LATITUDE|:,LONGITUDE|:,year|:), (/dimsizes(p_sm&lat)*dimsizes(p_sm&lon),\
                                     dimsizes(p_sm&year)/))
  
  p_col=reshape(p_sm_lin(LAT|:,LON|:,year|:), (/dimsizes(p_sm&lat)*dimsizes(p_sm&lon),\
                                   dimsizes(p_sm&year)/))
  printVarSummary(sst_col)
  printVarSummary(p_col)


  dsize=dimsizes(sst_col)

  ncols=dsize(0)
  nsvd=3
 
  homlft = new((/nsvd,ncols/),float)
  hetlft = new((/nsvd,ncols/),float)
  homrgt = new((/nsvd,ncols/),float)
  hetrgt = new((/nsvd,ncols/),float)

  x= svdstd(sst_col,p_col, nsvd, homlft, hetlft, homrgt, hetrgt)
  printVarSummary(x)

  end
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

begin
  ; yearstart = 200001
  ; yearend = 201912
    yearstart = 195001
  yearend = 201912

  lat1 = -30
  lat2 = 30
  lon1 = 0
  lon2 = 360

 
   ;;读入HadiSST海温数据
    ; diri="/home/yangsong3/data-observation/linshh/data/sst/"
    ; fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
     fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201912.and.date.ge.197901)

    sst  = f1->sst(timeind,:,:)
    ; sst  = short2flt(f1->sst(timeind,:,:))

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    dsst = center_finite_diff_n(ssta, 1, False, 1, 0)
    copy_VarCoords(sst, dsst)


    sst_nep = dim_avg_n_Wrap(ssta(:,{0:20},{240:280}), (/1,2/))
    dsst_nep = dim_avg_n_Wrap(dsst(:,{0:20},{240:280}), (/1,2/))
    


   mxlag    = 12
   
  x = sst_nep
  y = dsst_nep

  x_Lead_y = esccr(x,y,mxlag)
  y_Lead_x = esccr(y,x,mxlag)    ; switch the order of the series

  ccr = new ( 2*mxlag+1, float)    
  ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
  ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"


  wks=gsn_open_wks("pdf", "sst和dsst超前滞后相关")
  resp=True
  resp@gsnXRefLine = 0
  resp@trYMaxF=1.0
  resp@trYMinF=-1.0
  plot = gsn_csm_xy(wks,fspan(-12,12,25),ccr,resp)
  end 


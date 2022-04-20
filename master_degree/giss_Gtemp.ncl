begin
;;;读入hadisst
  diri="/mnt/f/share/GISS/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"air.2x2.1200.mon.anom.comb.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201212.and.date.ge.197901)
  sst=f1->air(timeind,:,:)

;;;
  ; sst=(/rmMonAnnCycTLL(sst)/)
  sst_annual=month_to_annual(sst, 1)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
;;;


;;;

    ;;;global mean temp

  sst_avg=dim_avg_n_Wrap(sst_annual, (/1,2/))
  sst_avg=(/dim_standardize_n_Wrap(sst_avg, 0, 0)/)
  sst_avg_sm=runave_n_Wrap(sst_avg, 3, 1, 0)
;;;
opt1=diri+"giss_gtemp_1979_2012.ascii"
asciiwrite(opt1, sprintf ("%6.2f,",sst_avg))
  ;;;mask掉陆地的温度
   fils2=systemfunc("ls /mnt/f/share/shixiong/landsea.nc")
   a=addfile(fils2, "r")
   lsdata = a->LSMASK
   lsm  = landsea_mask(lsdata,sst_annual&lat,sst_annual&lon)   ;;得到的结果，1表示land   
   sst_annual = mask(sst_annual,lsm.eq.1,False)   

;;;
  sst_ipo=(/dim_avg_n_Wrap(sst_annual(:,{5:-5},{210:270}), (/1,2/))\
    -  (dim_avg_n_Wrap(sst_annual(:,{45:25},{150:210}), (/1,2/))+dim_avg_n_Wrap(sst_annual(:,{-25:-45},{150:210}), (/1,2/)))/2\
    /)
    sst_ipo_sm=runave_n_Wrap(sst_ipo, 3, 1, 0)
    sst_ipo_sm=(/dim_standardize_n_Wrap(sst_ipo_sm, 0, 0)/)
    ; sst_ipo_sm=(/dtrend_n(sst_ipo_sm, False, 0)/)


;;;
 
  opt2=diri+"giss_ipo_1979_2012.ascii"
  asciiwrite(opt2, sprintf ("%6.2f,", sst_ipo))


  
;;;画图
  wks=gsn_open_wks("pdf", ""+diri_pi+"giss_GTMP")
  
  lineres=True 
  lineres@gsnDraw=False
  lineres@gsnFrame=False
  lineres@gsnYRefLine=0.0
  lineres@xyLineColors="blue"
  lineres@xyDashPatterns=0
 
  lineres@tmXTOn=False
  lineres@tmYROn=False
  lineres@tmXBMode="Explicit"
  lineres@tmXBMinorValues=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  lineres@tmXBValues=ispan(1,dimsizes(sst_annual(:,1,1)),5)
  lineres@tmXBLabels=ispan(1980,2012,5)
  lineres@tmYLLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  lineres@vpHeightF=0.25
  lineres@vpWidthF=0.75
  
  plot=gsn_csm_xy(wks, ispan(0,dimsizes(sst_avg_sm)-1,1),(/sst_avg,sst_ipo_sm/), \
                     lineres)
  draw(plot)
  frame(wks)
  
end
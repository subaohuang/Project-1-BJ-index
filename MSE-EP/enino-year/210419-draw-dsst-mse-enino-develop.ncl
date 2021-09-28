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

function mse_diff(data,year1,year2,year3,year4)
begin
  ; rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  ; copy_VarCoords(data(1,:,:), rc1)
  ; rc1 = rc1*(year2-year1+1)
  rc2 = dim_avg_n_Wrap(data({year1:year2},:,:),0)
  rc2 = dim_avg_n_Wrap(data({year3:year4},:,:),0) - dim_avg_n_Wrap(data({year1:year2},:,:),0)
  return rc2
end 


begin
 yearstart = 195001
 yearend = 201912
    ;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ; fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")
    ; ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind1=ind(date.le.yearend.and.date.ge.yearstart)
    ; timeind2=ind(date.le.201012.and.date.ge.197901)
    sst  =lonFlip(f1->sst(timeind1,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ; sst_clm = clmMonTLL(f1->sst(timeind2,:,:))
    ; ssta = calcMonAnomTLL(sst, sst_clm)
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)


    dsst = center_finite_diff_n(ssta, 1, False, 1, 0)
    copy_VarCoords(ssta, dsst)
    
    dsst_nino1_2 = dim_avg_n_Wrap(dsst(:,{0:-10},{270:280}), (/1,2/))


 
   ;;;;cal season trend 
   

    dsst_MAM = month_to_season(dsst_nino1_2, "MAM")
    dsst_JJA = month_to_season(dsst_nino1_2, "JJA")
    dsst_SON = month_to_season(dsst_nino1_2, "SON")
    dsst_DJF = month_to_season(dsst_nino1_2, "DJF")


    dsst_DJF(0:dimsizes(dsst_DJF&time)-2) = dsst_DJF(1:)
    dsst_DJF(dimsizes(dsst_DJF&time)-1) = dsst_nino1_2(dimsizes(dsst_nino1_2&time)-1)

    ; dsst_JJA = (/(dsst_JJA + dsst_SON)*0.5/)
    
    ; dsst_JJA = (/dim_standardize_n(dsst_JJA, 1, 0)/)


    ssta_DJF = month_to_season(ssta, "DJF")
    ssta_DJF(0:dimsizes(ssta_DJF&time)-2,:,:) = ssta_DJF(1:,:,:)
    ssta_DJF(dimsizes(ssta_DJF&time)-1,:,:) = ssta(dimsizes(ssta&time)-1,:,:)

    nino3_4 = dim_avg_n_Wrap(ssta_DJF(:,{-5:5},{190:240}), (/1,2/))

    std = dim_stddev_n_Wrap(nino3_4, 0)


;;;读入MSE的数据
  f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")
  date:=cd_calendar(f1->time, 1)
  timeind:=ind(date.le.yearend.and.date.ge.yearstart)
  MSE = f1->MSE(timeind,:,:,:)

  ;;;MSE 
  plev = MSE&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

  MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

  MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

  MSE_col = (/MSE_col/10^6/)

 
  MSE_col_MAM = month_to_season(MSE_col, "MAM")
  MSE_col_JJA = month_to_season(MSE_col, "JJA")
  MSE_col_SON = month_to_season(MSE_col, "SON")
  MSE_col_DJF = month_to_season(MSE_col, "DJF")

  MSE_nep_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{0:20},{210:280}), (/1,2/))

  ; MSE_nep_mam = (/dim_standardize_n(MSE_nep_mam, 1, 0)/)

 ;;;计算趋势
  dsst_MAM&time = ispan(1950,2019,1)
  dsst_JJA&time = ispan(1950,2019,1)
  dsst_SON&time = ispan(1950,2019,1)
  dsst_DJF&time = ispan(1950,2019,1)
  
  ; dsst_MAM&time = ispan(1900,2019,1)
  ; dsst_JJA&time = ispan(1900,2019,1)
  ; dsst_SON&time = ispan(1900,2019,1)
  ; dsst_DJF&time = ispan(1900,2019,1)

 
  year  = ispan(1950, 2019, 1)
  
  ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1
  ; ely = year(ind(nino3_4.ge.0.5*std))
   x = ispan(1,dimsizes(ely), 1)


  dsst_MAM_enino = dsst_MAM({ely})
  dsst_JJA_enino = dsst_JJA({ely})
  dsst_SON_enino = dsst_SON({ely})
  dsst_DJF_enino = dsst_DJF({ely})

  MSE_nep_mam&time = ispan(1950,2019,1)
  MSE_nep_mam_enino = MSE_nep_mam({ely})

 corr = escorc_n(MSE_nep_mam_enino, dsst_JJA_enino, 0, 0)
 print(corr)
 ; exit()
    
  wks   = gsn_open_wks("pdf", "HadIsst data-El Nnio发展年 NINO1+2 区 dsst&MSE的时间序列-1950-2019")

  res               = True                   ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False
  res@tiMainString  = " Time series of MAM MSE and Nino 1+2 JJA SSTA tendency"   ; title
  res@tiMainFontHeightF = 0.02
  ; res@tiYAxisFontHeightF = 0.02
  ; res@tiXAxisFontHeightF = 0.02
  ; res@tiYAxisString = "MSE anomaly"           ; y axis title
  res@tiXAxisString = "El Nino Year"

  res@xyLineColors      = (/"blue","red"/)  ; line colors
   res@xyMonoLineThickness = True
  res@xyLineThicknessF = (/4.0/)        ; line thicknesses
  res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
  res@gsnYRefLine = 0.0

  ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
  ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
  ; res@xyMarkers         = (/2,0/)          ; style

  res@tmXBMode="Explicit"

  res@trXMaxF = dimsizes(ely)+1
  ; res@trXMinF = 1950
  res@trXMinF = 0
  ; res@trYMinF = 0.2
  ; res@trYMaxF = 1.8

  res@tmXTOn = False
  res@tmYROn  = False

  res@tmXBValues= x
  res@tmXBLabels= ely
  res@tmXBLabelAngleF = -60.0
  res@tmXBLabelJust = "CenterCenter"
  ; res@tmXBLabelFontHeightF = 0.018
  res@tmXBLabelFontHeightF = 0.02
  res@tmXBMajorLengthF = 0.01
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmYRMode="Automatic"


  ; plot = gsn_csm_xy(wks, x, (/dsst_MAM_enino,dsst_JJA_enino,dsst_SON_enino,dsst_DJF_enino/), res)
  plot = gsn_csm_xy(wks, x, (/dsst_JJA_enino,MSE_nep_mam_enino/), res)


    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"JJA Nino1+2 SSTA tendency","MAM MSE"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot = simple_legend(wks,plot,gres,lineres,textres)   

draw(plot)
frame(wks)
end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
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
;;;;

function change(lhtfl_ct,lhtfl_sen, season)
begin 
    lhtfl_ct_annual = month_to_season(lhtfl_ct, season)
    lhtfl_sen_annual = month_to_season(lhtfl_sen, season)
    lhtfl_change = lhtfl_ct_annual(1,:,:)
    lhtfl_change  = (/dim_avg_n_Wrap(lhtfl_sen_annual, 0)-dim_avg_n_Wrap(lhtfl_ct_annual, 0)/)
    

return lhtfl_change
end     

;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =25601
    endmon = 28012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -30
   lat2 = 30
   lon1 = 0
   lon2 = 360

;;;;
 

;;;;model name 
   modelname = "B2000_alt_north_year_CTRL"
;;;;


;;;;读入控制实验结果
    diri="/home/ys17-19/lsh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".TEMP.h0.0251-0280.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    

    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})
    
;;;read NCEP-NCAR heat flux data
    fils3 = diri + modelname+".cam.h0.0251-0280_OLR_heatflux.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    
    lhtfl_ct = f3->LHFLX(timeind,{lat1:lat2},{lon1:lon2})
    shtfl_ct = f3->SHFLX(timeind,{lat1:lat2},{lon1:lon2})
    nswrs_ct = f3->FSNS(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs_ct = f3->FLNS(timeind,{lat1:lat2},{lon1:lon2})
;;;



;;;读入B22000-SST海温数据

    ;;;;model name 
     modelname = "B2000_alt_north_year"
    ;;;;
    diri="/home/ys17-19/lsh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".TEMP.h0.0251-0280.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    sst_sen = sst(:,{lat1:lat2},{lon1:lon2})
    
;;;read NCEP-NCAR heat flux data
    fils3 = diri+modelname+".cam.h0.0251-0280_OLR_heatflux.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    
    lhtfl_sen = f3->LHFLX(timeind,{lat1:lat2},{lon1:lon2})
    shtfl_sen = f3->SHFLX(timeind,{lat1:lat2},{lon1:lon2})
    nswrs_sen = f3->FSNS(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs_sen = f3->FLNS(timeind,{lat1:lat2},{lon1:lon2})
;;;

 printMinMax(lhtfl_sen, 1)
 printMinMax(shtfl_sen, 1)
 printMinMax(nswrs_sen, 1)
 printMinMax(nlwrs_sen, 1)


printVarSummary(lhtfl_sen)


  ;;;;choose season

  season = (/"MAM","JJA","SON","DJF"/)

  wks   = gsn_open_wks("pdf", "ATM-ocean-forcing-SST-annual-global-cesm-B2000")

  plot = new(12,graphic)
  str = "" + ispan(1,dimsizes(plot),1)
  str2 = "0" + ispan(1,dimsizes(plot),1)

 do i=0,3   
    sst_change = change(sst_ct,sst_sen, season(i))
    printVarSummary(sst_change)
    lhtfl_change = change(lhtfl_ct,lhtfl_sen, season(i))
    shtfl_change = change(shtfl_ct,shtfl_sen, season(i))
    nswrs_change = change(nswrs_ct,nswrs_sen, season(i))
    nlwrs_change = change(nlwrs_ct,nlwrs_sen, season(i))

    Do_change  = lhtfl_change
    Do_change = (/-1*( nswrs_change - nlwrs_change - lhtfl_change -shtfl_change)/) ;;;;Do = -Qnet
    

    lhtfl_clm  = dim_avg_n_Wrap(lhtfl_ct, 0)
    printVarSummary(lhtfl_clm)
    
    lhtfl_change_atm = lhtfl_clm
    lhtfl_change_atm  = (/lhtfl_change - 0.06*lhtfl_clm*sst_change/)

    Qa_change = lhtfl_change_atm 
    Qa_change = (/ nswrs_change - nlwrs_change -shtfl_change - lhtfl_change_atm/)


    T_change_atm = Do_change
    T_change_atm = (/Qa_change/(0.06*lhtfl_clm)/)

    T_change_ocean = Do_change
    T_change_ocean = (/Do_change/(0.06*lhtfl_clm)/)
    
    T_total = T_change_ocean
    T_total = (/T_change_atm +T_change_ocean/)
;;;;
   
print("*********开始画图********")
;;**************************************************************
;画图
;;***************************************************************


  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_jaisnd"
  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLineLabelInterval = 2
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillDrawOrder = "PostDraw"
  res@mpMinLatF             = -10
  res@mpMaxLatF             = 30
  res@mpMinLonF             =70
  res@mpMaxLonF             = 170
   ; lat1 = 0
   ; lat2 = 15
   ; lon1 = 110
   ; lon2 = 155

   ; lat1 = -30
   ; lat2 = 30
   ; lon1 = 0
   ; lon2 = 360

  ; res@mpMinLatF             = lat1
  ; res@mpMaxLatF             = lat2
  ; res@mpMinLonF             =lon1
  ; res@mpMaxLonF             = lon2

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
  res@gsnRightString = ""

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.1





  printVarSummary(sst_change)
  printVarSummary(T_change_atm)
  printVarSummary(T_change_ocean)

  res@gsnLeftString = "SST trend "
  plot(i*3) = gsn_csm_contour_map(wks, sst_change, res)
  
  ; delete(res@gsnSpreadColorStart)
  ; delete(res@gsnSpreadColorEnd)
  res@cnFillPalette = "BlRe"
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.2


  res@gsnLeftString = "Atm Forcing SST trend"
  plot(i*3+1) = gsn_csm_contour_map(wks, T_change_atm, res)
  res@gsnLeftString = "Ocean  Forcing SST trend"
  plot(i*3+2) = gsn_csm_contour_map(wks, T_change_ocean, res)


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
   


     plot@$str(i*3)$ = gsn_add_polyline(wks, plot(i*3), pgx, pgy, pgres)
     plot@$str(i*3+1)$ = gsn_add_polyline(wks, plot(i*3+1), pgx, pgy, pgres)
     plot@$str(i*3+2)$ = gsn_add_polyline(wks, plot(i*3+2), pgx, pgy, pgres)




  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
  plot@$str2(i*3)$ = gsn_add_polyline(wks, plot(i*3), pgx_pool, pgy_pool, pgres)
  plot@$str2(i*3+1)$ = gsn_add_polyline(wks, plot(i*3+1), pgx_pool, pgy_pool, pgres)
  plot@$str2(i*3+2)$ = gsn_add_polyline(wks, plot(i*3+2), pgx_pool, pgy_pool, pgres)
end do 

  
  pres = True
  pres@gsnPanelLabelBar  = False
  ; pres@gsnPanelFigureStrings = (/"a","b","c"/)
  gsn_panel(wks, plot, (/4,3/), pres)
 


end
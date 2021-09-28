; function  rc_test(ssta_rc)
; begin
;   tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
;    df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
;    b = tval    ; b must be same size as tval (and df)
;    b = 0.5
;    prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
;    copy_VarCoords(ssta_rc, prob)

; return prob
; end
;;;add ncl_to_py
function  rc_test(ssta_rc,id1,id2)
begin
   dim = dimsizes(id1)
   varname = getVarDimNames(id2)
   ;;有效自由度 
   r1_r = esacr_n(id1, 1,0)  
   r2_r = esacr_n(id2, 1,0)

   r1_o = r1_r(1)
   r2 = r2_r(1,:,:)  ;;;取lag1
   f = dimsizes(r2)

   r1  = onedtond(r1_o, dimsizes(r2))

   v = dim*(/1-r1*r2/)/(/1+r1*r2/) 
 
   ssta_rc_new = onedtond(ssta_rc, dimsizes(ssta_rc))
   tval = sqrt(v)*(/ssta_rc_new/)/(sqrt(1- ssta_rc_new^2))

   df   = v
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
   ;;或者用student_t
  ;; prob_1 = student_t(tval, df)x
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
    yearstart = 198001
    yearend = 201912

    year1 = 1980
    year2 = 2019

    lat1 = -30
    lat2 = 30
    lon1 = 0
    lon2 = 360
    
    runmean_time = 11

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

    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

    ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)

    ; MSE_col = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_col = (/MSE_col/10^6/)
    MSE_col_annual = month_to_annual(MSE_col, 1)
    MSE_col_annual!0 = "time"

    MSE_col_MAM = month_to_season(MSE_col, "MAM")
    MSE_col_JJA = month_to_season(MSE_col, "JJA")
    MSE_col_SON = month_to_season(MSE_col, "SON")
    MSE_col_DJF = month_to_season(MSE_col, "DJF")
    
    MSE_col_MAM = (/dim_standardize_n_Wrap(MSE_col_MAM, 0, 0)/)

    MSE_col_MAM = (/runave_n(MSE_col_MAM, runmean_time, 0, 0)/)

    MSE_nep_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{0:20},{240:280}), (/1,2/))

    ; MSE_tna_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{20:30},{270:300}), (/1,2/))

    ; MSE_ind_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{10:30},{50:100}), (/1,2/))
    
    
   
    ;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)

    sst  = lonFlip(f1->sst(timeind,{-50:50},:))

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    ; ssta = (/runave_n_Wrap(ssta, 9*12, 0, 0)/)
    ssta_mam = month_to_season(ssta, "MAM")

    ssta_mam = (/dim_standardize_n_Wrap(ssta_mam, 0, 0)/)

    ssta_mam = (/runave_n(ssta_mam, runmean_time, 0, 0)/)



    ;;;局地海温
    ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))

    ssta_tna_mam = dim_avg_n_Wrap(ssta_mam(:,{20:30},{270:300}), (/1,2/))

    ssta_ind_mam = dim_avg_n_Wrap(ssta_mam(:,{10:30},{50:100}), (/1,2/))
    ;;;

    printMinMax(MSE_nep_mam,1)
    

    corr1 = escorc_n(MSE_nep_mam, ssta_mam_nep, 0, 0)
    corr2 = escorc_n(MSE_nep_mam, ssta_tna_mam, 0, 0)
    corr3 = escorc_n(MSE_nep_mam, ssta_ind_mam, 0, 0)

   
    plot = new(4,graphic)
;;;;
    MSE_nep_mam = where(ismissing(MSE_nep_mam), ssta_mam_nep@_FillValue, MSE_nep_mam)
    MSE_nep_mam@_FillValue = ssta_mam_nep@_FillValue
    
; wks   = gsn_open_wks("pdf", "Nodetrend-MSE与SST的时间序列-MAM-"+runmean_time+"年滑动平均")
  wks   = gsn_open_wks("pdf", "标准化的MSE与SST的时间序列-MAM-"+runmean_time+"年滑动平均-1980-2020")
; wks   = gsn_open_wks("pdf", "MSE与SST的相关-MAM-去掉局地的影响")
;;**************************************************************
;画图
;;***************************************************************
      
    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "MSE(80W-150W,0-20N) & SST decadal change "   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02
    res@tiYAxisString = "MSE anomaly"           ; y axis title
    res@tiXAxisString = "Year"

    
    res@xyLineColors      = (/"black","black","red","blue"/)  ; line colors
   

    res@xyDashPatterns    = (/0.0,1.0,1.0,1.0/)        ; line patterns
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    res@gsnYRefLine = 0.0

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = year2
    res@trXMinF = year1
    res@trYMinF = -1.0
    res@trYMaxF = 1.0

    res@tmXTOn = False
    res@tmYROn = False
   

    res@tmXBValues=ispan(1950,2020,20)
    res@tmXBLabels=ispan(1950,2020,20)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(1950,2020,5)
    res@tmYRMode="Automatic"


    plot(0)  = gsn_csm_xy(wks, ispan(year1, year2, 1), (/MSE_nep_mam,ssta_mam_nep,ssta_tna_mam,ssta_ind_mam/), res)

    res@xyLineColors      = (/"black","black","red","blue"/)  ; line colors
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/1.0/)        ; line thicknesses
    
    res@xyDashPatterns    = (/0.0,1.0,1.0,1.0/)        ; line patterns

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

    textres@lgLabels = (/"MSE","SST(local) "+sprintf("%6.2f", corr1),"SST(TNA) "+sprintf("%6.2f", corr2),"SST(IND)"+sprintf("%6.2f", corr3)/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


  

    draw(plot)
    frame(wks)

end



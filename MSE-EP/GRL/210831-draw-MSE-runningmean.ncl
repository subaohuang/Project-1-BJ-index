;;;;for each month
undef("rmstd") ;;;滑动std
function rmstd(x,nave)
begin
dim = dimsizes(x)
length=dim(0)
std = new((/length-nave+1/),typeof(x))

;;;
do t=nave,length

   ; m = x((t-nave)*12:t*12-1)
   m = x(t-nave:t-1)
    ;;;rm annual cycle in each running period followed as kim and cai,2013 ncc
   ; m = (/rmMonAnnCycTLL(m)/)
    ;;; detrend 
   ; m = (/dtrend_msg_n(ispan(1,nave*12,1),m,False,True,0)/)
    ;;; 7-year running mean 
   ; m = (/ m - runave_n(m, 7*12,0,0)/)
   ; std_mv = dim_avg_n(m,0)
   std_mv = dim_stddev_n_Wrap(m,0)

   std(t-nave) = std_mv
   delete(m)
end do 
return std
end 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 


begin
    ;;;time span
    startmon =198001
    endmon = 202012

    running_time = 15

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1


    ;; 
    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_ncep1_50-20.nc", "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    MSE = f1->MSE(timeind,:,:,:)

    ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)

    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    ; MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

    ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)

    ; MSE_col = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_col = (/MSE_col/10^6/)


    MSE_nep = dim_avg_n_Wrap(MSE_col(:,{0:20},{260:280}), (/1,2/))

    MSE_nep_12 = split_12(MSE_nep)

  
    MSE_nep_mam = month_to_annual(MSE_nep, 1)
    
    MSE_nep_mam = (/dim_avg_n(MSE_nep_12(0:4,:),0)/)


;;;;;11 windows std 

      ; std_nino1_2 = rmstd(MSE_nep_mam,running_time)

      std_nino1_2 = runave_n_Wrap(MSE_nep_mam, 15, 0, 0)

      
;;;


  ; wks   = gsn_open_wks("pdf", "JFMAM MSE指数的滑动标准差"+""+running_time)
  wks   = gsn_open_wks("pdf", "JFMAM MSE指数的滑动平均"+""+running_time)

;;**************************************************************
;画图
;;***************************************************************
      
    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    ; res@tiMainString  = "15-years window running std for JFMAM MSE "   ; title
    res@tiMainString  = "15-years window running mean for JFMAM MSE "   ; title

    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02


    ; res@tiYAxisString = "Std"           ; y axis title
    res@tiYAxisString = "MSE anomaly"           ; y axis title
    res@tiXAxisString = "Year"

    
    res@xyLineColors      = (/"black","green","blue","red"/)  ; line colors
   

    res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    res@gsnYRefLine = 0.0

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = endyear
    res@trXMinF = startyear
    ; res@trYMinF = 0.0
    ; res@trYMaxF = 2.0

    res@trYMinF = -0.5
    res@trYMaxF = 0.5

    res@tmXTOn = False
    res@tmYROn = False

     
    res@vpHeightF = 0.4
    res@vpWidthF = 0.7

   

    res@tmXBValues=ispan(startyear,endyear,10)
    res@tmXBLabels=ispan(startyear,endyear,10)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(startyear,endyear,2)
    res@tmYRMode="Automatic"


    ; plot = gsn_csm_xy(wks, ispan(startyear, endyear-(running_time-1), 1), (/std_nino1_2/), res)
    plot = gsn_csm_xy(wks, ispan(startyear, endyear, 1), (/std_nino1_2/), res)

  
    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = res@xyLineThicknessF                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"MSE"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    draw(plot)
    frame(wks)

end


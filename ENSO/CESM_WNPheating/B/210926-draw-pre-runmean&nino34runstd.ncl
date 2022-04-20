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

 ; ;;;读入HadiSST海温数据
      diri="/home/yangsong3/data-observation/linshh/data/sst/"
      fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
      f1=addfile(fils1, "r")
      date:=cd_calendar(f1->time, 1)
      
      timeind:=ind(date.le.endmon.and.date.ge.startmon)

      sst =  lonFlip(f1->sst(timeind,:,:))
   
      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)

      ssta_de = ssta
      ; ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
      ssta_index = dim_avg_n_Wrap(ssta_de(:,{5:-5},{210:270}), (/1,2/))
      printVarSummary(ssta_index)

      nino3_4 = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
      nino3_4(1:) = nino3_4(0:dimsizes(nino3_4)-2)
      nino3_4(dimsizes(nino3_4)-1) = ssta_index(dimsizes(nino3_4)-1)
      printVarSummary(nino3_4)

      ; nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)

;;;;;11 windows std 
      std_nino1_2 = rmstd(nino3_4, 11)
;;;;


  ;;;
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    f2=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p=f2->precip(timeind,:,:)

    p = (/rmMonAnnCycTLL(p)/)

    p_wnp = dim_avg_n_Wrap(p(:,{0:20},{125:160}), (/1,2/))

    ; p_wnp = (/dim_standardize_n(p_wnp, 1, 0)/)

    p_wnp_MAM = month_to_season(p_wnp, "MAM")
    
    p_wnp_MAM = (/runave(p_wnp_MAM, 11, 0)/)

printVarSummary(std_nino1_2)
printVarSummary(p_wnp_MAM)


  wks   = gsn_open_wks("pdf", "降水滑动平均以及Nino滑动标准差"+""+running_time)

;;**************************************************************
;画图
;;***************************************************************
      
    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    ; res@tiMainString  = "15-years window running std for JFMAM MSE "   ; title
    res@tiMainString  = "11-years window running std for nino3_4 "   ; title

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

    ; res@trYMinF = -0.5
    ; res@trYMaxF = 0.5
    
       res@trYMinF = -2.0
    res@trYMaxF = 2.0


    res@tmXTOn = False
    res@tmYROn = False

     
    res@vpHeightF = 0.4
    res@vpWidthF = 0.7


    res@tmXBValues=ispan(startyear,endyear,10)
    res@tmXBLabels=ispan(startyear,endyear,10)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(startyear,endyear,2)
    res@tmYRMode="Automatic"



    plot = gsn_csm_xy(wks, ispan(startyear, endyear, 1), (/p_wnp_MAM/), res)
   ; plot = gsn_csm_xy(wks,  ispan(startyear+5, endyear-5, 1), std_nino1_2, res)   

    ployres = True
    ployres@gsLineColor = "red"
    ; ployres@gsLineThicknessF
    dim = gsn_add_polyline(wks, plot, ispan(startyear+5, endyear-5, 1), std_nino1_2, ployres)
    
    ; ;;添加legend 
    ; gres = True
    ; gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;       ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ; gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ; lineres = True
    ; ; lineres@lgLineColors = (/"black"，/) ; line colors
    ; lineres@lgLineColors       = res@xyLineColors ; colors chosen
    ; lineres@lgLineThicknesses = res@xyLineThicknessF                       ; line thicknesses
    ; lineres@LineLengthPercent = 9.
    ; lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    ; textres = True

    ; textres@lgLabels = (/"MSE"/)
    ; textres@lgLabelFontHeightF = 0.01
    ; ; textres@lgLabelOpacity = (/0/) 
   
    ; plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    draw(plot)
    frame(wks)

end


;;;;for each month
undef("rmstd") ;;;滑动std
function rmstd(x,nave)
begin
dim = dimsizes(x)
; length=dim(0)
length=dim(0)/12
std = new((/length-nave+1/),typeof(x))

;;;
do t=nave,length
   m = x((t-nave)*12:t*12-1)
   ; m = x(t-nave:t-1)
    ;;;rm annual cycle in each running period followed as kim and cai,2013 ncc
   ; m = (/rmMonAnnCycTLL(m)/)
    ;;; detrend 
   ; m = (/dtrend_msg_n(ispan(1,nave*12,1),m,False,True,0)/)
    ;;; 7-year running mean 
   ; m = (/ m - runave_n(m, 7*12,0,0)/)
   std_mv = dim_stddev_n_Wrap(m,0)
   std(t-nave) = std_mv
   delete(m)
end do 
return std
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
;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
   ;;
;;;

;;;
    ; ssta_mam = month_to_season(ssta, "DJF")
    ssta_mam = ssta

;;;DJF nino index 
     ; ;;局地海温
     ;  ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))

     ;;nino1+2:
      nino1_2 = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{270:280}), (/1,2/))
    ;;;
    
     ;;nino3:
      nino3 = dim_avg_n_Wrap(ssta_mam(:,{5:-5},{210:270}), (/1,2/))
    ;;;
     
     ;;;nino4 
       nino4 = dim_avg_n_Wrap(ssta_mam(:,{5:-5},{160:210}), (/1,2/))
    ;;;

     ;;;nino 3+4 
       nino3_4 = dim_avg_n_Wrap(ssta_mam(:,{5:-5},{190:240}), (/1,2/))
    ;;;
;;;;


;;;;;11 windows std 

      std_nino1_2 = rmstd(nino1_2,running_time)
      std_nino3= rmstd(nino3,running_time)
      std_nino4= rmstd(nino4,running_time)
      std_nino3_4 = rmstd(nino3_4,running_time)
;;;


  wks   = gsn_open_wks("pdf", "80-20-Nino指数的15滑动标准差")

;;**************************************************************
;画图
;;***************************************************************
      
    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "15-years window running std for El Nino index "   ; title

    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02


    res@tiYAxisString = "Std"           ; y axis title
    res@tiXAxisString = "Running window mid-year"

    
    res@xyLineColors      = (/"black","green","blue","red"/)  ; line colors
   

    res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    res@gsnYRefLine = 0.0

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    ; res@trXMaxF = endyear-(running_time-1)
    res@trXMaxF = endyear
    res@trXMinF = startyear
    res@trYMinF = 0.5
    res@trYMaxF = 1.4

    res@tmXTOn = False
    res@tmYROn = False
   

    res@tmXBValues=ispan(startyear,endyear,10)
    res@tmXBLabels=ispan(startyear,endyear,10)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(startyear,endyear,2)
    res@tmYRMode="Automatic"
    
    res@vpHeightF = 0.4
    res@vpWidthF = 0.7

    plot = gsn_csm_xy(wks, ispan(startyear+(running_time-1)/2, endyear-(running_time-1)/2, 1), (/std_nino1_2,std_nino3,std_nino4,std_nino3_4/), res)

  
    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 70    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = res@xyLineThicknessF                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"Nino1.2","Nino3","Nino4","Nino3.4"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    draw(plot)
    frame(wks)

end


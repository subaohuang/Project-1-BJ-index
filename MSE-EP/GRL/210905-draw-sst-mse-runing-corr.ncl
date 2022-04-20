;;;;for each month
undef("rmcorr") ;;;滑动std
function rmcorr(x,y,nave)
begin
dim = dimsizes(x)
; length=dim(0)
length=dim

std = new((/length-nave+1/),typeof(x))
;;;
do t=nave,length
   m = nave
   std_mv =escorc(x(t-m:t-1),y(t-m:t-1)) 
   std(t-nave) = std_mv
   delete(m)
end do 
return std
end 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
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
    ssta = sst
    ssta=(/rmMonAnnCycTLL(sst)/)


    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    copy_VarCoords(sst, ssta)
   ;;
;;;


;;;;;; MSE
 ; f1 = addfile("MSE_era5.nc", "r")
 data = "ncep1"
    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_"+data+"_50-20.nc", "r")

    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    startmon = 198001
    endmon = 202012

    startyear = 1980
    endyear = 2020

    timeind :=ind(date.le.endmon.and.date.ge.startmon)

    MSE = f1->MSE(timeind,:,:,:)
  ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)
    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)
    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)


    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

    ; MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_nep = MSE_col
    MSE_nep = (/MSE_col/10^6/)

;;;;;;split to 12 month  is a array size (12,yearnum)
  
     MSE_nep_12 = split_12(MSE_nep)
     ssta_12 = split_12(ssta)
   
;;;;;


;;;等到年平均的数据
    MSE_nep_annual = month_to_annual(MSE_nep, 1)
    MSE_nep_annual!0 = "time"
    MSE_nep_annual&time = ispan(startyear,endyear,1)

    ssta_annual = month_to_annual(ssta, 1)
    ssta_annual!0 = "time"
    ssta_annual&time = ispan(startyear,endyear,1)
    
    ssta_annual = (/dim_avg_n(ssta_12(0:4,:,:,:),0)/)
    MSE_nep_annual = (/dim_avg_n(MSE_nep_12(0:4,:,:,:),0)/)
;;;;;

;;;
    ; ssta_mam = month_to_season(ssta, "DJF")
    ssta_mam = ssta_annual

;;;DJF nino index 
     ; ;;局地海温
     ;  ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))

     ;;nino1+2:
      indian = dim_avg_n_Wrap(ssta_mam(:,{20:30},{300:330}), (/1,2/))
    ;;;
    
     ;;nino3:
      local_sst = dim_avg_n_Wrap(ssta_mam(:,{0:20},{250:280}), (/1,2/))
    ;;;
     
     ;;;nino4 
       atlantic = dim_avg_n_Wrap(ssta_mam(:,{20:-20},{40:110}), (/1,2/))
    ;;;
printVarSummary(MSE_nep_annual)
;;;;
    MSE_nep_annual_mean = dim_avg_n_Wrap(MSE_nep_annual(:,{0:20},{250:290}),(/1,2/))
   
;;;;;11 windows std 



      std_nino1_2 = rmcorr(indian,MSE_nep_annual_mean, running_time)
      std_nino3= rmcorr(local_sst,MSE_nep_annual_mean, running_time)
      std_nino4= rmcorr(atlantic,MSE_nep_annual_mean, running_time)
      
;;;


  wks   = gsn_open_wks("pdf", "MSE-SST-corr-running-window")

;;**************************************************************
;画图
;;***************************************************************
      
    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "15-years window running corr "   ; title

    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02


    res@tiYAxisString = "Correlation "           ; y axis title
    res@tiXAxisString = "Running window mid-year"

    
    res@xyLineColors      = (/"black","green","blue"/)  ; line colors
   

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
    res@trYMinF = 0.0
    res@trYMaxF = 1.0

    res@tmXTOn = False
    res@tmYROn = False
   

    res@tmXBValues=ispan(startyear,endyear,10)
    res@tmXBLabels=ispan(startyear,endyear,10)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(startyear,endyear,2)
    res@tmYRMode="Automatic"
    
    res@vpHeightF = 0.4
    res@vpWidthF = 0.7

    plot = gsn_csm_xy(wks, ispan(startyear+(running_time-1)/2, endyear-(running_time-1)/2, 1), (/std_nino1_2,std_nino3,std_nino4/), res)

  
    ;;添加legend 
    gres = True
    gres@YPosPercent = 15.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = res@xyLineThicknessF                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"Indian","Local","North Atlantic"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    draw(plot)
    frame(wks)

end



;;**************************************************************
;;;main program
;;***************************************************************
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


begin
;读入数据
model_name = (/ "ACCESS1-0","ACCESS1-3", "CanESM2","CMCC-CESM","CMCC-CM","CMCC-CMS",\
            "CSIRO-Mk3-6-0", "FGOALS-s2","GFDL-CM2p1","GFDL-CM3" ,"GFDL-ESM2G","GFDL-ESM2M",\
            "GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadCM3","HadGEM2-CC","HadGEM2-ES",\
            "inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR",\
            "MPI-ESM-P","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
; model_name = (/"GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC"/)
; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(dim_model,graphic)

; std_all = new((/dim_model,/))
wks = gsn_open_wks("pdf","30个CESM5模式的-STD-HadISST")         ; send graphics to PNG file
do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start read"+" "+model+"**********************************************")
    yearstart = 190001
    yearend = 200412


    fils1 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")
    ; print(fils1)
    f1 := addfiles(fils1, "r")

    time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)

    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))

    sst_r := f1[:]->tos(timeind,:,:)
    sst := lonFlip(sst_r)

    ssta :=(/rmMonAnnCycTLL(sst)/)
   
    copy_VarCoords(sst, ssta)
    printVarSummary(ssta)
    ;;;
;;;;

;;;;计算标准偏差(5S-5N)

        ssta_mean = dim_avg_n_Wrap(ssta(:,{-5:5},{120:280}), 1)

        std1 = dim_stddev_n(ssta_mean(0:49*12-1,:), 0)

        std2 = dim_stddev_n(ssta_mean(59*12:,:), 0)

        x = ssta_mean&lon
        delete(ssta_mean)
        printVarSummary(std1)
   
;;;;

        data_all = (/std1,std2/)
        delete(std1)
        delete(std2)

;;;;draw 
;;*************************************************
; create plot
;**************************************************
 

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = model   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisString = "SST anomaly s.d."           ; y axis title
    res@tiXAxisString = ""
    

    res@tmXBLabelFontHeightF=0.01
    res@tmYLLabelFontHeightF=0.01
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0


    res@xyLineColors      = (/"blue","red"/)  ; line colors
    res@xyLineThicknesses = (/4.0,4.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0/)        ; line patterns

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 280.0
    res@trXMinF = 120.0
    res@tmXTOn = False
    res@tmYROn  = False

    ; res@tmXBValues=(/0,90,180,270,357.5/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(0,360,10)
    res@vpWidthF = 0.6
    res@vpHeightF = 0.4

    res@tmXBValues=ispan(135,270,45)
    res@tmXBLabels=(/"135~S~o~N~E","180~S~o","135~S~o~N~W","90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(120,280,15)
    res@tmYRMode="Automatic"

     plot(k) = gsn_xy(wks,x,data_all,res)          ; Draw an XY plot with 1 curve.
end do 

pres= True
gsn_panel(wks, plot, (/6,5/), pres)

end

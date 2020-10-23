
;;**************************************************************
;;;main program
;;***************************************************************

function  ploy(std1,std2)
begin

dim  = dimsizes(std1)
mnmx  = new ((/6,dim(1)/), typeof(std1))

mnmx(0,0:) =  dim_avg_n_Wrap(std1, 0)
mnmx(1,0:) = dim_avg_n_Wrap(std2, 0)

;;;计算标准偏差 
sigma1  = dim_stddev_n(std1, 0)
sigma2  = dim_stddev_n(std2, 0)

; mnmx(2,0:) = dim_max_n_Wrap(std1, 0)
; mnmx(3,0:) = dim_min_n_Wrap(std1, 0)

; mnmx(4,0:) = dim_max_n_Wrap(std2, 0)
; mnmx(5,0:) = dim_min_n_Wrap(std2, 0)

mnmx(2,0:) = mnmx(0,:) + sigma1
mnmx(3,0:) = mnmx(0,:) - sigma1

mnmx(4,0:) = mnmx(1,:) + sigma2
mnmx(5,0:) = mnmx(1,:) - sigma2
 
return mnmx
end 

begin
;读入数据

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;     "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;     "GFDL-CM4", "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         "HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM"/)
model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-ESM1",\
              "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
   , "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
  "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
        ,"HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
        "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM","SAM0-UNICON","UKESM1-0-LL"/)

; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(1,graphic)

std_all_1 = new((/dim_model,64/),float)
std_all_2 = new((/dim_model,64/),float)
wks = gsn_open_wks("pdf","MME的CESM模式的-STD-HadISST")         ; send graphics to PNG file
do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start read"+" "+model+"**********************************************")
    yearstart = 190001
    yearend = 201412


    fils1 := systemfunc("ls /home/yangsong3/data-observation/cmip6/"+model+"/historical/tos_fixed/tos_Omon_"+"*")
    ; print(fils1)

   if dimsizes(fils1).eq.1 then 

        f1 := addfile(fils1, "r")

        time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)

        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        delete(date)     
        sst = lonFlip(f1->tos(timeind,:,:))
    else 
        f1 := addfiles(fils1, "r")

        time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)
        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        delete(date)     
        sst = lonFlip(f1[:]->tos(timeind,:,:))
    end if

        ssta :=(/rmMonAnnCycTLL(sst)/)
        copy_VarCoords(sst, ssta)
        printVarSummary(ssta)
    ;;;
;;;;

;;;;计算标准偏差(5S-5N)

        ssta_mean = dim_avg_n_Wrap(ssta(:,{-5:5},{120:280}), 1)

        std1 = dim_stddev_n(ssta_mean(0:49,:), 0)

        std2 = dim_stddev_n(ssta_mean(59:,:), 0)

        x = ssta_mean&lon
        delete(ssta_mean)


        std_all_1(k,:) = std1
        std_all_2(k,:) = std2
;;;;
end do 

    mnmx = ploy(std_all_1,std_all_2)



;;;;draw 
;;*************************************************
; create plot
;**************************************************
 

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "Standard deviations of SSTA in CMIP6 Historical run"   ; title
    res@tiMainFontHeightF = 0.015
    res@tiYAxisString = "SST anomaly s.d."           ; y axis title
    res@tiXAxisString = ""
    

    ; res@tmXBLabelFontHeightF=0.02
    ; res@tmYLLabelFontHeightF=0.02
    ; res@tmXBMajorThicknessF=1.0
    ; res@tmYLMajorThicknessF=1.0
    ; res@tmBorderThicknessF=1.0


    


    res@xyLineColors      = (/"blue","red"/)  ; line colors
    res@xyLineThicknesses = (/4.0,4.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0/)        ; line patterns

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 280.0
    res@trXMinF = 120.0
    res@trYMinF = 0.0
    res@trYMaxF = 1.8

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
    

    mid_plot = plot
    plot_std1 = new(1,graphic)
    plot_std2 = new(1,graphic)
    mid_plot= gsn_csm_xy(wks,x,mnmx(0:1,:),res)       ; create line plot

    
    
    delete(res@xyLineColors)
    res@gsnXYFillColors = "blue"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    res@gsnXYFillOpacities=0.1

    plot_std1= gsn_csm_xy(wks,x,mnmx(2:3,:),res)       ; create line plot
  
    

    res@gsnXYFillColors = "pink"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    ; res@gsnXYFillOpacities=0.1

    plot_std2= gsn_csm_xy(wks,x,mnmx(4:5,:),res)       ; create line plot
    
    overlay(mid_plot, plot_std1)
    overlay(mid_plot, plot_std2)

    
    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"Blue","red"/) ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"1900-1950","1960-2014"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 


    mid = simple_legend(wks,mid_plot,gres,lineres,textres)   
   




draw(mid_plot)
frame(wks)
end

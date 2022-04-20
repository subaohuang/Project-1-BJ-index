
;;**************************************************************
;;;main program
;;***************************************************************

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

std_all_1 = new((/dim_model,72,144/),float)
std_all_2 = new((/dim_model,72,144/),float)
wks = gsn_open_wks("pdf","201026-MME的CESM模式的STD-水平分布")         ; send graphics to PNG file
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

    ; ssta_mean = dim_avg_n_Wrap(ssta(:,{-5:5},{120:280}), 1)

        std1 = dim_stddev_n_Wrap(ssta(0:49,:,:), 0)

        std2 = dim_stddev_n_Wrap(ssta(59:,:,:), 0)
        
        printVarSummary(std1)
        
        std_all_1(k,:,:) = std1
        std_all_2(k,:,:) = std2
;;;;
end do 
   
   std_mme1 = dim_avg_n_Wrap(std_all_1, 0)
   std_mme2 = dim_avg_n_Wrap(std_all_2, 0)

   copy_VarCoords(std1, std_mme1)
   copy_VarCoords(std1, std_mme2)

   

   std_change = std_mme1
   std_change = (/std_mme2 - std_mme1/)



;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =True

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0



    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.06
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 110.0
    res@mpMaxLonF             = 300.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
    ;;


    ; ;;;LineArrow
    ; res@gsnScalarContour        = True

    ; res@gsnLeftString = "AGCM JJA alt warming induced 200hpa geopotensial Height and U,V diff (SEN - CTRL)"


    ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; ; res@vcGlyphStyle="CurlyVector"

    ; res@vcMinDistanceF           = 0.02           ; thin vectors
    ; res@vcMinMagnitudeF          =0.5
    ; ; res@vcMaxMagnitudeF = 7.0
    ; res@vcRefMagnitudeF          = 5.0
    ; res@vcRefLengthF             = 0.05     ; define length of v
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ; res@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
    ; res@vcVectorDrawOrder="PostDraw"
    res@gsnRightString = ""
    ; plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)
    plot = gsn_csm_contour_map(wks, std_change,res)

draw(plot)
frame(wks)
end 


 


  
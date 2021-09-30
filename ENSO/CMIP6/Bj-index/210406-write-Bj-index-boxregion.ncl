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

   ;  model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-CSM2-MR",\
   ;            "CAMS-CSM1-0","CanESM5", "CAS-ESM2-0" ,"CESM2","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","CNRM-CM6-1"\
   ; , "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
   ;   "GISS-E2-1-G","HadGEM3-GC31-LL","IPSL-CM6A-LR","MIROC-ES2L","MIROC6", \
   ;  "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3",  "NorESM2-LM"/)

    model_name = (/"CanESM5"/)
    dim_model = dimsizes(model_name)

    plot = new(dim_model,graphic)
    iplot= plot
    plot1 = plot

    yearstart = 190001
    yearend = 201412

    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)

    
    startyear_r = str_split_by_length(yearstart, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(yearend, 4 )

    yearnum = stringtoint(endyear_r(0)) - stringtoint(startyear_r(0)) + 1
    print(yearnum)
;;;


    wks = gsn_open_wks("pdf", "cmip6模式 Heat content pattern EOF")

    lat1 = -20
    lat2 = 20
    lon1 = 120
    lon2 = 280


do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon/thetao_Omon_"+model+"_historical_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    tos_or = lonFlip(f2[:]->thetao(timeind,{0:300},{lat1:lat2},:))

    depth = tofloat(tos_or&lev)
    dz = depth
    dz(1:)= depth(1:) - depth(0:dimsizes(depth)-2)
    dz(0)= depth(0) 
;;;
    
    tos =dim_avg_wgt_n_Wrap(tos_or, dz, 1, 1)  

    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos = (/rmMonAnnCycTLL(tos)/)

    
;;;输入数据

    diri = "/home/ys17-19/lsh/Project/ENSO/CMIP6/Bj-index/EOF/"

    f1 = addfile(diri+"EOF-ssta-"+model+".nc", "r")
    ts_model = f1->ts_model
    pattern_model = f1->pattern_model
    pcvar = f1->pcvar
;;;

    ts_sd = dim_standardize(ts_model, 1)

    
    tos_enino  =  regCoef_n(ts_sd, tos, 0, 0)

    copy_VarCoords(tos(1,:,:), tos_enino)
    printVarSummary(tos_enino)



;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn = True


    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0



    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -0.5
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.05
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -20.0
    res@mpMaxLatF             = 20.0
    res@mpMinLonF             = 125.0
    res@mpMaxLonF             = 280.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    res@gsnRightString = pcvar +""
    res@gsnLeftString = model 

    res@gsnLeftStringFontHeightF = 0.03
       
    res@vpHeightF = 0.3
    res@vpWidthF = 0.7
    res@mpShapeMode = "FreeAspect"
    plot_model = gsn_csm_contour_map(wks,tos_enino,res)
    plot(k) = plot_model
    isolines=get_isolines(plot_model@contour, 0.0)

    y=isolines(0,:)
    x=isolines(1,:)

    location = x(ind(y.eq.0))
    print(location)  ;;;145.2 E
 

 print("********************************finish calculate"+" "+model+"************************************")
end do
   
    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    pres@pmLabelBarHeightF = 0.05
    pres@pmLabelBarWidthF  = 0.6
    pres@pmLabelBarOrthogonalPosF = -0.02
    ; pres@lbOrientation = "Vertical"
    pres@lbLabelFontHeightF = 0.01
    gsn_panel(wks, plot, (/6,5/), pres)
end 


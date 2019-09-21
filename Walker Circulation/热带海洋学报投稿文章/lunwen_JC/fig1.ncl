
;;****************************************************************
;;main program
;;*************************************************************

;;测试
begin
    yearstart=197901 
    yearend=200812
    year1=1979
    year2=2008
  
    model="amip"
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    fils1=systemfunc("ls "+diri+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)


    timeind=ind(date.le.yearend.and.date.ge.yearstart)
   
    msf_time=f1->msf_time(:,{80:290},timeind)
    msf_time=(/msf_time/10^9/)
    msf_time&level@units="hpa"
    msf_time&time=ispan(0, dimsizes(timeind)-1, 1)

    delete(msf_time&time)
    msf_time&time=date(timeind)
     delete(date)
    printVarSummary(msf_time)

    msf_clm_amip=dim_avg_n_Wrap(msf_time, 2)
    delete(msf_time)
    delete(timeind)
    

    ;;;;读入CMIP集合平均的数据

    model="cmip"
    diri2="/mnt/f/share/CMIP5/msf_history/"
    fils2=systemfunc("ls "+diri2+"msf_time_"+model+"_cmip_global.nc")
    f2=addfile(fils2, "r")
    date=cd_calendar(f2->time, 1)

    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    msf_time=f2->msf_time(:,{80:290},timeind)
    msf_time=(/msf_time/10^9/)
    msf_time&level@units="hpa"
    msf_time&time=ispan(0, dimsizes(timeind)-1, 1)

    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    printVarSummary(msf_time)

    msf_clm_cmip=dim_avg_n_Wrap(msf_time, 2)
    msf_clm_cmip_inten=avg(dim_avg_n(msf_clm_cmip(:,{150:210}),0))
    delete(msf_time)
    delete(timeind)
    
    ;;;读入再分析
  
    f3=addfile("/mnt/f/share/rea_msf/msf_time_ens_global.nc", "r")
    date=cd_calendar(f3->time, 1)

    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    msf_time=f3->msf_time(:,{80:290},timeind)
    msf_time=(/msf_time/10^9/)
    msf_time&level@units="hpa"
    msf_time&time=ispan(0, dimsizes(timeind)-1, 1)

    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    printVarSummary(msf_time)
    msf_clm_era=dim_avg_n_Wrap(msf_time, 2)
    msf_clm_era_inten=avg(dim_avg_n(msf_clm_era(:,{150:210}),0))
    
    delete(msf_time)
    delete(timeind)

    ;;;相对大小
    inten_r=(msf_clm_cmip_inten - msf_clm_era_inten)/msf_clm_era_inten
    print("cmip气候态强度相对于观测的大小inten_r="+inten_r)

 ;;;

 ;;; 读入观测ERSST与HadIsst海温
   print("*****************读入海温场*********************")
   ;;;海温场
    fils4=systemfunc("ls /mnt/f/share/ersst/sst.mnmean.v4.nc")
    f4=addfile(fils4, "r")
    date=cd_calendar(f4->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_mre=f4->sst(timeind,{-20:20},{80:290})

    sst_mre_clm=dim_avg_n_Wrap(sst_mre, 0)
    delete(timeind)
    delete(date)
   

    print("*****************读入海温场*********************")
    ;;;海温场
    fils5=systemfunc("ls /mnt/f/share/ersst/HadIsst_sst.nc")
    f5=addfile(fils5, "r")
    date=cd_calendar(f5->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_had=lonFlip(f5->sst(timeind,{-20:20},:))
    sst_had_clm=dim_avg_n_Wrap(sst_had(:,:,{80:290}), 0) 
    delete(timeind)
    delete(date)



    print("*****************读入海温场*********************")
    ;;;海温场
    fils6=systemfunc("ls /mnt/f/share/CMIP5/tos_chazhi/sst_cmip_r1i1p1_197901-201212.nc")
    f6=addfile(fils6, "r")
    date=cd_calendar(f6->time, 1)


    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_cmip=f6->sst(timeind,{-20:20},{80:290})
    sst_cmip_clm=dim_avg_n_Wrap(sst_cmip, 0) 
    delete(timeind)
    delete(date)
 
  
   wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/AMIP,CMIP,再分析的MSF以及SST的气候态比较")
;********************************************************
;;draw
;********************************************************

    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
    plot=new(6,graphic)
    plot_sst=new(3,graphic)

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftStringFontHeightF=0.03
    res2@gsnRightStringFontHeightF=0.03
    res2@gsnLeftString="Reanalysis"
    res2@gsnRightString="Mean state"
    res2@cnLinesOn=True
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -6.0
    res2@cnMaxLevelValF = 6.0
    res2@cnLevelSpacingF = 1.0
    res2@cnFillPalette="BlueDarkRed18"
    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    res2@lbLabelBarOn=False
    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False

    res2@tmXBLabelFontHeightF=0.025
    res2@tmYLLabelFontHeightF=0.025
    res2@tmYLMajorThicknessF=2.0
    res2@tmXBMajorThicknessF=2.0
    res2@tmBorderThicknessF=4.0

    res2@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res2@tmXBValues=ispan(90,270,30)
    res2@tmXBLabels=(/"90E","120E","150E","180","150W","120W","90W"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(80,290,10)
   
    res2@tmYRMode="Automatic"

    res2@vpHeightF=0.4
    res2@vpWidthF=0.6

    res2@cnLineThicknessF = 2.
    res2@cnFillOn = True
    res2@gsnContourNegLineDashPattern=1
    res2@gsnContourZeroLineThicknessF=4.0

     
    plot(0) = gsn_csm_pres_hgt(wks, msf_clm_era , res2)
    res2@gsnLeftString="AMIP"
    res2@gsnLeftStringOrthogonalPosF=0.05
    res2@gsnRightStringOrthogonalPosF=0.05

    res2@tmYLOn=False
    ; res2@lbLabelBarOn=True
    ; res2@lbLabelFontHeightF=0.03
    ; res2@lbLeftMarginF=-4.5
    ; res2@lbRightMarginF=-4.5
    ; res2@pmLabelBarOrthogonalPosF=-0.01
    ; res2@pmLabelBarHeightF=0.14

    plot(1) = gsn_csm_pres_hgt(wks, msf_clm_amip , res2)
     res2@gsnLeftString="CMIP5"
     res2@tmYLOn=False
     res2@lbLabelBarOn=True
     res2@lbOrientation="Vertical"
    plot(2) = gsn_csm_pres_hgt(wks, msf_clm_cmip , res2)

   
    res_sst = True
    res_sst@gsnAddCyclic=False
    res_sst@gsnDraw=False
    res_sst@gsnFrame=False
    res_sst@gsnLeftString="ERSST"
    res_sst@gsnRightString="Mean state"
    res_sst@gsnLeftStringFontHeightF=0.03
    res_sst@gsnRightStringFontHeightF=0.03

    res_sst@cnLinesOn=True
    res_sst@cnLevelSelectionMode ="ManualLevels"
    res_sst@cnMinLevelValF = 20.0
    res_sst@cnMaxLevelValF = 29.0
    res_sst@cnLevelSpacingF = 1.0

    res_sst@cnFillPalette="BlueDarkRed18"
    res_sst@cnLineLabelFontHeightF=0.01

    res_sst@cnFillOn = True 
    res_sst@cnLineLabelsOn=False
    res_sst@cnInfoLabelOn=False
    res_sst@lbLabelBarOn=False

    res_sst@tmXTOn=False
    res_sst@tmYROn=False
    res_sst@tmXBLabelFontHeightF=0.02
    res_sst@tmYLLabelFontHeightF=0.02
    res_sst@tmYLMajorThicknessF=2.0
    res_sst@tmXBMajorThicknessF=2.0
    res_sst@tmBorderThicknessF=4.0
    
    res_sst@mpFillDrawOrder="PostDraw"
    res_sst@mpLandFillColor="gray"
    res_sst@mpFillOn=True

    res_sst@vpHeightF=0.4
    res_sst@vpWidthF=0.6


    res_sst@mpLimitMode           = "LatLon"
    res_sst@mpMinLatF             = -20
    res_sst@mpMaxLatF             = 20
    res_sst@mpMinLonF             = 80
    res_sst@mpMaxLonF             = 290
    res_sst@mpCenterLonF          = (80+290)/2.0
    res_sst@mpShapeMode ="FreeAspect"

    

    ; plot_sst(0) = gsn_csm_contour_map(wks, sst_mre_clm , res_sst)
    ; res_sst@gsnLeftString="HadISST"
    ; res_sst@tmYLOn=False
    ; plot_sst(1) = gsn_csm_contour_map(wks, sst_had_clm , res_sst)
    ; res_sst@gsnLeftString="CMIP5"
    ; res_sst@tmYLOn=False
    ; plot_sst(2) = gsn_csm_contour_map(wks, sst_cmip_clm , res_sst)

    
    plot(3) = gsn_csm_contour_map(wks, sst_mre_clm , res_sst)
    res_sst@gsnLeftString="HadISST"
    res_sst@tmYLOn=False
    ; res_sst@lbLabelBarOn=True
    ; res_sst@lbLeftMarginF=-4.5
    ; res_sst@lbRightMarginF=-4.5
    ; res_sst@pmLabelBarOrthogonalPosF=0.15
    ; res_sst@pmLabelBarHeightF=0.14
     
    plot(4) = gsn_csm_contour_map(wks, sst_had_clm , res_sst)
    res_sst@gsnLeftString="CMIP5"
    res_sst@tmYLOn=False
    ; res_sst@lbLabelBarOn=False
    res_sst@lbLabelBarOn=True
     res_sst@lbOrientation="Vertical"
     ; res_sst@pmLabelBarOrthogonalPosF=-0.14
    plot(5) = gsn_csm_contour_map(wks, sst_cmip_clm , res_sst)




        ;;;添加地形
        ployres                  = True          ; polyline mods desired
        ployres@gsLineThicknessF = 8.0            ; line thickness   
        ployres@gsLineColor      = "black"         ; color of lines
        ployres@gsLineDashPattern=1
        ;---Attach the polyline
        ploylevel=1000
        ploylon1=fspan(9, 39, 25)
        ploylon2=fspan(100, 149, 10)
        ploylon3=fspan(280, 315, 25)
        ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
        ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
        ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  

        dum1=new(4,graphic)
        dum2=dum1
        dum3=dum1

        do i=0,2
        ployres@gsLineDashPattern=0
        dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
        ployres@gsLineDashPattern=16
        dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
        ployres@gsLineDashPattern=0
        dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
        end do 

    
        
        pres=True
        pres@gsnFrame=False
        pres@gsnPanelLabelBar=False
        pres@gsnPanelFigureStringsFontHeightF=0.01
        ; pres@pmLabelBarParallelPosF=0.04
        pres@lbLabelStride=1
        pres@gsnPanelFigureStrings="("+(/"a","b","c","d","e","f"/)+")"
        pres@gsnPanelFigureStringsJust="TopLeft"
        pres@gsnPanelXF=(/0.05,0.35,0.65,0.05,0.35,0.65/)
        pres@gsnPanelYF=(/0.8,0.8,0.8,0.55,0.55,0.55/)
        pres@gsnPanelDebug=True

        

        gsn_panel(wks, plot, (/2,3/), pres)
        




  end 

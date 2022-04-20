function decade_trend(data)
    local data_a, data_a_new, rc
    begin
    ;-----------------------------------------------
    rc = regCoef_n(data&time,data,0,0)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data&level
    rc&lon=data&lon
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
    ; print(rc)
return(rc)
end

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
;;****************************************************************
;;main program
;;*************************************************************

begin

    yearstart=197901 
    yearend=200812
    year1=1979
    year2=2008

    f1=addfile("/mnt/f/share/CMIP5/msf_AMIP/msf_time_amip_amip_r1i1p1_197901-200812.nc", "r")
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    msf_time_amip_r=f1->msf_time(:,{90:290},timeind)
    msf_time_amip=tofloat(msf_time_amip_r/10^9)
    copy_VarMeta(msf_time_amip_r, msf_time_amip)
    delete(msf_time_amip_r)
    printVarSummary(msf_time_amip)
    delete(timeind)
    delete(f1)
    delete(date)
    
    f1=addfile("/mnt/f/share/rea_msf/msf_time_ens_global.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    msf_time_rea=f1->msf_time(:,{90:290},timeind)
    msf_time_rea=(/msf_time_rea/10^9/)
    printVarSummary(msf_time_rea)
    delete(timeind)
    delete(f1)
    delete(date)


    models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)


    xi_model=(/"GFDL-CM3","IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/)
    dong_model=(/"MPI-ESM-MR","MIROC5"/)
    dim=dimsizes(models)
    msf_total=new((/dim,12,81,360/), float)
     msf_total_amip=new((/dim,12,81,360/), float)



    ;******************************************************************************************************
    ;循环读入数据
    ;******************************************************************************************************
    do i=0,dimsizes(models)-1
    model=models(i)
    print("**********model*************="+model)
    ; model="ACCESS1-0"
    diri="/mnt/f/share/CMIP5/msf_history/"
    diri1="/mnt/f/share/CMIP5/msf_rcp45/"

    fils1=systemfunc("ls "+diri+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")
    fils2=systemfunc("ls "+diri1+"msf_time_"+model+"_rcp45_r1i1p1_200601-201212.nc")

    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    yearstart1=197901
    yearend1=200512

    yearstart2=200401   ;;;因为存在200512开始的
    yearend2=200812

    year1=1979
    year2=2008

    ;;;********************************************* 
    f=addfiles(fils, "r")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f[:]->time, 1)

    date1=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)


    timeind=ind(date.ge.197901.and.date.le.200812)
    timeind1=ind(date1.le.yearend1.and.date1.ge.yearstart1)
    timeind2=ind(date2.le.yearend2.and.date2.ge.yearstart2)

    delete([/date1 ,date2/])
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind1=get1Dindex(f1->level,z0)
    levind2=get1Dindex(f2->level,z0)
    year=year2-year1+1

    msf_time1_r=f1->msf_time(levind1,{90:290},timeind1)
    msf_time1=tofloat(msf_time1_r)
    copy_VarMeta(msf_time1_r, msf_time1)
    delete(msf_time1_r)

    delete(f1)
    delete(levind1)
    printVarSummary(msf_time1)
    dim_data=dimsizes(msf_time1)
    msf_time=new((/dim_data(0),dim_data(1),year*12/),typeof(msf_time1))

    msf_time2_r=f2->msf_time(levind2,{90:290},timeind2)
    msf_time2=tofloat(msf_time2_r)
    copy_VarMeta(msf_time2_r, msf_time2)
    delete(msf_time2_r)
    delete(f2)
    delete(levind2)
    printVarSummary(msf_time2)



    msf_time(:,:,0:dimsizes(timeind1)-1)=msf_time1
    msf_time(:,:,dimsizes(timeind1):)=msf_time2
    delete(msf_time1)
    delete(msf_time2)
    delete(timeind1)
    delete(timeind2)
    printVarSummary(msf_time)
    delete(f)
    msf_time=(/msf_time/10^9/)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    msf_total(i,:,:,:)=msf_time
    delete(msf_time) 
    delete(date)


    
    diri2="/mnt/f/share/CMIP5/msf_AMIP/"

    fils3=systemfunc("ls "+diri2+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")
    f3=addfile(fils3, "r")
    date=cd_calendar(f3->time, 1)
    timeind1=ind(date.le.200812.and.date.ge.197901)
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind1=get1Dindex(f3->level,z0)
    year=year2-year1+1
    msf_time1_r=f3->msf_time(levind1,{90:290},timeind1)
    msf_time_amip1=tofloat(msf_time1_r)
    copy_VarMeta(msf_time1_r, msf_time_amip1)
    delete(msf_time1_r)
    msf_time_amip1=(/msf_time_amip1/10^9/)
 
    delete(date)
    msf_total_amip(i,:,:,:)=msf_time_amip1
    delete(msf_time_amip1) 
    delete(timeind1)
    delete(f3)
    delete(fils3)
  end do

    indxi =get1Dindex(models,xi_model)
    inddong = get1Dindex(models,dong_model)

    msf_xi=dim_avg_n_Wrap(msf_total(indxi,:,:,:), 0)
    msf_dong=dim_avg_n_Wrap(msf_total(inddong,:,:,:), 0)
    
    msf_amip_xi=dim_avg_n_Wrap(msf_total_amip(indxi,:,:,:), 0)
    msf_amip_dong=dim_avg_n_Wrap(msf_total_amip(inddong,:,:,:), 0)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    msf_clm_rea=dim_avg_n_Wrap(msf_time_rea, 2)
    msf_clm_amip=dim_avg_n_Wrap(msf_time_rea, 2)
   
    
   
    msf_clm_amip_xi=dim_avg_n_Wrap(msf_amip_xi,2)
    msf_clm_amip_dong=dim_avg_n_Wrap(msf_amip_dong,2)

    msf_clm_cmip_xi=dim_avg_n_Wrap(msf_xi,2)
    msf_clm_cmip_dong=dim_avg_n_Wrap(msf_dong,2)

    msf_aunnal_rea=month_to_annual(msf_time_rea(time|:,level|:,lon|:), 1)
    msf_aunnal_rea!0="time"
    msf_aunnal_rea&time=ispan(year1, year2, 1)

    rc_rea=decade_trend(msf_aunnal_rea)
    prob_rea=rc_test(rc_rea)
    rc_rea=(/rc_rea*10.0/)

    
    msf_aunnal_amip=month_to_annual(msf_time_amip(time|:,level|:,lon|:), 1)
    msf_aunnal_amip!0="time"
    msf_aunnal_amip&time=ispan(year1, year2, 1)

    rc_amip=decade_trend(msf_aunnal_amip)
    prob_amip=rc_test(rc_amip)
    rc_amip=(/rc_amip*10.0/)
   
    ;;;;耦合模式分开

    msf_aunnal_xi=month_to_annual(msf_xi(time|:,level|:,lon|:), 1)
    msf_aunnal_xi!0="time"
    msf_aunnal_xi&time=ispan(year1, year2, 1)

    rc_xi=decade_trend(msf_aunnal_xi)
    prob_xi=rc_test(rc_xi)
    rc_xi=(/rc_xi*10.0/)

    
    msf_aunnal_dong=month_to_annual(msf_dong(time|:,level|:,lon|:), 1)
    msf_aunnal_dong!0="time"
    msf_aunnal_dong&time=ispan(year1, year2, 1)

    rc_dong=decade_trend(msf_aunnal_dong)
    prob_dong=rc_test(rc_dong)
    rc_dong=(/rc_dong*10.0/)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;AMIP模式分开
    
    msf_aunnal_amip_xi=month_to_annual(msf_amip_xi(time|:,level|:,lon|:), 1)
    msf_aunnal_amip_xi!0="time"
    msf_aunnal_amip_xi&time=ispan(year1, year2, 1)

    rc_xi_amip=decade_trend(msf_aunnal_amip_xi)
    prob_amip_xi=rc_test(rc_xi_amip)
    rc_xi_amip=(/rc_xi_amip*10.0/)

    
    msf_aunnal_amip_dong=month_to_annual(msf_amip_dong(time|:,level|:,lon|:), 1)
    msf_aunnal_amip_dong!0="time"
    msf_aunnal_amip_dong&time=ispan(year1, year2, 1)

    rc_dong_amip=decade_trend(msf_aunnal_amip_dong)
    prob_amip_dong=rc_test(rc_dong_amip)
    rc_dong_amip=(/rc_dong_amip*10.0/)










;;;画图
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析,amip及cmip加强减弱模式的趋势的MSF图")
;********************************************************
;;draw
;********************************************************

    cmap=read_colormap_file("MPL_PuOr")
    cmap=(/cmap(::-1,0:3)/)

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="Reanalysis"
    res2@gsnLeftStringFontHeightF=0.02
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -0.3
    res2@cnMaxLevelValF = 0.3
    res2@cnLevelSpacingF = 0.05
    res2@cnFillPalette=cmap

    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    res2@lbLabelBarOn=False
    
    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBLabelFontHeightF=0.02
    res2@tmYLLabelFontHeightF=0.02
    res2@tmYLMajorThicknessF=2.0
    res2@tmXBMajorThicknessF=2.0
    res2@tmBorderThicknessF=4.0



    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
    res2@tmXBValues=ispan(90, 270, 30)
    res2@tmXBLabels=(/"90E","120E","150E","180","150W","120W","90W"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(90,270,10)
    res2@tmYRMode="Automatic"
    res2@vpHeightF=0.4
    res2@vpWidthF=0.6
   
    plot=new(6,graphic)

    plot(0)=gsn_csm_pres_hgt(wks,rc_rea,res2)
    res2@gsnLeftString="All AMIP"
    plot(1)=gsn_csm_pres_hgt(wks,rc_amip,res2)
    res2@gsnLeftString="CMIP5(+PWCI)"
    plot(2)=gsn_csm_pres_hgt(wks,rc_xi,res2)
    res2@gsnLeftString="CMIP5(-PWCI)"
    plot(3)=gsn_csm_pres_hgt(wks,rc_dong,res2)
    res2@gsnLeftString="AMIP(+CMIP5)"
    plot(4)=gsn_csm_pres_hgt(wks,rc_xi_amip,res2)
    res2@gsnLeftString="AMIP(-CMIP5)"
    plot(5)=gsn_csm_pres_hgt(wks,rc_dong_amip,res2)

  
;;;gsn_contour_shade
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@gsnRightString=""
        res_pro@gsnLeftString=""
        res_pro@cnLevelSelectionMode = "ExplicitLevels"
        res_pro@cnLinesOn=False
        res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing
        res_pro@cnLineLabelsOn=False
        res_pro@tmYRMode = "Automatic"
        res_pro@cnInfoLabelOn=False
        

        iplot=new(6,graphic)

        iplot(0)=gsn_csm_pres_hgt(wks,prob_rea,res_pro)
        iplot(1)=gsn_csm_pres_hgt(wks,prob_amip,res_pro)
        
        iplot(2)=gsn_csm_pres_hgt(wks,prob_xi,res_pro)
        iplot(3)=gsn_csm_pres_hgt(wks,prob_dong,res_pro)
        iplot(4)=gsn_csm_pres_hgt(wks,prob_amip_xi,res_pro)
        iplot(5)=gsn_csm_pres_hgt(wks,prob_amip_dong,res_pro)

        

        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        opt@gsnShadeFillScaleF   = 1.0         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0025     ; make dots larger
        
        do i=0,5
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.1, opt)
       end do

 ;;;    
        res2@gsnRightString=""
        res2@gsnLeftString=""
        res2@cnLevelSelectionMode = "ExplicitLevels"
        res2@cnLinesOn=True
        res2@cnLevels = ispan(-6,6,1)
        res2@cnLineLabelsOn=True
        res2@cnLineLabelBackgroundColor="Transparent"
        res2@cnLineThicknessF = 2.
        res2@cnFillOn = False
        res2@gsnContourNegLineDashPattern=1
        res2@gsnContourZeroLineThicknessF=4.0
        
        plot1=new(6,graphic)

        plot1(0) = gsn_csm_pres_hgt(wks, msf_clm_rea , res2)
        plot1(1) = gsn_csm_pres_hgt(wks, msf_clm_amip , res2)

        plot1(2) = gsn_csm_pres_hgt(wks, msf_clm_cmip_xi , res2)
        plot1(3) = gsn_csm_pres_hgt(wks, msf_clm_cmip_dong , res2)

        plot1(4) = gsn_csm_pres_hgt(wks, msf_clm_amip_xi , res2)
        plot1(5) = gsn_csm_pres_hgt(wks, msf_clm_amip_dong , res2)
        
       
        ;;;添加地形
        ployres                  = True          ; polyline mods desired
        ployres@gsLineThicknessF = 15.0            ; line thickness   
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

        dum1=new(6,graphic)
        dum2=dum1
        dum3=dum1
    
        do i=0,5
        ployres@gsLineDashPattern=0
        dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
        ployres@gsLineDashPattern=16
        dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
        ployres@gsLineDashPattern=0
        dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
        end do

        do i=0,5
        overlay(plot(i),plot1(i))
        overlay(plot(i),iplot(i))
        end do
        pres=True
        pres@gsnPanelLabelBar=True
        pres@pmLabelBarWidthF=0.7
        pres@gsnPanelFigureStrings="("+(/"a","b","c","d","e","f"/)+")"
        pres@gsnPanelFigureStringsJust="TopLeft"
             gsn_panel(wks, plot,(/3,2/), pres)
        end 
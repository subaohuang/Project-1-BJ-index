
function one_test(pra_rc)
begin
    
    n=dimsizes(pra_rc)
    var_pra_rc    =  dim_variance_n(pra_rc(:,:,:), 0)    ; 13.477
    var_pra_rc  =  sqrt(var_pra_rc/n(0)) 
    var_pra_rc  =where(var_pra_rc.eq.0, var_pra_rc@_FillValue, var_pra_rc)
    
    pra_rc_mean=dim_avg_n(pra_rc, 0)
    t_pra_rc = pra_rc_mean/var_pra_rc
    confi_pr_rc=student_t(t_pra_rc , n(0)-1)
 return confi_pr_rc
 end

;;;********************************************************************
;;;主程序
;;;******************************************************************
begin



; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FGOALS-g2","FGOALS-s2","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
 ; dim=dimsizes(models)
;;;强度
 ;;str model all:

 diri="/mnt/f/share/CMIP5/msf_history/"


xi_model=(/"CSIRO-Mk3-6-0","INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
 ; xi_model=(/"CSIRO-Mk3-6-0"/)
; xi_model=(/"INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
dong_model=(/"MPI-ESM-MR","MIROC5"/)
;;;全球
latmin=-30
latmax=30
lonmin=0
lonmax=360
lat=fspan(30, -30, 25)
lon=fspan(0, 359.5, 144)
lat_pr=fspan(30, -30, 25)
lon_pr=fspan(0, 359.5, 144)
;;;



psl_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/),float)
psl_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/),float)


pr_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat_pr),dimsizes(lon_pr)/) ,float)
pr_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat_pr),dimsizes(lon_pr)/),float)

uwnd_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/) ,float)
uwnd_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/), float)

vwnd_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/), float)
vwnd_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/),float)


do i=0,1
    if i.eq.0 then 
     dim=dimsizes(xi_model)   
     models=xi_model
 else
    dim=dimsizes(dong_model)
    models=dong_model
   end if

 do j=0,dim-1
    
    model=models(j)
    yearstart=197901
    yearend=200812
    
    year1=1979
    year2=2008
    
   ;;;读入降水数据*************************************************************************************
    print("******************"+model+"*********************")
    fils1="/mnt/i/CMIP5/history/pr_chazhi/pr_"+model+"_historical_r1i1p1_197901-200512.nc"
    fils2= "/mnt/i/CMIP5/rcp45/pr_chazhi/pr_"+model+"_rcp45_r1i1p1_200601-201212.nc"

    printVarSummary(fils1)
    printVarSummary(fils2)
    
    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils1)
    delete(fils2)
     
    time=f1[:]->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    pr=f1[:]->pr(timeind,{30:-30},{0:360})
    printVarSummary(pr)
    ;;;换算成mm/day
    pr=(/pr*3600*24.0/)

    delete(f1)
    delete(pr&time)
    pr&time=date(timeind)

    delete(time)
    delete(date)
    printVarSummary(pr)
    delete(timeind)
    
    pr_annual=month_to_annual(pr, 1)
    pr_annual!0="time"
    pr_annual&time=ispan(year1,year2,1)
    delete(pr)


    ;***********************************************************************************************
    ;;
    ;;;**********************************************************************************************
    print("******************读入风场*********************")
    fils1=systemfunc("ls /mnt/i/CMIP5/history/uv/uv_"+model+"_history_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/rcp45/uv/uv_"+model+"_rcp45_r1i1p1_"+"*"+".nc")
    print(fils1)
    print(fils1)

    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils1)
    delete(fils2)
     
    time=f1[:]->time    
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd=f1[:]->uwnd(timeind,{30:-30},{0:360})  ;;;
    vwnd=f1[:]->vwnd(timeind,{30:-30},{0:360})  ;;;
    delete(f1)
   
;;;;
    delete(uwnd&time)
    delete(vwnd&time)
    uwnd&time=date(timeind)
    vwnd&time=date(timeind)
     
    uwnd_annual=month_to_annual(uwnd, 1)
    uwnd_annual!0="time"
    uwnd_annual&time=ispan(year1,year2,1)
    vwnd_annual=month_to_annual(vwnd, 1)
    vwnd_annual!0="time"
    vwnd_annual&time=ispan(year1,year2,1)
    delete(uwnd)
    delete(vwnd)
    delete(timeind)
    delete(time)
    delete(date)
;;;**************************************************************************************
;;;;**********************************************************************************
;;;;***************************************************************************


   print("******************读入海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/i/CMIP5/history/slp_chazhi/psl_"+model+"_historical_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/RCP45/slp_chazhi/psl_"+model+"_rcp45_r1i1p1_"+"*"+".nc")
    print(fils1)
    print(fils2)


    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils1)
    delete(fils2)
   
    time=f1[:]->time
    date=cd_calendar(time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    printVarSummary(timeind)
    psl=lonPivot(f1[:]->psl(0:359,:,:),0)   
    delete(f1)
    psl=(/psl/)/100.0
   
     
    psl_annual=month_to_annual(psl(:,{30:-30},{0:360}), 1)
    psl_annual!0="time"
    psl_annual&time=ispan(year1,year2,1)
    delete(timeind)
    delete(psl)
    delete(time)
    delete(date)

;;;***********************************************************************************************
;;;
;;;;'**********************************************************************************************'
    print("*****************计算趋势************************")  

    pr_rc1=regCoef_n(ispan(1,dimsizes(pr_annual&time),1),pr_annual(:,:,:),0,0)  
    copy_VarCoords(pr_annual(0,:,:), pr_rc1)
    delete(pr_annual)
    printVarSummary(pr_rc1)
 
    

    uwnd_rc1=regCoef_n(ispan(1,dimsizes(uwnd_annual&time),1),uwnd_annual(:,:,:),0,0)  
    copy_VarCoords(uwnd_annual(0,:,:), uwnd_rc1)
    delete(uwnd_annual)
    printVarSummary(uwnd_rc1)


    vwnd_rc1=regCoef_n(ispan(1,dimsizes(vwnd_annual&time),1),vwnd_annual(:,:,:),0,0)  
    copy_VarCoords(vwnd_annual(0,:,:), vwnd_rc1)
    delete(vwnd_annual)
    printVarSummary(vwnd_rc1)

    
    psl_rc1=regCoef_n(ispan(1,dimsizes(psl_annual&time),1),psl_annual(:,:,:),0,0)  
    copy_VarCoords(psl_annual(0,:,:), psl_rc1)
    delete(psl_annual)
    printVarSummary(psl_rc1)   

    
    if i.eq.0 then 
      psl_rc_ximodel(j,:,:)=psl_rc1
      pr_rc_ximodel(j,:,:)=pr_rc1
      uwnd_rc_ximodel(j,:,:)=uwnd_rc1
      vwnd_rc_ximodel(j,:,:)=vwnd_rc1
   else
      psl_rc_dongmodel(j,:,:)=psl_rc1
      pr_rc_dongmodel(j,:,:)=pr_rc1
      uwnd_rc_dongmodel(j,:,:)=uwnd_rc1
      vwnd_rc_dongmodel(j,:,:)=vwnd_rc1
    end if

    delete([/psl_rc1,pr_rc1,uwnd_rc1,vwnd_rc1/])
end do
delete(dim)
delete(models)
end do

psl_rc1_xi_mean=dim_avg_n(psl_rc_ximodel,0)
pr_rc1_xi_mean=dim_avg_n(pr_rc_ximodel,0)
uwnd_rc1_xi_mean=dim_avg_n(uwnd_rc_ximodel,0)
vwnd_rc1_xi_mean=dim_avg_n(vwnd_rc_ximodel,0)

psl_rc1_dong_mean=dim_avg_n(psl_rc_dongmodel,0)
pr_rc1_dong_mean=dim_avg_n(pr_rc_dongmodel,0)
uwnd_rc1_dong_mean=dim_avg_n(uwnd_rc_dongmodel,0)
vwnd_rc1_dong_mean=dim_avg_n(vwnd_rc_dongmodel,0)


psl_rc1_xi_mean!0="lat"
psl_rc1_xi_mean!1="lon"

psl_rc1_xi_mean&lat=lat
psl_rc1_xi_mean&lon=lon
psl_rc1_xi_mean&lat@units="degrees_north"
psl_rc1_xi_mean&lon@units="degrees_east"
printVarSummary(psl_rc1_xi_mean)

pr_rc1_xi_mean!0="lat"
pr_rc1_xi_mean!1="lon"

pr_rc1_xi_mean&lat=lat_pr
pr_rc1_xi_mean&lon=lon_pr
pr_rc1_xi_mean&lat@units="degrees_north"
pr_rc1_xi_mean&lon@units="degrees_east"

printVarSummary(psl_rc1_xi_mean)


copy_VarCoords(psl_rc1_xi_mean, psl_rc1_dong_mean)
copy_VarCoords(pr_rc1_xi_mean, pr_rc1_dong_mean)
copy_VarCoords(psl_rc1_xi_mean, uwnd_rc1_xi_mean)
copy_VarCoords(psl_rc1_xi_mean, uwnd_rc1_dong_mean)
copy_VarCoords(psl_rc1_xi_mean, vwnd_rc1_xi_mean)
copy_VarCoords(psl_rc1_xi_mean, vwnd_rc1_dong_mean)


prob1=one_test(psl_rc_ximodel)
copy_VarCoords(psl_rc1_xi_mean, prob1)


prob2=one_test(psl_rc_dongmodel)
copy_VarCoords(psl_rc1_dong_mean, prob2)



uwnd_prob1=one_test(uwnd_rc_ximodel)
vwnd_prob1=one_test(vwnd_rc_ximodel)
uwnd_rc1_xi_mean=where(uwnd_prob1.le.0.1.or.vwnd_prob1.le.0.1, uwnd_rc1_xi_mean, uwnd_rc1_xi_mean@_FillValue)

uwnd_prob2=one_test(uwnd_rc_ximodel)
vwnd_prob2=one_test(vwnd_rc_dongmodel)
uwnd_rc1_dong_mean=where(uwnd_prob2.le.0.1.or.vwnd_prob2.le.0.1, uwnd_rc1_dong_mean, uwnd_rc1_dong_mean@_FillValue)




psl_rc1_xi_mean=(/psl_rc1_xi_mean/)*10.0
psl_rc1_dong_mean=(/psl_rc1_dong_mean/)*10.0

pr_rc1_xi_mean=(/pr_rc1_xi_mean/)*10.0
pr_rc1_dong_mean=(/pr_rc1_dong_mean/)*10.0

uwnd_rc1_xi_mean=(/uwnd_rc1_xi_mean/)*10.0
uwnd_rc1_dong_mean=(/uwnd_rc1_dong_mean/)*10.0

vwnd_rc1_xi_mean=(/vwnd_rc1_xi_mean/)*10.0
vwnd_rc1_dong_mean=(/vwnd_rc1_dong_mean/)*10.0
;;*********************************************************************************************
;;读入观测的海温和风场的情况
;;;*********************************************************************************************
    yearstart=197901
    yearend=200812

;;;
    print("*****************读入MRE的风场和海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_u_v_chazhi_1000.nc")

    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd_mre_r=f1->uwnd(timeind,{30:-30},{0:360})
    vwnd_mre_r=f1->vwnd(timeind,{30:-30},{0:360})

    delete(uwnd_mre_r&time)
    delete(vwnd_mre_r&time)
    uwnd_mre_r&time=date(timeind)
    vwnd_mre_r&time=date(timeind)
    
    uwnd_mre=month_to_annual(uwnd_mre_r, 1)
    vwnd_mre=month_to_annual(vwnd_mre_r, 1)

    uwnd_mre!0="time"
    vwnd_mre!0="time"
    uwnd_mre&time=ispan(year1, year2,1)
    vwnd_mre&time=ispan(year1, year2,1)


    uwnd_mre&lat@units="degree_north"
    uwnd_mre&lon@units="degree_east"
    vwnd_mre&lat@units="degree_north"
    vwnd_mre&lon@units="degree_east"

    delete(uwnd_mre_r)
    delete(vwnd_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)

  
    ;;;气压场
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_slp_chazhi.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    slp_r=f1->slp(timeind,{30:-30},{0:360})
    delete(slp_r&time)
    slp_r&time=date(timeind)

    slp_mre=month_to_annual(slp_r, 1)

    slp_mre!0="time"  
    slp_mre&time=ispan(year1, year2,1)

    slp_mre&lat@units="degree_north"
    slp_mre&lon@units="degree_east"
    delete(slp_r)
    delete(date)
    delete(f1)
    delete(fils1)
;;;
  ;;;降水场
    fils1=systemfunc("ls /mnt/f/share/GPCP/precip.mon.mean.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    pr_mre_r=f1->precip(timeind,{30:-30},{0:360})

    delete(pr_mre_r&time)
    printVarSummary(pr_mre_r)
    pr_mre_r&time =  date(timeind)
    pr_mre=month_to_annual(pr_mre_r, 1)
    
    pr_mre!0="time"
    pr_mre&time=ispan(year1, year2,1)
    pr_mre&lat@units="degree_north"
    pr_mre&lon@units="degree_east"

    delete(pr_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)
  
    ;;;
    print("*****************计算趋势************************")  

    pr_rc1_mre=regCoef_n(ispan(1,dimsizes(pr_mre&time),1),pr_mre(:,:,:),0,0)  
    copy_VarCoords(pr_mre(0,:,:), pr_rc1_mre)
    printVarSummary(pr_rc1_mre)
    
    pr_rc1_mre=(/pr_rc1_mre/)*10

    uwnd_rc1_mre=regCoef_n(ispan(1,dimsizes(uwnd_mre&time),1),uwnd_mre(:,:,:),0,0)  
    copy_VarCoords(uwnd_mre(0,:,:), uwnd_rc1_mre)
    printVarSummary(uwnd_rc1_mre)
  
    uwnd_rc1_mre=(/uwnd_rc1_mre/)*10
      
    vwnd_rc1_mre=regCoef_n(ispan(1,dimsizes(vwnd_mre&time),1),vwnd_mre(:,:,:),0,0)  
    copy_VarCoords(vwnd_mre(0,:,:), vwnd_rc1_mre)
    printVarSummary(vwnd_rc1_mre)
    vwnd_rc1_mre=(/vwnd_rc1_mre/)*10
    ;;;***********************************
    slp_rc1_mre=regCoef_n(ispan(1,dimsizes(slp_mre&time),1),slp_mre(:,:,:),0,0)  
    copy_VarCoords(slp_mre(0,:,:), slp_rc1_mre)
    printVarSummary(slp_rc1_mre)
    slp_rc1_mre=(/slp_rc1_mre/)*10


 print("**********************画图*****************************")

wks1=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/tos_group/变强变弱模式的集合以及再分析平均结果全球")
 

    ; colormap=read_colormap_file("GMT_jet")
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)

    ; cmap(0:6,:)=colormap(2:92:15,:)
    ; ; cmap(0:6,:)=colormap(2:50:8,:)
    ; cmap(7:17,:)=colormap(105:255:15,:)

    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterString ="CMIP5(+PWC trend)"
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
  
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnFillDrawOrder="PreDraw"
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False
    res@mpFillDrawOrder="PostDraw"
    ; res@mpLandFillColor="Transparent"

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@vpWidthF=0.7
    res@vpHeightF=0.35
  

    ;;;LineArrow
    res@gsnScalarContour        = True
    res@vcRefMagnitudeF          = 0.4 
    res@vcRefLengthF             = 0.03     ; define length of vec ref
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.03
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01


    ;;;CurlyVector

    res@vcLineArrowThicknessF=1.5
    res@vcMinFracLengthF=0.3
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
        


    plot=new(3,graphic)
    plot(0) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_xi_mean,vwnd_rc1_xi_mean,psl_rc1_xi_mean,res)
     res@gsnCenterString ="CMIP5(-PWC trend)"
    plot(1) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_dong_mean,vwnd_rc1_dong_mean,psl_rc1_dong_mean,res)
     res@gsnCenterString ="Observation "
     res@vcRefMagnitudeF          = 0.6
    plot(2) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_mre,vwnd_rc1_mre,slp_rc1_mre,res)

         ;;;打点
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@cnLinesOn=False
        res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

        res_pro@cnLineLabelsOn=False
        res_pro@tmYRMode = "Automatic"
        res_pro@cnInfoLabelOn=False

        iplot=new(2,graphic)

        iplot(0)=gsn_csm_contour(wks1,prob1,res_pro)
        iplot(1)=gsn_csm_contour(wks1,prob2,res_pro)
     
        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 0.2         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
        do i=0,1
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
        ;;;  
        overlay(plot(i), iplot(i))
        end do

        ;;;画slp的等值线
        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -1.2        ; set min contour level
        res_slp@cnMaxLevelValF      = 1.2       ; set max contour level
        res_slp@cnLevelSpacingF     = 0.6       ; set contour spacing
        res_slp@cnLineLabelsOn     =True
        res_slp@cnLineLabelBackgroundColor="Transparent"
        res_slp@gsnContourNegLineDashPattern=2
        res_slp@cnLineLabelInterval        = 2             ; default = 2
        res_slp@cnLineLabelPlacementMode   = "constant"    ; choose constant label method
        res_slp@cnLineLabelFormat= "0@;*.1f" 
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineColor="black"
        res_slp@cnLineThicknessF=1.5

        plot_slp=new(3,graphic)

        
        plot_slp(0)=gsn_csm_contour(wks1,pr_rc1_xi_mean,res_slp)
        plot_slp(1)=gsn_csm_contour(wks1,pr_rc1_dong_mean,res_slp)
        res_slp@cnLevelSpacingF     = 0.1      ; set contour spacing
        plot_slp(2)=gsn_csm_contour(wks1,pr_rc1_mre ,res_slp)

        do i=0,2
        overlay(plot(i), plot_slp(i))
    end do




pres=True
pres@gsnPanelLabelBar=True
pres@gsnPanelFigureStrings="("+(/"a","b","c"/)+")"
pres@gsnPanelFigureStringsJust="TopLeft"
gsn_panel(wks1, plot, (/3,1/), pres)


end
    
  

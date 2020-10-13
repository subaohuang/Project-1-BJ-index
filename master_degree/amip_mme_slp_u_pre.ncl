function  one_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end


;;;********************************************************************
;;;主程序
;;;******************************************************************
begin



diri="/mnt/f/share/CMIP5/msf_history/"
 ; xi_model = (/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
 ;         "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
 ;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
 xi_model = (/"ACCESS1-3"/)
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
;;;;AMIP
psl_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/),float)
pr_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat_pr),dimsizes(lon_pr)/) ,float)
uwnd_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/) ,float)
vwnd_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/), float)

;;;
   dim=dimsizes(xi_model)   
     models=xi_model


 do j=0,dim-1
    
    model=models(j)
    yearstart=197901
    yearend=200812
    
    year1=1979
    year2=2008
    
   ;;;读入降水数据************************************************************************************
    fils3= "/mnt/f/share/CMIP5/AMIP/pr_chazhi/pr_"+model+"_AMIP_r1i1p1_197901-200812.nc"
    f3=addfile(fils3, "r")
    time=f3->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    pr_amip=f3->pr(timeind,{30:-30},{0:360})
    delete(f3)
    delete(fils3)
    ;;;换算成mm/day
    pr_amip=(/pr_amip*3600*24.0/)

    pr_amip_djf  = month_to_season(pr_amip, "DJF")
    pr_amip_mam  = month_to_season(pr_amip, "MAM")
    
    pr_amip_annual=month_to_annual(pr_amip, 1)
    pr_amip_annual = (/(pr_amip_djf + pr_amip_mam)/2.0/)

    pr_amip_annual!0="time"
    pr_amip_annual&time=ispan(year1,year2,1)
    delete(pr_amip)
    delete(time)
    delete(date)

    ;;;AMIP的风场数据
    fils3=systemfunc("ls /mnt/f/share/CMIP5/AMIP/uv/uv_"+model+"_amip_r1i1p1_"+"*"+".nc") 
    f3=addfile(fils3, "r")
    time=f3->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    uwnd_amip=f3->uwnd(timeind,{30:-30},{0:360})
    vwnd_amip=f3->vwnd(timeind,{30:-30},{0:360})
    delete(f3)
    delete(fils3)
    ;;;换算成mm/day
    uwnd_amip_djf  = month_to_season(uwnd_amip, "DJF")
    uwnd_amip_mam  = month_to_season(uwnd_amip, "MAM")
    
    uwnd_amip_annual = month_to_annual(uwnd_amip, 1)
    uwnd_amip_annual = (/(uwnd_amip_djf + uwnd_amip_mam)/2.0/)

    uwnd_amip_annual!0="time"
    uwnd_amip_annual&time=ispan(year1,year2,1)
    delete(uwnd_amip)

    vwnd_amip_djf  = month_to_season(vwnd_amip, "DJF")
    vwnd_amip_mam  = month_to_season(vwnd_amip, "MAM")

    vwnd_amip_annual = month_to_annual(vwnd_amip, 1)
    vwnd_amip_annual = (/(vwnd_amip_djf + vwnd_amip_mam)/2.0/)

    vwnd_amip_annual!0="time"
    vwnd_amip_annual&time=ispan(year1,year2,1)
    delete(vwnd_amip)
    delete(timeind)
    delete(time)
    delete(date)
    ;;;;

    ;;;气压的数据
    fils3= "/mnt/f/share/CMIP5/AMIP/slp_chazhi/psl_"+model+"_amip_r1i1p1_197901-200812.nc"
    f3=addfile(fils3, "r")
    time=f3->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    psl_amip=lonPivot(f3->psl(0:359,:,:),0) 
    psl_amip=(/psl_amip/100.0/)
    delete(f3)
    delete(fils3)

    psl_amip_djf  = month_to_season(psl_amip(:,{30:-30},:), "DJF")
    psl_amip_mam  = month_to_season(psl_amip(:,{30:-30},:), "MAM")
    
    psl_amip_annual=month_to_annual(psl_amip(:,{30:-30},:), 1)
    psl_amip_annual = (/(psl_amip_djf + psl_amip_mam)/2.0/)

    psl_amip_annual!0="time"
    psl_amip_annual&time=ispan(year1,year2,1)
    delete(psl_amip)
    delete(time)
    delete(timeind)
    delete(date)
    ;;;
    ;;;;

;;;***********************************************************************************************
;;;
;;;;'**********************************************************************************************'
 printVarSummary(psl_amip_annual)


      psl_amip_ximodel(j,:,:,:)=psl_amip_annual
      pr_amip_ximodel(j,:,:,:)=pr_amip_annual
      uwnd_amip_ximodel(j,:,:,:)=uwnd_amip_annual
      vwnd_amip_ximodel(j,:,:,:)=vwnd_amip_annual

    delete([/psl_amip_annual,pr_amip_annual,uwnd_amip_annual,vwnd_amip_annual/])
end do
delete(dim)
delete(models)



psl_xi_mean_amip=dim_avg_n(psl_amip_ximodel,0)
pr_xi_mean_amip=dim_avg_n(pr_amip_ximodel,0)
uwnd_xi_mean_amip=dim_avg_n(uwnd_amip_ximodel,0)
vwnd_xi_mean_amip=dim_avg_n(vwnd_amip_ximodel,0)


psl_xi_mean_amip!0="time"
psl_xi_mean_amip!1="lat"
psl_xi_mean_amip!2="lon"

psl_xi_mean_amip&time=ispan(year1, year2, 1)
psl_xi_mean_amip&lat=lat
psl_xi_mean_amip&lon=lon
psl_xi_mean_amip&lat@units="degrees_north"
psl_xi_mean_amip&lon@units="degrees_east"
printVarSummary(psl_xi_mean_amip)

copy_VarCoords(psl_xi_mean_amip, pr_xi_mean_amip)
copy_VarCoords(psl_xi_mean_amip, uwnd_xi_mean_amip)
copy_VarCoords(psl_xi_mean_amip, vwnd_xi_mean_amip)




    print("*****************计算趋势************************")  
  ;;;AMIP实验
    pr_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(pr_xi_mean_amip&time),1),pr_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(pr_xi_mean_amip(0,:,:), pr_rc1_xi_mean_amip)
    delete(pr_xi_mean_amip)
    printVarSummary(pr_rc1_xi_mean_amip)


    uwnd_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(uwnd_xi_mean_amip&time),1),uwnd_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(uwnd_xi_mean_amip(0,:,:), uwnd_rc1_xi_mean_amip)
    delete(uwnd_xi_mean_amip)
    printVarSummary(uwnd_rc1_xi_mean_amip)

  
    vwnd_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(vwnd_xi_mean_amip&time),1),vwnd_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(vwnd_xi_mean_amip(0,:,:), vwnd_rc1_xi_mean_amip)
    delete(vwnd_xi_mean_amip)
    printVarSummary(vwnd_rc1_xi_mean_amip)

    
    psl_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(psl_xi_mean_amip&time),1),psl_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(psl_xi_mean_amip(0,:,:), psl_rc1_xi_mean_amip)
    printVarSummary(pr_rc1_xi_mean_amip)
    delete(psl_xi_mean_amip)

    prob1_amip=one_test(psl_rc1_xi_mean_amip)
    
    uwnd_xi_prob1=one_test(uwnd_rc1_xi_mean_amip)
    vwnd_xi_prob1=one_test(vwnd_rc1_xi_mean_amip)
    

    uwnd_rc1_xi_mean_amip=where(uwnd_xi_prob1.le.0.1, uwnd_rc1_xi_mean_amip, uwnd_rc1_xi_mean_amip@_FillValue)
    psl_rc1_xi_mean_amip=(/psl_rc1_xi_mean_amip/)*10.0
    pr_rc1_xi_mean_amip=(/pr_rc1_xi_mean_amip/)*10.0
    uwnd_rc1_xi_mean_amip=(/uwnd_rc1_xi_mean_amip/)*10.0
    vwnd_rc1_xi_mean_amip=(/vwnd_rc1_xi_mean_amip/)*10.0




 print("**********************画图*****************************")

     wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/amip_mme_slp_u_pre")
 

 
      cmap=read_colormap_file("MPL_PuOr")
      cmap=(/cmap(::-1,0:3)/)

    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
  
    res@cnFillPalette="GMT_polar"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnFillDrawOrder="PreDraw"
    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4
    res@cnLevelSpacingF      =  0.05
    
    res@tmXTOn=True
    res@tmYROn=True
    res@tmXUseBottom=False
    res@tmXBLabelsOn=False
    res@tmXTMinorOn=False
    res@tmXBMinorOn=False
    res@tmYLMinorOn=False

    res@tmXBLabelFontHeightF=0.012
    res@tmYLLabelFontHeightF=0.012

    res@tmXBMajorOutwardLengthF=0
    res@tmXBMajorLengthF=0.007

    res@tmXTMajorOutwardLengthF=0
    res@tmXTMajorLengthF=0.007

    res@tmYLMajorOutwardLengthF=0
    res@tmYLMajorLengthF=0.01

    res@tmYLMajorThicknessF=2.0
    res@tmXTMajorThicknessF=2.0
    res@tmBorderThicknessF=2.5

    
    res@tmXTMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXTValues=ispan(0,360,60)
    res@tmXTLabelsOn=True
    res@tmXTLabels=(/"0E","60~S~o","120~S~o","180~S~o","120~S~o","60~S~o","0E"/)
    res@tmXTLabelDeltaF=-0.6
    
    
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLValues=ispan(30,-30,10)
    res@tmYLLabels=(/"30~S~o","20~S~o","10~S~o","0","10~S~o","20~S~o","30~S~o"/)
    res@tmYLLabelDeltaF=-0.6


    res@mpFillDrawOrder="PostDraw"
    res@mpFillOn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@vpWidthF=0.65
    res@vpHeightF=0.35
  

    ;;;LineArrow
    res@gsnScalarContour        = True
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.03
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"
    ; res@vcRefAnnoBackgroundColor="Transparent"
    ; res@vcRefAnnoPerimOn=False
    
    res@vcRefAnnoString1="0.6m/s"
    res@vcRefAnnoFontHeightF=0.01

    ;;;CurlyVector
    res@vcLineArrowThicknessF=1.5
    res@vcMinFracLengthF=0.3
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    

  
    res@vcRefMagnitudeF          = 0.4 
    res@vcRefLengthF             = 0.04     ; define length of vec ref
    res@vcRefAnnoString1="0.4m/s"
    res@tmYLLabelsOn=False
    plot = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_xi_mean_amip,vwnd_rc1_xi_mean_amip,psl_rc1_xi_mean_amip ,res)
  

     

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


        iplot=gsn_csm_contour(wks1,prob1_amip,res_pro)
      
        
     
        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 1.0         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0025      ; make dots larger
        i=0
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
        ;;;  
        overlay(plot(i), iplot(i))
        

        ;;;画slp的等值线
        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -1.0        ; set min contour level
        res_slp@cnMaxLevelValF      = 1.0      ; set max contour level
        res_slp@cnLevelSpacingF     = 0.2      ; set contour spacin

        res_slp@cnLineLabelsOn     =False
        res_slp@cnLineLabelBackgroundColor="Transparent"

        
        colors=new(11,string)
    
    
        colors(0:4)="blueviolet"
 
        colors(5)="Transparent"
   
        colors(6:)="green"
     

        res_slp@cnLineColors=colors
        res_slp@cnMonoLineColor=False
        res_slp@gsnContourNegLineDashPattern=12
        res_slp@cnLineLabelInterval        = 2             ; default = 2
        res_slp@cnLineLabelPlacementMode   = "constant"    ; choose constant label method
        res_slp@cnLineLabelFormat= "0@;*.1f" 
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineThicknessF=2.5
        res_slp@gsnContourZeroLineThicknessF=2.0


        plot_slp=gsn_csm_contour(wks1,pr_rc1_xi_mean_amip,res_slp)

        i = 0
        overlay(plot(i), plot_slp(i))



 ; drawNDCGrid(wks1)
      draw(plot)
      frame(wks1)
end
    
  

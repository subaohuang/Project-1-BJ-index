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

function  sst_cold_annual(sst,timeind,year1,year2)
begin 

    ;;;先得到秋季和冬季的SONDJF海温
    sst_son=month_to_season(sst(time|:,lat|:,lon|:), "SON")
    sst_djf=month_to_season(sst(time|:,lat|:,lon|:), "DJF")
    printVarSummary(sst_son)

    sst_son&time=ispan(year1,year2,1)
    sst_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195lonmin1,195lonmin2的平均
    ;;注意点
    sst_djf(0:dimsizes(sst_djf&time)-2,:,:)=(/sst_djf(1:,:,:)/)
    sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/sst(dimsizes(timeind)-1,:,:)/)
    
    ;;;
    sst_mean_annual=new((/dimsizes(timeind)/12,dimsizes(sst&lat),dimsizes(sst&lon)/),typeof(sst))
    printVarSummary(sst_mean_annual)
    ; sst_mean_annual=month_to_annual(sst_mean, 1)
     ; sst_mean_annual=(/sst_son(:,:,:)+sst_djf(:,:,:)/)*0.5
     sst_mean_annual=(/sst_djf(:,:,:)/)*1.0   ;取冬季

    delete(sst_son)
    delete(sst_djf)

    sst_mean_annual!0="time"
    sst_mean_annual!1="lat"
    sst_mean_annual!2="lon"
    sst_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    sst_mean_annual&lon=sst&lon
    sst_mean_annual&lat=sst&lat
    printVarSummary(sst_mean_annual)
return sst_mean_annual
end
   
function draw_cmip_sst(xi_model,wks)
 begin
;;读入数据
  
  
    yearstart=197901
    yearend=200512
    
    year1=1979
    year2=2005
    
    latmin=-50
    latmax=50
    lonmin=90
    lonmax=300

     if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC".or.xi_model.eq."HadGEM2-ES" then

    yearstart=197901
    yearend=200412

    year1=1979
    year2=2004
    end if

    fils1=systemfunc("ls /mnt/j/CMIP5/history/tos_new/"+xi_model+"/r1i1p1/tos_Omon_"+xi_model+"_historical_r1i1p1_"+"*"+".nc")
    printVarSummary(fils1)
    print("******************"+xi_model+"*********************")
    f1=addfiles(fils1, "join")
    delete(fils1)
    vNames= getfilevarnames(f1[:]) ; get variable names of all groups on file
    print(vNames)
    
    if (any(vNames.eq."LONN179_180")) then 


    time=f1[:]->TIME
 
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    if xi_model.eq."FGOALS-s2" then
    sst_r=f1[:]->TOS(timeind,:,:)
    else
    sst_r=lonFlip(f1[:]->TOS(timeind,:,:))
    end if    
    sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r)
    dNames = getvardims(f1[:]->TOS)
      delete(f1)
    print(dNames)
    sst_r!0="time"
    sst_r!1="lat"
    sst_r!2="lon"
    end if

  
    delete(vNames)
    delete(dNames)

    delete(sst_r&time)
    sst_r&time=date(timeind)
      if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC".or.xi_model.eq."HadGEM2-ES" then
    sst=sst_cold_annual(sst_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2004)
    else 
    sst=sst_cold_annual(sst_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    end if
    delete(sst_r)
    delete(time)
    delete(date)
    printVarSummary(sst)
    delete(timeind)
 
    print("******************读入风场*********************")
    fils1=systemfunc("ls /mnt/j/CMIP5/history/uv/uv_"+xi_model+"_history_r1i1p1_"+"*"+".nc")

    print(fils1)
    f1=addfile(fils1, "r")

    time=f1->time    
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd_r=f1->uwnd(timeind,:,:)  ;;;
    vwnd_r=f1->vwnd(timeind,:,:)  ;;;
    delete(f1)

    
;;;;
    delete(uwnd_r&time)
    delete(vwnd_r&time)
    uwnd_r&time=date(timeind)
    vwnd_r&time=date(timeind)
     if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC".or.xi_model.eq."HadGEM2-ES" then
    uwnd=sst_cold_annual(uwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2004)
    vwnd=sst_cold_annual(vwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2004)
    else 
    uwnd=sst_cold_annual(uwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    vwnd=sst_cold_annual(vwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    end if
    delete(timeind)
    delete(uwnd_r)
    delete(vwnd_r)
    delete(time)
    delete(date)

   print("******************读入海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/j/CMIP5/history/slp_chazhi/psl_"+xi_model+"_historical_r1i1p1_"+"*"+".nc")
    print(fils1)
    f1=addfile(fils1, "r")
   
    time=f1->time

    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    

      ; psl_r=f1->PSL(timeind,:,:)
      psl_r=lonPivot(f1->psl(timeind,:,:),0)
 

     print(psl_r&lon)

      psl_r=(/psl_r/)/100.0
    printVarSummary(psl_r)
    delete(psl_r&time)
    psl_r&time=date(timeind)
    
    lonind=ind(psl_r&lon.le.lonmax.and.psl_r&lon.ge.lonmin)

     if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC".or.xi_model.eq."HadGEM2-ES" then
    psl=sst_cold_annual(psl_r(:,{latmax:latmin},lonind),timeind,1979,2004)
    else
    psl=sst_cold_annual(psl_r(:,{latmax:latmin},lonind),timeind,1979,2005)
end if
    delete(timeind)
    delete(psl_r)
    delete(time)
    delete(date)
    delete(lonind)




    print("*****************计算趋势************************")  

    sst_rc1=regCoef_n(ispan(1,dimsizes(sst&time),1),sst(:,:,:),0,0)  
    copy_VarCoords(sst(0,:,:), sst_rc1)
    delete(sst)
    printVarSummary(sst_rc1)
    prob1=rc_test(sst_rc1)
    copy_VarCoords(sst_rc1, prob1)
    printVarSummary(prob1)
    sst_rc1=(/sst_rc1/)*10
    

    uwnd_rc1=regCoef_n(ispan(1,dimsizes(uwnd&time),1),uwnd(:,:,:),0,0)  
    copy_VarCoords(uwnd(0,:,:), uwnd_rc1)
    delete(uwnd)
    printVarSummary(uwnd_rc1)
    uwnd_prob1=rc_test(uwnd_rc1)
    
   

    vwnd_rc1=regCoef_n(ispan(1,dimsizes(vwnd&time),1),vwnd(:,:,:),0,0)  
    copy_VarCoords(vwnd(0,:,:), vwnd_rc1)
    delete(vwnd)
    printVarSummary(vwnd_rc1)
    vwnd_prob1=rc_test(vwnd_rc1)

    uwnd_rc1=where(uwnd_prob1.le.0.3.or.vwnd_prob1.le.0.3, uwnd_rc1, uwnd_rc1@_FillValue)
    uwnd_rc1=(/uwnd_rc1/)*10

    ; vwnd_rc1=where(uwnd_prob1.le.0.2.or.vwnd_prob1.le.0.2, vwnd_rc1, vwnd_rc1@_FillValue)
    vwnd_rc1=(/vwnd_rc1/)*10
    
    psl_rc1=regCoef_n(ispan(1,dimsizes(psl&time),1),psl(:,:,:),0,0)  
    copy_VarCoords(psl(0,:,:), psl_rc1)
    delete(psl)
    printVarSummary(psl_rc1)   
    psl_rc1=(/psl_rc1/)*10

    print("**********************画图*****************************")


    ; colormap=read_colormap_file("BlWhRe")
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)
    ; cmap(0:4,:)=colormap(0:40:10,:)
    ; cmap(5:17,:)=colormap(52:100:4,:)

    ; colormap=read_colormap_file("GMT_jet")
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)
    ; cmap(0,:)=colormap(5,:)
    ; cmap(1:4,:)=colormap(22:92:23,:)
    ; cmap(5:17,:)=colormap(106:2latmax:12,:)

     
  colormap=read_colormap_file("GMT_jet")
    printVarSummary(colormap)
    cmap=new((/20,4/),typeof(colormap))
    printVarSummary(colormap)
    cmap(0:6,:)=colormap(2:92:15,:)
    cmap(7:19,:)=colormap(106:238:11,:)


    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterString =xi_model
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
  
    res@cnFillPalette=cmap
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnFillDrawOrder="PreDraw"
   
    
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False
    
    res@mpFillDrawOrder="PostDraw"
    
    res@mpLandFillColor="gray40"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
     
    res@mpShapeMode="FreeAspect"
      res@vpWidthF=0.65
    res@vpHeightF=0.35
  


    res@gsnScalarContour        = True
    res@vcRefMagnitudeF          = 0.4 
    res@vcRefLengthF             = 0.03     ; define length of vec ref
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
   
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    
    plot = gsn_csm_vector_scalar_map(wks,uwnd_rc1,vwnd_rc1,sst_rc1,res)

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

        iplot=new(1,graphic)

        iplot=gsn_csm_contour(wks,prob1,res_pro)
     
        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 0.2         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

        iplot = gsn_contour_shade(iplot, 0,0.05, opt)
        ;;;  
        overlay(plot, iplot)

        ;;;画slp的等值线
        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -0.8        ; set min contour level
        res_slp@cnMaxLevelValF      = 0.8        ; set max contour level
        res_slp@cnLevelSpacingF     = 0.2       ; set contour spacing
        res_slp@cnLineLabelsOn     =True
        res_slp@cnLineLabelBackgroundColor="Transparent"
        res_slp@gsnContourNegLineDashPattern=2
        res_slp@cnLineLabelInterval        = 1             ; default = 2
        res_slp@cnLineLabelPlacementMode   = "constant"    ; choose constant label method
        res_slp@cnLineLabelFormat= "0@;*.1f" 
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineColor="black"
        res_slp@cnLineThicknessF=1.0
        plot_slp=gsn_csm_contour(wks,psl_rc1,res_slp)
        overlay(plot, plot_slp)

  return plot
end 

;;;********************************************************************
;;;主程序
;;;******************************************************************
begin


; xi_model=(/"ACCESS1-0","CCSM4","CSIRO-Mk3-6-0","GFDL-ESM2M"\
; ,"INMCM4","IPSL-CM5A-MR","IPSL-CM5B-LR","MPI-ESM-LR","NorESM1-M"/)


; dong_model=(/"BCC-CSM1-1","CanESM2","CMCC-CM","GFDL-CM2p1","GFDL-CM3","GFDL-ESM2G","GISS-E2-H-CC","GISS-E2-R"\
; ,"HadGEM2-AO","IPSL-CM5A-LR","MIROC5","MPI-ESM-MR","MRI-CGCM3"/)


; xi_model=(/"GFDL-CM3","GFDL-ESM2M","IPSL-CM5A-MR",\
; "IPSL-CM5B-LR","INMCM4","MPI-ESM-LR","NorESM1-M"/)

; dong_model=(/"ACCESS1-0","BCC-CSM1-1","CMCC-CM","GFDL-CM2p1","GISS-E2-R","IPSL-CM5A-LR","MIROC5","MPI-ESM-MR","MRI-CGCM3"/)


xi_model=(/"GFDL-ESM2M","FIO-ESM","GFDL-CM3","IPSL-CM5A-MR",\
"IPSL-CM5B-LR","INMCM4","NorESM1-M"/)


dong_model=(/"BCC-CSM1-1","CESM1-BGC","BCC-CSM1-1-M","MPI-ESM-MR","CMCC-CMS","HadGEM2-CC","MIROC5"/)


plot_xi=new(dimsizes(xi_model),graphic)
plot_dong=new(dimsizes(dong_model),graphic)

wks1=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/tos_group/变宽模式的sst和风场_psl趋势变化_NEW")
wks2=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/tos_group/变窄模式的sst和风场_psl趋势变化_NEW")


do i=0,dimsizes(xi_model)-1
    plot_xi(i)=draw_cmip_sst(xi_model(i),wks1)
end do

do i=0,dimsizes(dong_model)-1
    plot_dong(i)=draw_cmip_sst(dong_model(i),wks2)
end do

pres=True
pres@gsnPanelLabelBar=True
pres@lbOrientation="Vertical"

gsn_panel(wks1, plot_xi, (/5,2/), pres)
gsn_panel(wks2, plot_dong,(/5,3/), pres)

end
    
  

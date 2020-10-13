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
;;读入数据
   
 
    model="FGOALS-g2"
    diri_output="/mnt/f/share/CMIP5/msf_history/"
   ; system("if ! test -d " + diri_output +" ; then mkdir " + diri_output + " ; fi")   ;;;make dir
   ; diri_output="/mnt/f/share/CMIP5/msf_history/"+model
    yearstart=197901
    yearend=200512
    
    year1=1979
    year2=2005
    fils1=systemfunc("ls /mnt/j/CMIP5/history/tos/"+model+"/r1i1p1/tos_Omon_"+model+"_historical_r1i1p1_"+"*"+".nc")
    print(fils1)
    f1=addfiles(fils1, "join")

    time=f1[:]->time    
    date=cd_calendar(time, 1)
    print("******************日期*********************")

    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    ; print(date(timeind))

    ; sst_r=f1[:]->tos(timeind,:,:)  ;;;195801-201212
    sst_r=lonPivot(f1[:]->tos(347:670,:,:),0.5)  ;;;195801-201212
    sst_r=(/sst_r-273.15/)

    dNames = getvardims(f1[:]->tos)        ; NCL 6.4.0 and earlier
    print(dNames)
    
    if (any(dNames.eq."rlat")) then 
    sst_r!0="time"
    sst_r!1="lat"
    sst_r!2="lon"

    end if 

    sst=sst_r
    ; print("*****************开始水平插值************************")  
    ; xo=fspan(0,357.5,144)
    ; yo=fspan(-87.5, 87.5, 71)
   
    ; sst=linint2_Wrap(sst_r&lon, sst_r&lat, sst_r, \
    ;                   True, xo, yo, 0)

     sst!1="lat" 
     sst!2="lon"
    ; sst&lat@units="degrees_north"
    ; sst&lon@units="degrees_east"
    ; delete(sst_r)
    ; printVarSummary(sst)
   
    
    ; lat=sst&lat
    ; lon=sst&lon 
    ; lon@units="degree_east"
    
    ; delete(sst&time)
    ; sst&time=date(timeind)
    sst&lat@units="degrees_north"
    sst&lon@units="degrees_east"
    
     print("*****************计算趋势************************")  
    
    printVarSummary(sst)
    ;;;先得到秋季和冬季的SONDJF海温
    sst_son=month_to_season(sst(time|:,lat|:,lon|:), "SON")
    sst_djf=month_to_season(sst(time|:,lat|:,lon|:), "DJF")
    printVarSummary(sst_son)

    sst_son&time=ispan(year1,year2,1)
    sst_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    sst_djf(0:dimsizes(sst_djf&time)-2,:,:)=(/sst_djf(1:,:,:)/)
    ; sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/sst({200512},:,:)/)
    sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/sst(323,:,:)/)
    ;;;
       sst_mean_annual=new((/324/12,dimsizes(sst&lat),dimsizes(sst&lon)/),typeof(sst))
    printVarSummary(sst_mean_annual)

    ; sst_mean_annual=new((/dimsizes(timeind)/12,dimsizes(sst&lat),dimsizes(sst&lon)/),typeof(sst))
    ; printVarSummary(sst_mean_annual)

    ; sst_mean_annual=month_to_annual(sst_mean, 1)
     sst_mean_annual=(/sst_son(:,:,:)+sst_djf(:,:,:)/)*0.5
    delete(sst_son)
    delete(sst_djf)

    sst_mean_annual!0="time"
    sst_mean_annual!1="lat"
    sst_mean_annual!2="lon"
    ; sst_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    sst_mean_annual&time=ispan(0,324/12-1,1)
    sst_mean_annual&lon=sst&lon
    sst_mean_annual&lat=sst&lat
   
    printVarSummary(sst_mean_annual)

   ; sst_mean_annual=where(sst_mean_annual.eq.sst_mean_annual@_FillValue, 0, sst_mean_annual)
 

    sst_rc1=regCoef_n(ispan(1,dimsizes(sst_mean_annual&time),1),sst_mean_annual(:,{30:-30},{90:300}),0,0)  
    copy_VarCoords(sst_mean_annual(0,{30:-30},{90:300}), sst_rc1)
    printVarSummary(sst_rc1)
    ; sst_rc1=where(sst_rc1.eq.0, sst_rc1@_FillValue, sst_rc1)
   

;************************************************
; 显著性检验
;************************************************
    ; prob1=rc_test(sst_rc1)
    ; sst_rc1=(/sst_rc1/)*10
    ; printVarSummary(sst_rc1)
    ; copy_VarCoords(sst_mean_annual(1,:,:), sst_rc1)
    
    
    print("**********************海温*****************************")

 wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/"+model+"_sst趋势变化")



    colormap=read_colormap_file("BlWhRe")
    printVarSummary(colormap)

    cmap=new((/18,4/),typeof(colormap))
     
    printVarSummary(colormap)
  
    cmap(0:4,:)=colormap(0:40:10,:)
    cmap(5:17,:)=colormap(52:100:4,:)


  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False

  res@cnLineColors="Transparent"  ;
  res@lbLabelBarOn=False
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  ; res@cnFillMode="RasterFill"
  ; res@cnRasterSmoothingOn=True
  
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  
  res@cnMinLevelValF       =   -0.02
  res@cnMaxLevelValF       =   0.06
  res@cnLevelSpacingF      =  0.005
  res@cnFillPalette        =cmap
  res@tmXTOn=False
  res@tmYROn=False
  res@mpProjection="CylindricalEquidistant"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(sst_rc1&lat)
  res@mpMaxLatF             = max(sst_rc1&lat)
  res@mpMinLonF             = min(sst_rc1&lon)
  res@mpMaxLonF             = max(sst_rc1&lon)
  res@mpCenterLonF          = 180.0

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3

  
  plot=new(1,graphic)
  plot(0)= gsn_csm_contour_map(wks,sst_rc1,res)


draw(plot)
frame(wks)
end
    
  

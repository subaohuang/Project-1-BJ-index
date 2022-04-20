;;;海温场
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
    ; sst_mean_annual=(/month_to_annual(sst(time|:,lat|:,lon|:), 1)/)
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
begin

    ; fils1=systemfunc("ls /mnt/f/share/MRE/haderscobe_sst_chazhi.nc")
    fils1=systemfunc("ls /mnt/f/share/ersst/sst.mnmean.v4.nc")
      ; fils1=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc")
    
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    yearstart=197901
    yearend=201312
    ;;;全球
    latmin=-60
    latmax=60
    lonmin=0
    lonmax=360
    ;;;
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_mre_r=f1->sst(timeind,:,:)

    delete(sst_mre_r&time)
    printVarSummary(sst_mre_r)
    sst_mre_r&time =  date(timeind)
    sst_djf_r=month_to_season(sst_mre_r(time|:,lat|:,lon|:), "DJF")
    sst_mam_r = month_to_season(sst_mre_r(time|:,lat|:,lon|:), "MAM")

    sst_mre=sst_djf_r(1:,:,:)
    sst_mre=(/(sst_djf_r(1:,:,:) + sst_mam_r(1:,:,:))/2.0/)

    printVarSummary(sst_mre)
    delete(sst_mre&time)
    sst_mre!0="time"
    sst_mre!1="lat"
    sst_mre!2="lon"
    sst_mre&time=ispan(1, 34, 1)
    sst_mre&lon=sst_mre_r&lon
    sst_mre&lat=sst_mre_r&lat
    ; sst_mre=sst_cold_annual(sst_mre_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2012)
    
    sst_mre&lat@units="degree_north"
    sst_mre&lon@units="degree_east"
    delete(sst_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)
    
    sst_mre_pi=lonPivot(sst_mre, 30)
  


    print("*****************计算趋势************************")  

    sst_rc1_mre=regCoef_n(ispan(1,dimsizes(sst_mre_pi&time),1),sst_mre_pi(:,:,:),0,0)  
    copy_VarCoords(sst_mre_pi(0,:,:), sst_rc1_mre)
    printVarSummary(sst_rc1_mre)
    
    sst_rc1_mre=(/sst_rc1_mre/)*10
    prob1=rc_test(sst_rc1_mre)
    printVarSummary(prob1)
    ; exit()


;;;
    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/海温的趋势")
    colormap=read_colormap_file("MPL_RdYlBu")
    printVarSummary(colormap)
    colormap=(/colormap(::-1,:)/)
    ; colormap=(/colormap()/)
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)

    ; cmap(0:6,:)=colormap(2:92:15,:)
    ; ; cmap(0:6,:)=colormap(2:50:8,:)
    ; cmap(7:17,:)=colormap(105:255:15,:)
    ;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterString ="SST trend"
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
  
    res@cnFillPalette=colormap
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnFillDrawOrder="PreDraw"
   
    
    ; res@cnMinLevelValF       =   -0.3
    ; res@cnMaxLevelValF       =   0.5
    ; res@cnLevelSpacingF      =  0.05

    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.025
    
    res@tmXTOn=False
    res@tmYROn=False

    res@tmXBValues=(/30,90,150,210,270,330,388/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

    res@tmYLValues=(/60,30,0,-30,-60/)
    res@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N","0"\
    ,"30~S~o~N~S","60~S~o~N~S"/)


    res@tmXBLabelFontHeightF=0.015
    res@tmYLLabelFontHeightF=0.015
    res@tmXBMajorThicknessF=3.0
    res@tmYLMajorThicknessF=3.0
    res@tmBorderThicknessF=4.0


    res@mpFillDrawOrder="PostDraw"
    
    res@mpLandFillColor="white"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -60.0
    res@mpMaxLatF             = 60.0
    res@mpMinLonF             = 30
    res@mpMaxLonF             = 390
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
    
    res@mpShapeMode="FreeAspect"
    res@vpWidthF=0.7
    res@vpHeightF=0.35

    res@pmLabelBarOrthogonalPosF=0.1
  
    plot=new(1,graphic)
    plot = gsn_csm_contour_map(wks1,sst_rc1_mre,res)
    
    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(1,graphic)

    iplot(0)=gsn_csm_contour(wks1,prob1,res_pro)

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 2.0        ; add extra density
    opt@gsnShadeFillDotSizeF = 0.005      ; make dots larger

    i=0
    iplot(i) = gsn_contour_shade(iplot(i), 0,0.1, opt)

    ;;;  

    overlay(plot, iplot)
    draw(plot)
    frame(wks1)
    end 
    
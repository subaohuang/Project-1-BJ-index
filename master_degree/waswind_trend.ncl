  
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
   diri="/mnt/f/share/WASwind/"
    fils1=systemfunc("ls "+diri+"waswind_v1_0_1.monthly.nc")
    ;   diri="/mnt/f/share/ncep/ncep2/"
    ; fils1=systemfunc("ls "+diri+"uwnd.10m.mon.mean.nc")

    year1=1979
    year2=2011
    yearend=201112
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    olr_r=f1->tx(timeind,{-30:30},:)
    
    if isatt(olr_r, "add_offset") then 
      olr=olr_r@scale_factor*olr_r + olr_r@add_offset
      print("**************"+"jius"+"***")
    else 
      olr=olr_r
    end if
    copy_VarCoords(olr_r, olr)
    copy_VarAtts(olr_r, olr)
    olr!0="time"
    olr!1="lat"
    olr!2="lon"
   printVarSummary(olr)
   
    delete(olr&time)
    olr&time=date(timeind)
    printVarSummary(olr)
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_djf=month_to_season(olr(time|:,lat|:,lon|:), "DJF")
    printVarSummary(olr_djf)
    delete(olr_djf&time)
    olr_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    olr_djf(0:dimsizes(olr_djf&time)-2,:,:)=(/olr_djf(1:,:,:)/)
    olr_djf(dimsizes(olr_djf&time)-1,:,:)=(/olr({yearend},:,:)/)
    ;;;

    olr_djf1=dim_avg_n_Wrap(olr_djf({1979:1994},{-30:30},:), 0)
    olr_djf2=dim_avg_n_Wrap(olr_djf({1995:},{-30:30},:), 0)
    olr_diff=olr_djf1
    olr_diff=(/olr_djf2-olr_djf1/)

    ; olr_rc=regCoef_n(ispan(1, 33, 1), olr_djf, 0, 0)
    ; copy_VarCoords(olr_djf(1,:,:), olr_rc)
    ; prob=rc_test(olr_rc)
    ; olr_rc=(/olr_rc*10.0/)
   wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/waswind比较")
   
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    ; res@gsnLeftString=LeftString
    res@cnFillPalette="BlueRed"
    res@cnFillOn=True
    res@cnLinesOn=False
    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -0.8
    ; res@cnMaxLevelValF       =   0.8
    ; res@cnLevelSpacingF      =  0.2
    res@lbLabelBarOn=False
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             =0
    res@mpMaxLonF             = 360
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpShapeMode="FreeAspect"
    res@vpHeightF=0.3
    res@vpWidthF=0.6
    plot=gsn_csm_contour_map(wks, olr_diff, res)
    

    ; ;;;打点
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.05      ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False

    ; iplot=new(1,graphic)

    ; iplot(0)=gsn_csm_contour(wks,prob,res_pro)
    ; delete(prob)

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; opt@gsnShadeMid     = 17         ; stipple
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.8         ; add extra density
    ; opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger


    ; iplot= gsn_contour_shade(iplot, 0,0.05, opt)

    ; ;;;  

    ; overlay(plot, iplot)
    draw(plot)
    frame(wks)
    end 
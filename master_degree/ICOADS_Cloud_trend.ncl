    
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
   diri="/mnt/f/share/ICOADS/"
   fils1=systemfunc("ls "+diri+"cldc.mean.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    year1=1979
    year2=2012
    yearend=201212
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    olr_r=f1->cldc(timeind,:,:)
    olr=(/olr_r@add_offset+olr_r@scale_factor*olr_r/)
    copy_VarCoords(olr_r, olr)
    copy_VarAtts(olr_r, olr)
    
    ; olr_clm=dim_avg_n_Wrap(olr, 0)
    olr_clm=olr(1,:,:)
   
    olr_pi=lonPivot(olr, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_djf_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "DJF")
    olr_djf=olr_djf_r(1:,:,:)

    ;;;计算趋势
    olr_rc=regCoef_n(ispan(1, 34, 1), olr_djf, 0, 0)
    copy_VarCoords(olr_djf(1,:,:), olr_rc)
    prob=rc_test(olr_rc)
    olr_rc=(/olr_rc*10.0/)


    olr_time=dim_avg_n_Wrap(olr_djf(:,{10:-10},{80:150}), (/1,2/))
   
   olr_time=(/runave(olr_time, 3, 0)/)
   olr_line_rc=regCoef_n(ispan(1, 34,1), olr_time,0,0)
   olr_line=olr_line_rc*ispan(1, 34,1)+olr_line_rc@yintercept







    wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/ICOADS的云量的趋势")
    wks2=gsn_open_wks("pdf", "/mnt/f/share/era/data/ICOADS的云量时间序列")

    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillPalette="NCV_blu_red"
    res@cnFillOn=True
    res@cnLinesOn=False
    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -8.0
    ; res@cnMaxLevelValF       =   8.0
    ; res@cnLevelSpacingF      =  1.0
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbOrientation="Vertical"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             =0
    res@mpMaxLonF             = 360
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpShapeMode="FreeAspect"
    res@vpHeightF=0.3
    res@vpWidthF=0.7
    plot=gsn_csm_contour_map(wks, olr_clm, res)


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

    iplot(0)=gsn_csm_contour(wks,prob,res_pro)
 


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.6         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

    i=0
    iplot(i)= gsn_contour_shade(iplot(i), 0,0.05,opt)

     
draw(plot)
frame(wks)
;;;
resxy=True

plot2=gsn_csm_xy(wks2, ispan(1, 34, 1), (/olr_time,olr_line/), resxy)







  end 
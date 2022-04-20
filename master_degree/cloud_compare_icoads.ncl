    ;;;Icoads 云量
    begin
    f4=addfile("/mnt/f/share/ICOADS/cldc.mean.nc", "r")
    date=cd_calendar(f4->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    cldc_r=f4->cldc(timeind,:,:)
    cldc=cldc_r*cldc_r@scale_factor+cldc_r@add_offset
    copy_VarCoords(cldc_r, cldc)
    delete(cldc&time)
    cldc&time=date(timeind)
    printVarSummary(cldc)
    
    ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
    opt=(/0,5/)
    ;;先插值lon方向的
    cldc=(/linmsg_n(cldc, opt, 2)/)
    ;;插值lat方向的
    opt=(/0,1/)
    cldc=(/linmsg_n(cldc, opt, 1)/)
    delete(opt)

    cldc_pi=lonPivot(cldc, 30)
    cldc_djf_r=month_to_season(cldc_pi(time|:,lat|:,lon|:), "DJF")
    ;;注意点
    cldc_djf=cldc_djf_r(0:,:,:)
    delete(cldc_djf&time)
    cldc_djf&time=ispan(1979,2012,1)
    ; cldc_djf_mean=dim_avg_n_Wrap(cldc_djf(:,{30:-30},:), (/1,2/))
    ; do i=0,dimsizes(cldc_djf&time)-1
    ;     cldc_djf(i,:,:)=(/cldc_djf(i,:,:)-cldc_djf_mean(i)/)
    ; end do
    ;;;
    cldc_djf1=dim_avg_n_Wrap(cldc_djf({1979:1994},:,:), 0)
    cldc_djf2=dim_avg_n_Wrap(cldc_djf({1995:2012},:,:), 0)
    cldc_diff=cldc_djf1
    cldc_diff=(/cldc_djf2- cldc_djf1/)


  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/icoads云量的时间平均比较")


  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftString="(a) 1979-1994 mean" 
  res@cnFillPalette="NCV_blu_red"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  .1
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  ; res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             =30
  res@mpMaxLonF             = 390
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  res@mpShapeMode="FreeAspect"
  res@vpHeightF=0.3
  res@vpWidthF=0.7  
  plot=new(2,graphic)
  plot(0)=gsn_csm_contour_map(wks, cldc_diff, res)
  res@gsnLeftString="(a) 1995-2012 mean" 
  plot(1)=gsn_csm_contour_map(wks, cldc_djf2, res)
draw(plot(0))
frame(wks)
  ; pres=True
  ; pres@gsnPanelLabelBar=True
  ; gsn_panel(wks, plot, (/2,1/), pres)

end 



  

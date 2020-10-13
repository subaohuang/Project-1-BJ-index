;得到的各相关系数的滤波
;*****************************************************
;滤波
;****************************************************
function clm_ano(x[*]:numeric)  ;去除一维时间序列的年循环
begin

   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check; Compute all 12 monthly averages first. 
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do
   copy_VarAtts (x,aveMonth)
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   xAnom=x
   do yr=0,ntim-1,nmos
      xAnom(yr:yr+nmos-1) = (/ x(yr:yr+nmos-1)- aveMonth /)
   end do
   return (xAnom)
end
;*************************************************************
begin
diri="/mnt/f/share/era/data/"
diri_pi="/mnt/f/share/era/"
fils1=systemfunc("ls "+diri+"msf_era_level_location1979-2017.nc")
f1=addfile(fils1, "r")
locationxi_r=f1->locationxi_a(0:455)
locationxi_ano=clm_ano(locationxi_r)





;**********************************************************************
;Lanczos滤波
;*********************************************************************


  ihp   = 2                            ; band_pass
  sigma = 1.0                           ; Lanczos sigma
  nWgt  = 121                           ; loose 60 months each end
  fcb   = 1./24.                        ; 2 years
  fca   = 1./84.                         ;7年
  wgts  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )

   locationxi_band=wgt_runave_Wrap(locationxi_ano, wgts, 0)
   locationxi_low=locationxi_ano- locationxi_band

    copy_VarCoords(locationxi_ano, locationxi_low)
   

  wks=gsn_open_wks("png", ""+diri_pi+"上升支位置滤波")

  res              = True            ; plot mods desired
  res@gsnDraw      = False           ; don't draw 
  res@gsnFrame     = False           ; don't advance frame yet      

  res@vpHeightF    = 0.4             ; change aspect ratio of plot
  res@vpWidthF     = 0.8                  
  
  ; res@trYMinF      = -3.0            ; min value on y-axis
  ; res@trYMaxF      =  3.0            ; max value on y-axis
  ; res@trXMinF      = 1860            ; min value on x-axis
  ; res@trXMaxF      = 2016            ; max value on x-axis
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(1,456,60)
    res@tmXBLabels=ispan(1979,2016,5)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(1,456,12)
    
  res@vpXF         = 0.1             ; start plot at x ndc coord 
  res@gsnYRefLine  = 0.0             ; create a reference line   
       
  res@xyMonoDashPattern= True
  res@xyLineThicknessF = 2
  res@xyLineColors     = (/ "blue", "red","green","black","yellow"/)
  ; res@gsnCenterString  = "Low Pass Filtered: 2-year and 10-year"     
  yplot=new((/3,dimsizes(locationxi_ano)/),typeof(locationxi_ano))

  yplot(0,:)=locationxi_ano
  yplot(1,:)=locationxi_band
  yplot(2,:)=locationxi_low
  
 time=ispan(61,396,1)

  rc1=regCoef(locationxi_low&time,locationxi_low)
  rc2=regCoef(locationxi_band&time,locationxi_band)
   y1 = rc1*time+rc1@yintercept
   y2 = rc2*time+rc2@yintercept

  plot= gsn_csm_xy (wks,locationxi_low&time,yplot,res)     
  
  ployres=True
  ployres1=True
  ployres@xyLineColor="black"
  ployres1@xyLineColor="black"
  plot1=gsn_add_polyline(wks, plot, time, y1, ployres)
  plot2=gsn_add_polyline(wks, plot, time, y2, ployres1)



  draw(plot)
  frame(wks)




end

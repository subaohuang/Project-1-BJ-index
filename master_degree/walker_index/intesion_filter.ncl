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
fils1=systemfunc("ls "+diri+"era_intension.nc")
f1=addfile(fils1, "r")
msf_500=f1->msf_500
STRF=f1->STRF
slp=f1->slp
ud=f1->ud
u10=f1->u10

msf_500_ano=clm_ano(msf_500)
STRF_ano=clm_ano(STRF)
slp_ano=clm_ano(slp)
ud_ano=clm_ano(ud)
u10_ano=clm_ano(u10)
;**********************************************************************
;Lanczos滤波
;*********************************************************************


  ihp   = 2                            ; band_pass
  sigma = 1.0                           ; Lanczos sigma
  nWgt  = 121                           ; loose 60 months each end
  fcb   = 1./24.                        ; 2 years
  fca   = 1./84.                         ;7年
  wgts  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )

   msf_500_band=wgt_runave_Wrap(msf_500_ano, wgts, 0)
   msf_low=msf_500_ano-msf_500_band

   copy_VarCoords(msf_500, msf_low)
   

 wks=gsn_open_wks("png", ""+diri+"滤波msf_500")

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
  res@xyLineColors     = (/ "blue", "red","green","black"/)
  res@gsnCenterString  = "Low Pass Filtered: 2-year and 10-year"     
  yplot=new((/4,dimsizes(msf_500_ano)/),typeof(msf_500_ano))

  yplot(0,:)=msf_500_ano
  yplot(1,:)=msf_500_band
  yplot(2,:)=msf_low

  rc1=regCoef(msf_low&time,msf_low)
  yplot(3,:) = rc1*msf_low&time+rc1@yintercept
  
  plot= gsn_csm_xy (wks,msf_500&time,yplot,res)     

printVarSummary(msf_low)
printVarSummary(msf_500)
  
draw(plot)
frame(wks)
printMinMax(msf_low, 1)
printMinMax(msf_500, 1)



end

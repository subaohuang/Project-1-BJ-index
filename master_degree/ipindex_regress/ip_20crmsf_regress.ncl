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
;;;
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
;;;
begin
ipindex=asciiread("/mnt/f/share/era/data/ip_index_1950_2010.ascii", 732, "float")
ipindex_de=dim_standardize_n(dtrend(clm_ano(ipindex),True),0,0)
ipindex_de!0="time"
ipindex_de&time=ispan(1,732,1)



fils=systemfunc("ls /mnt/f/share/ncep/20CR/msf_time_20cr_global_10_-10.nc") 
f1=addfile(fils, "r")


 level=f1->level(:)
 lon=f1->lon(:)

 msf_time=f1->msf_time(:,:,0:731)
 msf_time=(/msf_time/10^9/)
 ud_time=f1->ud_time(:,:,0:731)
 w_time=f1->w_time(:,:,0:731)
 msf_time&time=ispan(1, 732, 1)
 copy_VarCoords(msf_time, ud_time)
 copy_VarCoords(msf_time, w_time)
 


 msfa=dtrend_msg_n(ispan(0,dimsizes(msf_time&time)-1,1),rmMonAnnCycLLT(msf_time),\
                     False,False,2)
  copy_VarCoords(msf_time, msfa)

 uda=dtrend_msg_n(ispan(0,dimsizes(ud_time&time)-1,1),rmMonAnnCycLLT(ud_time),\
                     False,False,2)
  copy_VarCoords(ud_time, uda)

   wa=dtrend_msg_n(ispan(0,dimsizes(w_time&time)-1,1),rmMonAnnCycLLT(w_time),\
                     False,False,2)
  copy_VarCoords(w_time, wa)

 printVarSummary(uda) 

;;;回归p场

  msfa_rc=regCoef_n(ipindex_de,msfa,0,2)  
  copy_VarCoords(msfa(:,:,0), msfa_rc)
  
  uda_rc=regCoef_n(ipindex_de,uda,0,2)  
  copy_VarCoords(uda(:,:,0), uda_rc)

  wa_rc=regCoef_n(ipindex_de,wa,0,2)  
  copy_VarCoords(wa(:,:,0), wa_rc)

print("回归结束")



;;;


;;;显著性检验
  prob_msfa=rc_test(msfa_rc)
  prob_uda=rc_test(uda_rc)
  prob_wa=rc_test(wa_rc)
  
;;;
print("检验结束")

;;;风场显著性

  uda_rc=where(prob_uda.le.0.05, uda_rc, uda_rc@_FillValue)
  wa_rc=where(prob_wa.le.0.05, wa_rc, wa_rc@_FillValue)

  uda_rc=(/uda_rc*100/)
  wa_rc=(/wa_rc*100/)

;;;
print("风场结束")


;;;画图 
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/20crmsf与IP指数回归")

  gsn_define_colormap(wks, "BlueWhiteOrangeRed")
  ;;;画图
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    
    ; res@lbLabelStrings=ispan(-6,6,1)
    ; res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    ; res@cnLevelSelectionMode ="ManualLevels"
    ; res@cnMinLevelValF = -6.0
    ; res@cnMaxLevelValF = 6.0
    ; res@cnLevelSpacingF = 0.5
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
    
    res@cnLineLabelsOn=True
    res@cnLineLabelBackgroundColor="Transparent"
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False
      
    res@lbOrientation="Horizontal"
 
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
    res@tmYRMode="Automatic"

    ; res@vcRefMagnitudeF = 2.0                ; define vector ref mag
    ; res@vcRefLengthF    = 0.04        ; define length of vec ref
    ; res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    ; res@vcMinDistanceF  = 0.02       ; thin out vectors
    ; res@vcMapDirection  = False
    
    ; ; res@vcFillArrowHeadXF=0.1
    ; res@vcFillArrowMinFracWidthF=1.0
    ; res@vcFillArrowHeadMinFracYF=1.0
    ; ; res@vcFillArrowHeadXF=1.0
    ; res@vcFillArrowWidthF=0.05
    ; res@vcFillArrowHeadYF=0.15

    ; res@vcRefAnnoOrthogonalPosF  = -1.15
    ; res@vcRefAnnoPerimOn=False
    
    res@vpHeightF=0.3
    res@vpWidthF=0.7
    
    ; plot=gsn_csm_pres_hgt_vector(wks,msfa_rc,uda_rc,wa_rc,res)
    plot=gsn_csm_pres_hgt(wks, msfa_rc, res)
    print("画图结束")

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
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True
 
  

  iplot=gsn_csm_pres_hgt(wks,prob_msfa,res_pro)
 

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
  iplot = gsn_contour_shade(iplot, 0,0.05, opt)


;;---Resources for the polyline  添加下垫面的陆地
  ployres                  = True          ; polyline mods desired
  ployres@gsLineThicknessF = 15.0            ; line thickness   
  ployres@gsLineColor      = "black"         ; color of lines
  ployres@gsLineDashPattern=15
 ;---Attach the polyline
  ploylevel=1000
  ploylon1=fspan(9, 39, 25)
  ploylon2=fspan(100, 149, 25)
  ploylon3=fspan(280, 315, 25)
  ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
  ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
  ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)

  ployres@gsLineDashPattern=0
  
  dum1= gsn_add_polyline(wks, plot,  ploylon1,ploylevel1,  ployres)

  ployres@gsLineDashPattern=16
  dum2= gsn_add_polyline(wks, plot,  ploylon2,ploylevel2,  ployres)
  ployres@gsLineDashPattern=0
  dum3= gsn_add_polyline(wks, plot,  ploylon3,ploylevel3,  ployres)




;;;  

overlay(plot, iplot)
draw(plot)
frame(wks)



end 
function one_test(pv_ano,east_time,pv_clm)
begin

    pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)

    t_pv_ano_east = (pv_east-pv_clm)/var_pv_ano_east

    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)

    confi_pv_east_se=confi_pv_east(::4,::4)

 return confi_pv_east
 end
;************************************************************
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
;*********************************************************************
;画上升支位置的bar图
begin

; diri="/mnt/f/share/era/data/shixiong/"
; diri3=diri+"locationxi.ascii"
; diri4=diri+"locationdong.ascii"
; locationxi=asciiread(diri3, 432, "float") 
; locationdong=asciiread(diri4, 432, "float")

 diri_ceshi="/mnt/f/share/era/data/"
 fils0=systemfunc("ls "+diri_ceshi+"msf_era_level_location1979-2017.nc")
 f0=addfile(fils0, "r")
 locationxi_r=f0->locationxi_a(0:431)
 locationxi=clm_ano(locationxi_r)

diri="/mnt/f/share/era/data/shixiong/"

locationxi_anu=month_to_annual(locationxi, 1)
locationxi_anu_std=dim_standardize(locationxi_anu, 1)
locationxi_anu_std!0="time"
printVarSummary(locationxi_anu_std)
locationxi_anu_std&time=ispan(1,36,1)
wks4=gsn_open_wks("pdf", ""+diri+"WL年变化的bar图")
;;;提取数据
   
 east_time = ind(locationxi_anu_std.gt.0.5)
 west_time=ind(locationxi_anu_std.lt.0.5)
 climate_mean_state=ind(locationxi_anu_std.ge.-0.5.and.locationxi_anu_std.le.0.5)

fils=systemfunc("ls "+diri+"era_iod_pacific_w_u.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)
w_r=f1->w(tind,{100:1000},{10:-10},:)
u_r=f1->u(tind,{100:1000},{10:-10},:)
w_month=(w_r@scale_factor*w_r+w_r@add_offset)*(100.0)
u_month=(u_r@scale_factor*u_r+u_r@add_offset)/10.0
copy_VarCoords(w_r, w_month)
copy_VarCoords(u_r, u_month)

w=month_to_annual(w_month, 1)
w!0="time"
w&time=ispan(1, 36, 1)


u=month_to_annual(u_month, 1)
u!0="time"
u&time=ispan(1, 36, 1)

w_latavg=dim_avg_n_Wrap(w, 2)
u_latavg=dim_avg_n_Wrap(u, 2)
w_latavg=-1*w_latavg
u_latavg=-1*u_latavg

w_east= dim_avg_n_Wrap(w(east_time,:,:,:),2)
w_west= dim_avg_n_Wrap(w(west_time,:,:,:),2)
u_east= dim_avg_n_Wrap(u(east_time,:,:,:),2)
u_west= dim_avg_n_Wrap(u(west_time,:,:,:),2)

w_clm=dim_avg_n_Wrap(w,2)
u_clm=dim_avg_n_Wrap(u, 2)

w_east_avg=dim_avg_n_Wrap(w_east(level|:,longitude|:,time|:),2)
u_east_avg=dim_avg_n_Wrap(u_east(level|:,longitude|:,time|:),2)

w_west_avg=dim_avg_n_Wrap(w_west(level|:,longitude|:,time|:), 2)
u_west_avg=dim_avg_n_Wrap(u_west(level|:,longitude|:,time|:), 2)

w_clm_avg=dim_avg_n_Wrap(w_clm(level|:,longitude|:,time|:), 2)
u_clm_avg=dim_avg_n_Wrap(u_clm(level|:,longitude|:,time|:), 2)

w_east_avg_ne=-1*w_east_avg
w_west_avg_ne=-1*w_west_avg
w_clm_avg_ne=-1*w_clm_avg
copy_VarCoords(w_east_avg,w_east_avg_ne)
copy_VarCoords(w_west_avg,w_west_avg_ne)
copy_VarCoords(w_clm_avg,w_clm_avg_ne)

;;;w的t检验
    ; confi_w_east = one_test(w_latavg,east_time,w_clm_avg_ne)
    ; confi_w_west = one_test(w_latavg,west_time,w_clm_avg_ne)

    ; printVarSummary(confi_w_east)
    ; ; confi_w_east_se=confi_w_east(::2,::12)


    ; prob1d=ndtooned(confi_w_east)
    ; i1d=ind(prob1d.lt.0.05)
    ; i=ind_resolve(i1d, dimsizes(confi_w_east))

    ; ilat=confi_w_east&level(i(:,0))
    ; ilon=confi_w_east&longitude(i(:,1))


;;;w的t检验
  confi_w_east = one_test(w_latavg,east_time,w_clm_avg_ne)
  confi_w_west = one_test(w_latavg,west_time,w_clm_avg_ne)
; 




wks=gsn_open_wks("pdf", ""+diri+"流场图")
;;;画bar图

    resxy=True
    resxy@gsnLeftString="(c)"
    resxy@gsnXYBarChart = True            ; Create bar plot
    resxy@tmYLMode="Explicit"
    resxy@trYMaxF=2.0
    resxy@trYMinF=-2.0
    
    resxy@tmYLValues=fspan(-2, 2, 9)
    resxy@tmYLLabels=fspan(-2, 2, 9)
    resxy@gsnYRefLine  = 0.   
    resxy@gsnBelowYRefLineBarColors="blue"
    resxy@gsnAboveYRefLineBarColors="red"
    resxy@tiYAxisString="WL(index)"
    resxy@tmYLMinorOn=False
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(locationxi_anu_std&time)
    resxy@tmXBValues=ispan(1,36,5)
    resxy@tmXBLabels=ispan(1980,2015,5)
    resxy@tmXBMinorValues=ispan(1,36,1)
    
    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9  
    plot = gsn_csm_xy (wks4,locationxi_anu_std&time , locationxi_anu_std,resxy)
;;--------------------------------------------------------------------------------
;;;画w和流场图
 gsn_define_colormap(wks,"hotcolr_19lev")
 ; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillOn          = True                ; turn on color fill
    res@cnLinesOn=False
    ; res@cnFillPalette     = "BlAqGrYeOrReVi200" ; choose color map
    res@cnLineLabelsOn    = False               ; turn off line labels
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -10
    res@cnMaxLevelValF = 10
    res@cnLevelSpacingF = 2
    res@lbLabelStride     = 2                   ; label every other box
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(30,290,30)
    res@tmXBLabels=(/"30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
                    "180~S~o~N~E","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W"/)
    res@tmXBMinorOn=False

    res@tiYAxisString="pressure(hpa)"
    res@stMinArrowSpacingF = 0.008             ; arrow spacing.
    ; res@stArrowLengthF     = 0.008             ; arrow length
    
    res@stStepSizeF= 0.001   ;ontrols the basic step size used to create the streamline

    res@vpHeightF=0.26
    res@vpWidthF=0.7
    res@lbLabelBarOn=False
    res@cnInfoLabelOn=False
    res@gsnLeftString="(d)"
    res@gsnRightString="East phase"

delete(plot)
;---draw plot from pole to pole at 170E
plot=new(3,graphic)
 plot(0)  = gsn_csm_pres_hgt_streamline(wks,w_east_avg_ne,u_east_avg,w_east_avg,res )  
 delete(res@gsnLeftString)
 res@gsnRightString="climate_mean_state "
 plot(1)  = gsn_csm_pres_hgt_streamline(wks,w_clm_avg_ne,u_clm_avg,w_clm_avg,res )   
  res@gsnRightString="West phase"
 plot(2)  = gsn_csm_pres_hgt_streamline(wks,w_west_avg_ne,u_west_avg,w_west_avg,res ) 


   ; polyres                   = True
   ; polyres@gsMarkerIndex     = 16          ; polymarker style
   ; polyres@gsMarkerSizeF     = 3.         ; polymarker size
   ; polyres@gsMarkerColor   ="black"
   ;  dum1=gsn_add_polymarker(wks, plot(0), ilon(:), ilat(:), polyres)
  
;;过检验的shade图
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

    iplot=gsn_csm_pres_hgt(wks,confi_w_east,res_pro)
 
    iplot1=gsn_csm_pres_hgt(wks,confi_w_west,res_pro)

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    opt@gsnShadeFillDotSizeF = 3.0     ; make dots larger
    iplot = gsn_contour_shade(iplot, 0,0.05, opt)
    iplot1 = gsn_contour_shade(iplot1, 0,0.05, opt)
    
    overlay(plot(0), iplot)
    overlay(plot(2), iplot1)
;;;panel
    pres=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.5
    gsn_panel(wks, plot, (/3,1/), pres)

  end

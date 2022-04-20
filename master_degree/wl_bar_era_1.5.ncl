function one_test(pv_ano,east_time,pv_clm)
begin
 n = dimsizes(east_time)
    pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)

    t_pv_ano_east = (pv_east-pv_clm)/var_pv_ano_east

    confi_pv_east = student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)

    confi_pv_east_se=confi_pv_east(::4,::4)

 return confi_pv_east
 end
;************************************************************
;************************************************************
;画上升支位置的bar图
begin

; diri="/mnt/f/share/era/data/shixiong/"
; diri3=diri+"locationxi.ascii"
; diri4=diri+"locationdong.ascii"
; locationxi=asciiread(diri3, 432, "float") 
; locationdong=asciiread(diri4, 432, "float")

;;************************************
;;read data
;**************************************
locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")
locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)

locationxi!0 = "time"
locationxi&time = ispan(1,456,1)
locationxi_anu = month_to_annual(locationxi, 1)
 
; ;;注意点
;locationxi_anu = (/month_to_season(locationxi, "DJF")/)
; locationxi_anu(0:36) = locationxi_anu(1:)
; locationxi_anu(37) = locationxi(455)
; ;;;


locationxi_anu_std=dim_standardize(locationxi_anu, 1)
locationxi_anu_std!0="time"
printVarSummary(locationxi_anu_std)
locationxi_anu_std&time=ispan(1,38,1)
;;************************************
;;相位数据
;**************************************
   
 east_time = ind(locationxi_anu_std.lt.-1)
 west_time=ind(locationxi_anu_std.gt.1)
 climate_mean_state=ind(locationxi_anu_std.ge.-0.5.and.locationxi_anu_std.le.0.5)
;;;
;;************************************
;;read  w data
;**************************************
    fils=systemfunc("ls /mnt/f/share/era/data/data_series/era_pre_global.nc")
    f1=addfile(fils, "r")
    time=f1->time
    date=cd_calendar(time, 1)
    tind=ind(date.le.201612.and.date.ge.197901)
    w_r=f1->w(tind,{100:1000},{10:-10},:)
    u_r=f1->u(tind,{100:1000},{10:-10},:)
    w_month=(w_r@scale_factor*w_r+w_r@add_offset)*(100.0)
    u_month=(u_r@scale_factor*u_r+u_r@add_offset)/10.0
    copy_VarCoords(w_r, w_month)
    copy_VarCoords(u_r, u_month)

    w=month_to_annual(w_month, 1)
    w!0="time"
    w&time=ispan(1, 38, 1)


    u=month_to_annual(u_month, 1)
    u!0="time"
    u&time=ispan(1, 38, 1)

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
;;;



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
;;;


    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/流场图")
    wks4=gsn_open_wks("pdf", "/mnt/f/share/era/data/WL年变化的bar图")
    ; wks4=gsn_open_wks("pdf", "/mnt/f/share/era/data/WL冬季变化的bar图")

;;;画bar图

    resxy=True
    resxy@gsnDraw = False
    resxy@gsnFrame = False
    resxy@gsnLeftString="(c)"
    resxy@gsnXYBarChart = False            ; Create bar plot
    resxy@tmYLMode="Explicit"
    resxy@trYMaxF=2.0
    resxy@trYMinF=-2.0
    resxy@trXMinF =1
    resxy@trXMaxF = 38


    resxy@tmYLValues=fspan(-2, 2, 9)
    resxy@tmYLLabels=fspan(-2, 2, 9)
    
    resxy@gsnYRefLine = 0.0
    ; resxy@gsnBelowYRefLineBarColors="blue"
    ; resxy@gsnAboveYRefLineBarColors="red"

    resxy@gsnBelowYRefLineColor="blue"
    resxy@gsnAboveYRefLineColor="red"

    resxy@tiYAxisString="WL(index)"
    resxy@tmYLMinorOn=False
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(locationxi_anu_std&time)
    resxy@tmXBValues=ispan(2,38,5)
    resxy@tmXBLabels=ispan(1980,2016,5)
    resxy@tmXBMinorValues=ispan(1,38,1)
    
    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9  
    plot = gsn_csm_xy (wks4,locationxi_anu_std&time , locationxi_anu_std,resxy)
    
    lineres =True
    lineres@gsLineDashPattern = 1
    
    y1 = conform(ispan(1,38,1), -1, 0)
    y3 = conform(ispan(1,38,1), 1, 0)
    dum1 = gsn_add_polyline(wks4, plot, ispan(1,38,1), y1, lineres)
    dum3 = gsn_add_polyline(wks4, plot, ispan(1,38,1), y3, lineres)
    draw(plot)
    frame(wks4)
    delete(plot)
;;;;
    
;;--------------------------------------------------------------------------------
;;;画w和流场图
   ; gsn_define_colormap(wks,"hotcolr_19lev")
 ; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillOn          = True                ; turn on color fill
    res@cnLinesOn=False
    res@cnFillPalette     = "ViBlGrWhYeOrRe" ; choose color map
    res@cnLinesOn = True
    res@cnLineLabelsOn    = False               ; turn off line labels

    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -5
    res@cnMaxLevelValF = 5
    res@cnLevelSpacingF = 1
    
    ; res@cnLevelSelectionMode ="ExplicitLevels"
    ; res@cnLevels = (/-5,-4,-3,-2,2,3,4,5/)

    res@lbLabelStride     = 1                   ; label every other box
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
    res@tmXBMinorOn=False

    res@tiYAxisString="pressure(hpa)"
  
   
    res@lbLabelBarOn =False
    
    res@vpHeightF=0.26
    res@vpWidthF=0.7
    res@lbLabelBarOn=False
    res@cnInfoLabelOn=False
    res@gsnLeftString="(a)"
    res@gsnRightString="west phase"
 

    ;---draw plot 
    plot=new(3,graphic)
    plot(0)  = gsn_csm_pres_hgt(wks,w_east_avg_ne,res)  
    res@gsnLeftString="(b)"
    res@gsnRightString="East phase"
    plot(1)  = gsn_csm_pres_hgt(wks,w_west_avg_ne,res )
    res@gsnLeftString="(c)"
    res@gsnRightString="climate_mean_state "
    plot(2)  = gsn_csm_pres_hgt(wks,w_clm_avg_ne,res) 
   
    vcres =res
    vcres@gsnLeftString = ""
    vcres@gsnRightString = ""
    vcres@vcRefMagnitudeF         = 3.0             ; define vector ref magnitude
    vcres@vcRefLengthF            = 0.05             ; define length of vec ref
    vcres@vcMinFracLengthF        = 0.2               ;0.0完全按长度，1.0长度全都一样
    vcres@vcMinDistanceF          = 0.025             ;the distance between the adjacent vector 

    vcres@vcLineArrowThicknessF   =1.5
    vcres@vcGlyphStyle =        "CurlyVector"         ;矢量弯曲
   
    vcres@vcLevelColors= (/ "gray53","black"/)
    ; vcres@vcLevelColors= (/ "green","black"/)
    vcres@vcMonoLineArrowColor    = False
    vcres@vcLevelSelectionMode= "ExplicitLevels"
    vcres@vcLevels= (/0.05/)

    vcres@vcRefAnnoOrthogonalPosF = -1.0            ;设置参考矢量高度，越负越高
    vcres@vcRefAnnoFontHeightF    = 0.01            ;设置参考矢量大小
    vcres@vcRefAnnoPerimOn      =False

    confi_clm = w_clm_avg_ne
    confi_clm = 0.06
    plot_vc=new(3,graphic)
    plot_vc(0)=gsn_csm_vector_scalar(wks,u_east_avg,w_east_avg,confi_w_east,vcres)
    plot_vc(1)=gsn_csm_vector_scalar(wks,u_west_avg,w_west_avg,confi_w_west,vcres)
    plot_vc(2)=gsn_csm_vector_scalar(wks,u_clm_avg,w_clm_avg,confi_clm, vcres)
    
    overlay(plot(0), plot_vc(0))
    overlay(plot(1), plot_vc(1))
    overlay(plot(2), plot_vc(2))
;;;;
    

;;;panel
    pres=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.5
    gsn_panel(wks, plot, (/3,1/), pres)

  end

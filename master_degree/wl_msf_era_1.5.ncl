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
    
 return confi_pv_east
 end
;************************************************************
;*********************************************************************
;画上升支位置的bar图
begin

; diri="/mnt/f/share/era/data/shixiong/"
; diri3=diri+"locationxi.ascii"
; diri4=diri+"locationdong.ascii"
; locationxi=asciiread(diri3, 432, "float") 
; locationdong=asciiread(diri4, 432, "float")


locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")


locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)

locationxi_anu=month_to_annual(locationxi, 1)



locationxi_anu_std=dim_standardize(locationxi_anu, 1)
locationxi_anu_std!0="time"
printVarSummary(locationxi_anu_std)
locationxi_anu_std&time=ispan(1,38,1)
;;;提取数据
   
 east_time = ind(locationxi_anu_std.lt.-1)
 west_time=ind(locationxi_anu_std.gt.1)
 climate_mean_state=ind(locationxi_anu_std.ge.-0.5.and.locationxi_anu_std.le.0.5)

fils = systemfunc("ls /mnt/f/share/rea_msf/msf_time_era_global.nc")
f1 = addfile(fils, "r")
time = f1->time
date = cd_calendar(time, 1)
timeind = ind(date.le.201612.and.date.ge.197901)
w_r = f1->w_time({100:1000},:,timeind)
u_r = f1->ud_time({100:1000},:,timeind)
msf_r = f1->msf_time({100:1000},:,timeind)
msf_r = (/msf_r/10^9/)


w = month_to_annual(w_r(time|:,level|:,lon|:), 1)
w!0 ="time"
w&time = ispan(1, 38, 1)


u = month_to_annual(u_r(time|:,level|:,lon|:), 1)
u!0 = "time"
u&time = ispan(1, 38, 1)

msf_time = month_to_annual(msf_r(time|:,level|:,lon|:), 1)
msf_time!0 = "time"
msf_time&time = ispan(1, 38, 1)


w_east = w(east_time,:,:)
w_west = w(west_time,:,:)

u_east = u(east_time,:,:)
u_west = u(west_time,:,:)

msf_east = msf_time(east_time,:,:)
msf_west = msf_time(west_time,:,:)


w_clm = dim_avg_n_Wrap(w,0)
u_clm = dim_avg_n_Wrap(u, 0)
msf_clm = dim_avg_n_Wrap(msf_time, 0)

w_east_avg=dim_avg_n_Wrap(w_east,0)
u_east_avg=dim_avg_n_Wrap(u_east,0)

w_west_avg=dim_avg_n_Wrap(w_west, 0)
u_west_avg=dim_avg_n_Wrap(u_west, 0)

msf_west_avg=dim_avg_n_Wrap(msf_west, 0)
msf_east_avg=dim_avg_n_Wrap(msf_east, 0)



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
  confi_w_east = one_test(w,east_time,w_clm)
  confi_w_west = one_test(w,west_time,w_clm)

; 
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/WL指数定义的msf图")

;;--------------------------------------------------------------------------------
;;;画w和流场图
   gsn_define_colormap(wks,"hotcolr_19lev")
 ; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillOn          = True                ; turn on color fill
    res@cnLinesOn=False
    res@cnLineLabelsOn    = False               ; turn off line labels
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1
    res@cnInfoLabelOn=False

    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
    

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)

    res@tmXBMinorOn=False

    res@tiYAxisString="pressure(hpa)"
    res@gsnRightString="West phase"
    res@gsnLeftString="(a)"

    res@lbLabelBarOn =False
    res@lbLabelStride     = 1                   ; label every other box
    res@lbLabelBarOn=False
  
   
    res@vpHeightF=0.4
    res@vpWidthF=0.7

    ;---draw plot 
        
  
    plot=new(3,graphic)
    plot(0)  = gsn_csm_pres_hgt(wks,msf_east_avg,res)  
     res@gsnLeftString="(b)"
    res@gsnRightString="East phase"
    
    plot(1)  = gsn_csm_pres_hgt(wks,msf_west_avg,res)
   res@gsnLeftString="(c)"
    res@gsnRightString="climate_mean_state "
    plot(2)  = gsn_csm_pres_hgt(wks,msf_clm,res) 
    
    
    res@cnFillOn=False
    res@cnLinesOn = True
    res@cnMonoLineColor=True
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLineThicknessF =4.0
    res@cnLineColor="black"
    res@cnLineDashPattern = 0
    res@cnLevels=(/0/)
    res@gsnRightString=""
    res@gsnLeftString=""

    plot1=gsn_csm_pres_hgt(wks,msf_east_avg,res)
    plot2=gsn_csm_pres_hgt(wks,msf_west_avg,res)
    plot3=gsn_csm_pres_hgt(wks,msf_clm,res)

    overlay(plot(0), plot1)
    overlay(plot(1), plot2)
    overlay(plot(2), plot3)



    vcres =res
    vcres@gsnLeftString = ""
    vcres@gsnRightString = ""
    vcres@vcRefMagnitudeF         = 2.0             ; define vector ref magnitude
    vcres@vcRefLengthF            = 0.04             ; define length of vec ref
    vcres@vcMinFracLengthF        = 0.2               ;0.0完全按长度，1.0长度全都一样
    vcres@vcMinDistanceF          = 0.03            ;the distance between the adjacent vector 
    vcres@vcMapDirection = False

    ; vcres@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    vcres@vcGlyphStyle = "CurlyVector"
    vcres@vcLineArrowThicknessF   =1.5   
    vcres@vcFillArrowMinFracWidthF=1.0
    vcres@vcFillArrowHeadMinFracYF=1.0
    vcres@vcFillArrowWidthF=0.05
    vcres@vcFillArrowHeadYF=0.15
   
    ; vcres@vcLevelColors= (/ "gray53","black"/)
    vcres@vcLevelColors= (/ "gray53","black"/)
    vcres@vcMonoLineArrowColor    = False
    vcres@vcLevelSelectionMode= "ExplicitLevels"
    vcres@vcLevels= (/0.05/)

    vcres@vcRefAnnoOrthogonalPosF = -1.0            ;设置参考矢量高度，越负越高
    vcres@vcRefAnnoFontHeightF    = 0.01            ;设置参考矢量大小
    vcres@vcRefAnnoPerimOn      =False

    confi_clm = w_clm
    confi_clm = 0.06
    plot_vc=new(3,graphic)
    plot_vc(0)=gsn_csm_vector_scalar(wks,u_east_avg,w_east_avg,confi_w_east,vcres)
    plot_vc(1)=gsn_csm_vector_scalar(wks,u_west_avg,w_west_avg,confi_w_west,vcres)
    plot_vc(2)=gsn_csm_vector_scalar(wks,u_clm,w_clm,confi_clm, vcres)

    overlay(plot(0), plot_vc(0))
    overlay(plot(1), plot_vc(1))
    overlay(plot(2), plot_vc(2))
;;;

;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 8.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 10)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
    
     ployres@gsLineDashPattern=0
    do i=0,2
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")

   
    plot@$str1$= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    plot@$str2$= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
   end do




;;;panel
     pres=True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar=True
    pres@lbOrientation ="Vertical"
    ; pres@pmLabelBarWidthF=0.5
    gsn_panel(wks, plot, (/1,3/), pres)
  end

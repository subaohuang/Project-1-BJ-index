;;;************************************************************
;; 查看F-pday控制试验以及去掉印度洋海温增温和大西洋增温,太平洋增温的敏感性试验的
;; SLP变化的趋势空间分布。
;; linshuheng 2019.10.22
;;;*********************************************************

;************************************************************
;;       function of getting trend pattern
;************************************************************
function  get_trend_pattern(fils,timestart,timeend)
  begin
    ;;--------------read data-------- ------------
      f1  =  addfile(fils, "r")

      date  =  cd_calendar(f1->time, 1)
      timeind  =  ind(date.le.timeend.and.date.ge.timestart)

      startyear_r = str_split_by_length(timestart, 4)  ;; 将时间string分解成年和月
      endyear_r = str_split_by_length(timeend, 4 )

      startyear = stringtoint(startyear_r(0))
      endyear = stringtoint(endyear_r(0))
      
      PSL_CT  =  f1->PSL(timeind,:,:)
      PSL_CT = (/PSL_CT/100/)
      PSL_CT_a  =  rmMonAnnCycTLL(PSL_CT)

      u_CT  =  f1->U(timeind,{850},:,:)
      
      u_CT_a  =  rmMonAnnCycTLL(u_CT)
  
      v_CT  =  f1->V(timeind,{850},:,:)
      
      v_CT_a  =  rmMonAnnCycTLL(v_CT)

    ;;  

    ;;---------get trend pattern-------------------
      PSL_CT_trend = regCoef_n(ispan(1,dimsizes(PSL_CT_a&time),1),PSL_CT_a(:,:,:),0,0)  
      copy_VarCoords(PSL_CT_a(1,:,:), PSL_CT_trend)

      PSL_CT_trend = (/PSL_CT_trend*10*12/)  ;;; convert to decade trend 

      u_CT_trend = regCoef_n(ispan(1,dimsizes(u_CT_a&time),1),u_CT_a(:,:,:),0,0)  
      copy_VarCoords(u_CT_a(1,:,:), u_CT_trend)
      u_CT_trend = (/u_CT_trend*10*12/)  ;;; convert to decade trend 
      v_CT_trend = regCoef_n(ispan(1,dimsizes(v_CT_a&time),1),v_CT_a(:,:,:),0,0)  
      copy_VarCoords(v_CT_a(1,:,:), v_CT_trend)
      v_CT_trend = (/v_CT_trend*10*12/)  ;;; convert to decade trend 
    ;;

  return [/PSL_CT_trend,u_CT_trend,v_CT_trend/]
 end


;************************************************************
;;       function of t test for trend 
;************************************************************
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


;************************************************************
;;                     main program
;************************************************************


begin 
 ;;------------------read data---------------------------------
	

    fils1 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY/F_1850-PDAY.cam.h0.197801-201212_chazhi.nc"
    
    fils2 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmindian/F_1850-PDAY-rmindian.cam.h0.197801-201212_chazhi.nc"
   
    fils3 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmatlantic/F_1850-PDAY-rmatlan.cam.h0.197801-201212_chazhi.nc"

    ; fils4 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-cam4-rmpacific/F_1850-PDAY-rmpacific.cam.h0.197801-201212_chazhi.nc"
    fils4 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-cam4-rmallocean/F_1850-PDAY-rmallocean.cam.h0.197801-201212_chazhi.nc"

    timestart = 197901
    timeend = 201212
    
    ;;;;get trend pattern
    trend_ct = get_trend_pattern(fils1,timestart,timeend)
    trend_id = get_trend_pattern(fils2,timestart,timeend)
    trend_al = get_trend_pattern(fils3,timestart,timeend)
    trend_pa = get_trend_pattern(fils4,timestart,timeend)

    psl_trend_ct = trend_ct[0]
    psl_trend_id = trend_id[0]
    psl_trend_al = trend_al[0]
    psl_trend_pa = trend_pa[0]
    
    u_trend_ct = trend_ct[1]
    u_trend_id = trend_id[1]
    u_trend_al = trend_al[1]
    u_trend_pa = trend_pa[1]

    v_trend_ct = trend_ct[2]
    v_trend_id = trend_id[2]
    v_trend_al = trend_al[2]
    v_trend_pa = trend_pa[2]
    ;;;;

    ;;;;显著性检验
    pro_trend_ct = rc_test(psl_trend_ct)
    pro_trend_id = rc_test(psl_trend_id)
    pro_trend_al = rc_test(psl_trend_al)
    pro_trend_pa = rc_test(psl_trend_pa)

    pro_u_trend_ct = rc_test(u_trend_ct)
    pro_u_trend_id = rc_test(u_trend_id)
    pro_u_trend_al = rc_test(u_trend_al)
    pro_u_trend_pa = rc_test(u_trend_pa)

    pro_v_trend_ct = rc_test(v_trend_ct)
    pro_v_trend_id = rc_test(v_trend_id)
    pro_v_trend_al = rc_test(v_trend_al)
    pro_v_trend_pa = rc_test(v_trend_pa)
    ;;;
     
    u_trend_ct = where(pro_u_trend_ct.le.0.05, u_trend_ct, u_trend_ct@_FillValue)
    ; v_trend_ct = where(pro_v_trend_ct.le.0.05, v_trend_ct, v_trend_ct@_FillValue)
    
    u_trend_id = where(pro_u_trend_id.le.0.05, u_trend_id, u_trend_id@_FillValue)
    ; v_trend_id = where(pro_v_trend_id.le.0.05, v_trend_id, v_trend_id@_FillValue)
    
    u_trend_al = where(pro_u_trend_al.le.0.05, u_trend_al, u_trend_al@_FillValue)
    ; v_trend_al = where(pro_v_trend_al.le.0.05, v_trend_al, v_trend_al@_FillValue)
    
    u_trend_pa = where(pro_u_trend_pa.le.0.05, u_trend_pa, u_trend_pa@_FillValue)
    ; v_trend_pa = where(pro_v_trend_pa.le.0.05, v_trend_pa, v_trend_pa@_FillValue)
    


  ;;---------------画趋势图--------------------------------------
    ; wks=gsn_open_wks("pdf", "CESM-slp_u_v_pattern-F-pday-rmindian-rmalt-rmpaci")
    wks=gsn_open_wks("pdf", "CESM-slp_u_v_pattern-F-pday-rmindian-rmalt-rmallocean")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4
    res@cnLevelSpacingF      =  0.1

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -60.0
    res@mpMaxLatF             = 60.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.8
    ; res@vpHeightF = 0.2

    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.03           ; thin vectors
    res@vcMinMagnitudeF          =0.04
    res@vcRefMagnitudeF          = 0.6
    res@vcRefLengthF             = 0.04     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"

    ;;;CurlyVector
    res@vcLineArrowThicknessF=1.5
    res@vcMinFracLengthF=0.3
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ;;;
  

    plot=new(4,graphic)
    res@gsnLeftString = "CT"
    plot(0) = gsn_csm_vector_scalar_map(wks,u_trend_ct,v_trend_ct, psl_trend_ct,res)
    res@gsnLeftString = "reindian"
    plot(1) =  gsn_csm_vector_scalar_map(wks,u_trend_id,v_trend_id, psl_trend_id,res)
    res@gsnLeftString = "remove atlantic"
    plot(2) =  gsn_csm_vector_scalar_map(wks,u_trend_al,v_trend_al, psl_trend_al,res)
    res@gsnLeftString = "remove pacific"
    plot(3) =  gsn_csm_vector_scalar_map(wks,u_trend_pa,v_trend_pa, psl_trend_pa,res)

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

    iplot = new(4,graphic)

    iplot(0) = gsn_csm_contour(wks,pro_trend_ct ,res_pro)
    iplot(1) = gsn_csm_contour(wks,pro_trend_id ,res_pro)
    iplot(2) = gsn_csm_contour(wks,pro_trend_al ,res_pro)
    iplot(3) = gsn_csm_contour(wks,pro_trend_pa ,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

    do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
    end do

    ;;;  
    do i = 0,3 
    overlay(plot(i), iplot(i))
    end do 

    pres = True
    pres@gsnPanelLabelBar = True
    
    gsn_panel(wks, plot, (/2,2/), pres)
  end 


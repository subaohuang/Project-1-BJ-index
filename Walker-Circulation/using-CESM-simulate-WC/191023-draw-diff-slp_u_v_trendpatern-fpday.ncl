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
      
      ;;;read data
        vnames =  getfilevarnames(f1)
      
          if any(vnames.eq."U") 
            varid1 = ind(vnames.eq."U")
            varid2 = ind(vnames.eq."V")
            varid3 = ind(vnames.eq."PSL")
          else 
            varid1 = ind(vnames.eq."uwnd")
            varid2 = ind(vnames.eq."vwnd")
            varid3 = ind(vnames.eq."slp")
          end if 
         

          PSL_CT  =  f1->$vnames(varid3)$(timeind,:,:)
        
         if any(vnames.eq."PSL") 
           PSL_CT = (/PSL_CT/100/)
         end if
        PSL_CT_a  =  rmMonAnnCycTLL(PSL_CT)

        u_CT  =  f1->$vnames(varid1)$(timeind,{700},:,:)
        ; u_CT = f1->TAUX(timeind,:,:)

        u_CT_a  =  rmMonAnnCycTLL(u_CT)

        v_CT  =  f1->$vnames(varid2)$(timeind,{700},:,:)
        ; v_CT = F1->TAUY(timeind,:,:)

        v_CT_a  =  rmMonAnnCycTLL(v_CT)

        PSL_CT_a_year = month_to_annual(PSL_CT_a, 1)
        u_CT_a_year = month_to_annual(u_CT_a, 1)
        v_CT_a_year = month_to_annual(v_CT_a, 1)

        timedim = dimsizes(PSL_CT_a_year(:,1,1))
      ;;;

    ;;---------get trend pattern-------------------
      PSL_CT_trend = regCoef_n(ispan(1,timedim,1),PSL_CT_a_year(:,:,:),0,0)  
      copy_VarCoords(PSL_CT_a_year(1,:,:), PSL_CT_trend)

      PSL_CT_trend = (/PSL_CT_trend*10/)  ;;; convert to decade trend 

      u_CT_trend = regCoef_n(ispan(1,timedim,1),u_CT_a_year(:,:,:),0,0)  
      copy_VarCoords(u_CT_a_year(1,:,:), u_CT_trend)
      u_CT_trend = (/u_CT_trend*10/)  ;;; convert to decade trend 
      v_CT_trend = regCoef_n(ispan(1,timedim,1),v_CT_a_year(:,:,:),0,0)  
      copy_VarCoords(v_CT_a_year(1,:,:), v_CT_trend)
      v_CT_trend = (/v_CT_trend*10/)  ;;; convert to decade trend 
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

    fils4 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmpacific/F_1850-PDAY-rmpacific.cam.h0.197801-201212_chazhi.nc"

    fils5 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rafixed/F_1850-PDAY-rafixed.cam.h0.197801-201212_chazhi.nc"

    fils6 = "/home/ys17-19/lsh/data/MRE/mre_u_v_slp_chazhi.nc"

    timestart = 197901
    timeend = 201212
    
    ;;;;get trend pattern
    trend_ct = get_trend_pattern(fils1,timestart,timeend)
    trend_id = get_trend_pattern(fils2,timestart,timeend)
    trend_al = get_trend_pattern(fils3,timestart,timeend)
    trend_pa = get_trend_pattern(fils4,timestart,timeend)
    trend_ghg = get_trend_pattern(fils5,timestart,timeend)
    trend_mre = get_trend_pattern(fils6,timestart,timeend)

    psl_trend_ct = trend_ct[0]
    psl_trend_id = trend_id[0]
    psl_trend_al = trend_al[0]
    psl_trend_pa = trend_pa[0]
    psl_trend_ghg = trend_ghg[0]
    psl_trend_mre = trend_mre[0]

    
    u_trend_ct = trend_ct[1]
    u_trend_id = trend_id[1]
    u_trend_al = trend_al[1]
    u_trend_pa = trend_pa[1]
    u_trend_ghg = trend_ghg[1]
    u_trend_mre = trend_mre[1]

    v_trend_ct = trend_ct[2]
    v_trend_id = trend_id[2]
    v_trend_al = trend_al[2]
    v_trend_pa = trend_pa[2]
    v_trend_ghg = trend_ghg[2]
    v_trend_mre = trend_mre[2]
    ;;;;

    ;;;diff with Control run
    psl_trend_id = (/psl_trend_ct - psl_trend_id/)
    psl_trend_pa = (/psl_trend_ct - psl_trend_pa/)
    psl_trend_al = (/psl_trend_ct - psl_trend_al/)
    psl_trend_ghg = (/psl_trend_ct - psl_trend_ghg/)


    u_trend_id = (/u_trend_ct - u_trend_id/)
    u_trend_pa = (/u_trend_ct - u_trend_pa/)
    u_trend_al = (/u_trend_ct - u_trend_al/)
    u_trend_ghg = (/u_trend_ct - u_trend_ghg/)

    v_trend_id = (/v_trend_ct - v_trend_id/)
    v_trend_pa = (/v_trend_ct - v_trend_pa/)
    v_trend_al = (/v_trend_ct - v_trend_al/)
    v_trend_ghg = (/v_trend_ct - v_trend_ghg/)
    ;;;
   
    
    u_trend_ct_wsp = avg(dim_avg_n(u_trend_ct({-5:5},{150:210}), (/1/)))
    u_trend_id_wsp = avg(dim_avg_n(u_trend_id({-5:5},{150:210}), (/1/)))
    u_trend_al_wsp = avg(dim_avg_n(u_trend_al({-5:5},{150:210}), (/1/)))
    u_trend_pa_wsp = avg(dim_avg_n(u_trend_pa({-5:5},{150:210}), (/1/)))
    u_trend_ghg_wsp = avg(dim_avg_n(u_trend_ghg({-5:5},{150:210}), (/1/)))

    u_trend_ct_esp = avg(dim_avg_n(u_trend_ct({5:5},{240:270}), (/1/)))
    u_trend_id_esp = avg(dim_avg_n(u_trend_id({5:5},{240:270}), (/1/)))
    u_trend_al_esp = avg(dim_avg_n(u_trend_al({5:5},{240:270}), (/1/)))
    u_trend_pa_esp = avg(dim_avg_n(u_trend_pa({5:5},{240:270}), (/1/)))


    print("u_trend_ct_wsp="+u_trend_ct_wsp)
    print("u_trend_id_wsp="+u_trend_id_wsp)
    print("u_trend_al_wsp="+u_trend_al_wsp)
    print("u_trend_pa_wsp="+u_trend_pa_wsp)
    print("u_trend_ghg_wsp="+u_trend_ghg_wsp)

    print("u_trend_ct_esp="+u_trend_ct_esp)
    print("u_trend_id_esp="+u_trend_id_esp)
    print("u_trend_al_esp="+u_trend_al_esp)
    print("u_trend_pa_esp="+u_trend_pa_esp)
exit()

  ;;---------------画趋势图--------------------------------------
    wks=gsn_open_wks("pdf", "DIFF-slp_u_v_pattern-F-pday-rmssttrend-rafixed")
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

    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.6
    res@vpHeightF = 0.3

    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.03           ; thin vectors
    res@vcMinMagnitudeF          =0.05
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
  

    plot=new(6,graphic)
    res@gsnLeftStringFontHeightF = 20.0
    res@gsnLeftString = "Observation"
    plot(0) = gsn_csm_vector_scalar_map(wks,u_trend_mre,v_trend_mre , psl_trend_mre,res)
    res@gsnLeftString = "Control run"
    plot(1) = gsn_csm_vector_scalar_map(wks,u_trend_ct,v_trend_ct, psl_trend_ct,res)
    ; res@gsnLeftString = "Control - Indian warming removed" 
    res@gsnLeftString = "Indian warming "
    plot(2) =  gsn_csm_vector_scalar_map(wks,u_trend_id,v_trend_id, psl_trend_id,res)
    ; res@gsnLeftString = "Control - altlantic warming removed"
    res@gsnLeftString = "Atlantic warming"
    plot(3) =  gsn_csm_vector_scalar_map(wks,u_trend_al,v_trend_al, psl_trend_al,res)
    ; res@gsnLeftString = "Control - pacific SST trend removed"
    res@gsnLeftString = "Pacific Lanila-like SST trend "
    plot(4) =  gsn_csm_vector_scalar_map(wks,u_trend_pa,v_trend_pa, psl_trend_pa,res)
    ; res@gsnLeftString = "Control - GHG fixed"
    res@gsnLeftString = "GHG effect"
    plot(5) =  gsn_csm_vector_scalar_map(wks,u_trend_ghg,v_trend_ghg, psl_trend_ghg,res)

    pres = True
    pres@gsnPanelLabelBar = True
    pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
    
    gsn_panel(wks, plot, (/3,2/), pres)



  end 


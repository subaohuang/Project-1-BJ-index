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

      print("****************"+fils+"*******************8")
      
      if fils.eq."/home/ys17-19/lsh/data/rea_msf/msf_time_mre_global.nc"
        msf_r = f1->msf_time(:,:,timeind)
        msf_CT = msf_r(time|:,level|:,lon|:)
        
        w_r = f1->w_time(:,:,timeind)
        w_CT = w_r(time|:,level|:,lon|:)
        
        ud_r = f1->ud_time(:,:,timeind)
        ud_CT = ud_r(time|:,level|:,lon|:)
        delete(msf_r)
        delete(w_r)
        delete(ud_r)

      else 
        msf_CT  =  f1->msf_time(timeind,{1000:100},:)
        w_CT  =  f1->w_time(timeind,{1000:100},:)
        ud_CT  =  f1->ud_time(timeind,{1000:100},:)
      end if

      ;;;get anomaly       
      msf_CT = (/msf_CT/10^9/)
      msf_CT_a  =  rmMonAnnCycTLL(msf_CT)
      w_CT_a = rmMonAnnCycTLL(w_CT)
      ud_CT_a = rmMonAnnCycTLL(ud_CT)


      ;;;取年平均
      msf_CT_a_year = month_to_annual(msf_CT_a, 1)
      ud_CT_a_year = month_to_annual(ud_CT_a, 1)
      w_CT_a_year = month_to_annual(w_CT_a, 1)
      ; msf_CT_a_year = month_to_season(msf_CT_a, "SON")  
      ; ud_CT_a_year = month_to_season(ud_CT_a, "SON")
      ; w_CT_a_year = month_to_season(w_CT_a, "SON")

      timedim = dimsizes(msf_CT_a_year(:,1,1))

      
    
    ;;  
    ;;---------get trend pattern-------------------
      msf_CT_trend = regCoef_n(ispan(1,timedim,1),msf_CT_a_year(:,:,:),0,0)  
      copy_VarCoords(msf_CT_a_year(1,:,:), msf_CT_trend)
      printVarSummary(msf_CT_trend)

      msf_CT_trend = (/msf_CT_trend*10/)  ;;; convert to decade trend 

      w_CT_trend = regCoef_n(ispan(1,timedim,1),w_CT_a_year(:,:,:),0,0)  
      copy_VarCoords(w_CT_a_year(1,:,:), w_CT_trend)
      printVarSummary(w_CT_trend)
      w_CT_trend = (/w_CT_trend*10/)  ;;; convert to decade trend 

      ud_CT_trend = regCoef_n(ispan(1,timedim,1),ud_CT_a_year(:,:,:),0,0)  
      copy_VarCoords(ud_CT_a_year(1,:,:), ud_CT_trend)
      printVarSummary(ud_CT_trend)
      ud_CT_trend = (/ud_CT_trend*10/)  ;;; convert to decade trend 

      msf_CT_clm = dim_avg_n_Wrap(msf_CT,0)
    ;;
    ;;
  return  [/msf_CT_trend,w_CT_trend,ud_CT_trend,msf_CT_clm/]
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

function  get_topo(plot_diff,wks)
begin

  ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 30, 25)
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon6=fspan(30, 39, 25)


      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel6=conform_dims(dimsizes(ploylon6), ploylevel, 0)
    

      ployres@gsLineDashPattern=0
    
    
      str1 = unique_string("dum$m(i)$")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
      str2 = unique_string("dum2$m(i)$")
      str3 = unique_string("dum3$m(i)$")

      plot_diff@$str1$= gsn_add_polyline(wks,plot_diff,  ploylon1,ploylevel1,  ployres)

      ployres@gsLineDashPattern=16

      plot_diff@$str2$= gsn_add_polyline(wks,plot_diff,  ploylon2,ploylevel2,  ployres)

      ployres@gsLineDashPattern=0
      plot_diff@$str3$= gsn_add_polyline(wks,plot_diff,  ploylon3,ploylevel3,  ployres)
return plot_diff
end 


;************************************************************
;;                     main program
;************************************************************


begin 
 ;;------------------read data---------------------------------
	

  fils1  =  "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY_1979-201212.nc"

  fils2 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmindian_197901-201212.nc"

  fils3 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmatlan_197912-201212.nc"

  fils4 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmpacific_197912-201212.nc"

  fils5 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rafixed_197912-201212.nc"
  
  fils6 = "/home/ys17-19/lsh/data/rea_msf/msf_time_mre_global.nc"
  

  

    timestart = 197901  
    timeend = 201212
    
    ;;;;get trend pattern
    trend_ct = get_trend_pattern(fils1,timestart,timeend)
    trend_id = get_trend_pattern(fils2,timestart,timeend)
    trend_al = get_trend_pattern(fils3,timestart,timeend)
    trend_pa = get_trend_pattern(fils4,timestart,timeend)
    trend_ghg = get_trend_pattern(fils5,timestart,timeend)
    trend_mre = get_trend_pattern(fils6,timestart,timeend)
    ;;;;
    print("**************"+"done"+"*******************")

    ;;;;get trend pattern
    msf_trend_ct = trend_ct[0]
    msf_trend_id = trend_id[0]
    msf_trend_al = trend_al[0]
    msf_trend_pa = trend_pa[0]
    msf_trend_ghg = trend_ghg[0]
    msf_trend_mre = trend_mre[0]
    ;;;;

    ;;;;get trend pattern
    w_trend_ct = trend_ct[1]
    w_trend_id = trend_id[1]
    w_trend_al = trend_al[1]
    w_trend_pa = trend_pa[1]
    w_trend_ghg = trend_ghg[1]
    w_trend_mre = trend_mre[1]
    
    ud_trend_ct = trend_ct[2]
    ud_trend_id = trend_id[2]
    ud_trend_al = trend_al[2]
    ud_trend_pa = trend_pa[2]
    ud_trend_ghg = trend_ghg[2]
    ud_trend_mre = trend_mre[2]

    ;;;; 

    ;;;;get clm pattern
    msf_clm_ct = trend_ct[3]
    msf_clm_id = trend_id[3]
    msf_clm_al = trend_al[3]
    msf_clm_pa = trend_pa[3]
    msf_clm_ghg = trend_ghg[3]
    msf_clm_mre = trend_mre[3]
    ;;;;

    ; ;;;;显著性检验
    ; pro_trend_ct = rc_test(msf_trend_ct)
    ; pro_trend_id = rc_test(msf_trend_id)
    ; pro_trend_al = rc_test(msf_trend_al)
    ; pro_trend_pa = rc_test(msf_trend_pa)
    ; ;;;

    ;;;diff from control run
    diff_id = (/msf_trend_ct - msf_trend_id/)
    diff_al = (/msf_trend_ct - msf_trend_al/)
    diff_pa = (/msf_trend_ct - msf_trend_pa/)
    diff_ghg = (/msf_trend_ct - msf_trend_ghg/)
    ;;;

    copy_VarCoords(msf_trend_id, diff_id)
    copy_VarCoords(msf_trend_al, diff_al)
    copy_VarCoords(msf_trend_pa, diff_pa)
    copy_VarCoords(msf_trend_ghg, diff_ghg)
    ;;;


    ;;;diff w from control run
    w_diff_id = (/w_trend_ct - w_trend_id/)
    w_diff_al = (/w_trend_ct - w_trend_al/)
    w_diff_pa = (/w_trend_ct - w_trend_pa/)
    w_diff_ghg = (/w_trend_ct - w_trend_ghg/)
    ;;;
  
    copy_VarCoords(w_trend_id, w_diff_id)
    copy_VarCoords(w_trend_al, w_diff_al)
    copy_VarCoords(w_trend_pa, w_diff_pa)
    copy_VarCoords(w_trend_ghg, w_diff_ghg)
    ;;;
    

     ;;;diff ud from control run
    ud_diff_id = (/ud_trend_ct - ud_trend_id/)
    ud_diff_al = (/ud_trend_ct - ud_trend_al/)
    ud_diff_pa = (/ud_trend_ct - ud_trend_pa/)
    ud_diff_ghg = (/ud_trend_ct - ud_trend_ghg/)
    ;;;
  
    copy_VarCoords(ud_trend_id, ud_diff_id)
    copy_VarCoords(ud_trend_al, ud_diff_al)
    copy_VarCoords(ud_trend_pa, ud_diff_pa)
    copy_VarCoords(ud_trend_ghg, ud_diff_ghg)
    ;;;

    ud_diff_ct_wsp = avg(ud_trend_ct({700},{150:210}))
    ud_diff_id_wsp = avg(ud_diff_id({700},{150:210}))
    ud_diff_al_wsp = avg(ud_diff_al({700},{150:210}))
    ud_diff_pa_wsp = avg(ud_diff_pa({700},{150:210}))
    ud_diff_ghg_wsp = avg(ud_diff_ghg({700},{150:210}))

    print("ud_diff_ct_wsp="+ud_diff_ct_wsp)
    print("ud_diff_id_wsp="+ud_diff_id_wsp)
    print("ud_diff_al_wsp="+ud_diff_al_wsp)
    print("ud_diff_pa_wsp="+ud_diff_pa_wsp)
    print("ud_diff_ghg_wsp="+ud_diff_ghg_wsp)
; exit()
    
  ;;---------------画趋势图--------------------------------------
      print("**************"+"画图"+"*******************")

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnFillPalette = "BlRe"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.5
    res@cnMaxLevelValF = 0.5
    res@cnLevelSpacingF = 0.1
    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False


    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
  
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)

    ; res@tmXBValues=ispan(100,310,40)
    ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","60~S~o~N~W"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(100,310,10)
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7


    ; res@vcRefMagnitudeF = 0.2               ; define vector ref mag
    ; res@vcRefLengthF    = 0.04       ; define length of vec ref
    ; res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    ; res@vcMinDistanceF  = 0.03      ; thin out vectors
    ; res@vcMapDirection  = False
    ; res@vcLineArrowThicknessF=2.0
    ; ; res@vcFillArrowHeadXF=0.1
    ; res@vcFillArrowMinFracWidthF=1.0
    ; res@vcFillArrowHeadMinFracYF=1.0
    ; ; res@vcFillArrowHeadXF=1.0
    ; res@vcFillArrowWidthF=0.05
    ; res@vcFillArrowHeadYF=0.15
    
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
    res@vcRefMagnitudeF = 0.1              ; define vector ref mag
    res@vcRefLengthF    = 0.05        ; define length of vec ref

    res@vcGlyphStyle    = "CurlyVector"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.03       ; thin out vectors
    res@vcMapDirection  = False
    res@vcLineArrowThicknessF=2.0
    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15
    
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005

    res@vcVectorDrawOrder="PostDraw"

    res@vcRefAnnoOrthogonalPosF  = -1.14

    res@vcRefAnnoPerimOn=False
    res@vcRefAnnoString2On=False

    ;;;;


   
    ; res@cnFillPalette = "MPL_seismic"
     wks2 = gsn_open_wks("pdf", "年平均-CESM-diff_msf_pattern-F-pday-rmindian-rmalt-rmpacific")
     ; wks2 = gsn_open_wks("pdf", "秋季-CESM-diff_msf_pattern-F-pday-rmindian-rmalt-rmpacific")
    
    plot_diff = new(6,graphic)
    res@gsnLeftString = "Observation"
    ; print(w_trend_ct)

  
    plot_diff(0) = gsn_csm_pres_hgt_vector(wks2,msf_trend_mre,ud_trend_mre,w_trend_ct,res)
    ; plot_diff(0) = gsn_csm_pres_hgt(wks2,msf_trend_mre,res)
    
    print("**************"+"画图1完成"+"*******************")
    res@gsnLeftString = "Control run"
    plot_diff(1) = gsn_csm_pres_hgt_vector(wks2,msf_trend_ct,ud_trend_ct,w_trend_ct,res)
    ; res@gsnLeftString = "Control - Indian warming removed"
    res@gsnLeftString = "Indian warming"
    plot_diff(2) = gsn_csm_pres_hgt_vector(wks2,diff_id,ud_diff_id ,w_diff_id ,res)
    ; res@gsnLeftString = "Control - atlantic warming removed"
    res@gsnLeftString = "Atlantic warming "
    plot_diff(3) = gsn_csm_pres_hgt_vector(wks2, diff_al,ud_diff_al,w_diff_al,res)
     ; res@gsnLeftString = "Control - pacific SST trend removed"
     res@gsnLeftString = "Pacific Lanina-like SST trend "
    plot_diff(4) = gsn_csm_pres_hgt_vector(wks2,diff_pa,ud_diff_pa,w_diff_pa,res)
     res@gsnLeftString = "GHG effect "
    plot_diff(5) = gsn_csm_pres_hgt_vector(wks2,diff_ghg,ud_diff_ghg,w_diff_ghg,res)

    res@cnMinLevelValF = -6.0
    res@cnMaxLevelValF = 6.0
    res@cnLevelSpacingF = 1.0
    res@cnFillOn =False
    res@cnLinesOn = True
    res@gsnLeftString = ""
    res@gsnContourNegLineDashPattern = 1
    res@cnLineThicknesses = 3.0

     
    plot_clm = new(6,graphic)
    plot_clm(0) = gsn_csm_pres_hgt(wks2, msf_clm_mre, res)
    plot_clm(1) = gsn_csm_pres_hgt(wks2, msf_clm_ct, res)
    plot_clm(2) = gsn_csm_pres_hgt(wks2, msf_clm_id, res)
    plot_clm(3) = gsn_csm_pres_hgt(wks2, msf_clm_al, res)
    plot_clm(4) = gsn_csm_pres_hgt(wks2, msf_clm_pa, res)
    plot_clm(5) = gsn_csm_pres_hgt(wks2, msf_clm_ghg, res)
   
    do i = 0, 5
     overlay(plot_diff(i), plot_clm(i))
    end do 

    do i = 0,5
      plot_diff(i) = get_topo(plot_diff(i),wks2)
    end do

    pres =True
    pres@gsnPanelLabelBar  = True
    pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
    gsn_panel(wks2, plot_diff, (/3,2/), pres)

  end 


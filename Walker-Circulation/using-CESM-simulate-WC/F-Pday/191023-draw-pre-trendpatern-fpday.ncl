;;;************************************************************
;; 查看F-pday控制试验以及去掉印度洋海温增温和大西洋增温,太平洋增温的敏感性试验的
;; 降水的趋势空间分布。
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
      if fils.ne."/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc"

      PRECL_CT  =  f1->PRECL(timeind,:,:)
      PRECC_CT  =  f1->PRECC(timeind,:,:)
     
      PRE_CT = PRECC_CT
      PRE_CT = (/PRECC_CT+PRECL_CT/)
      PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day
      else 
      PRE_CT = f1->precip(timeind,:,:)
      end if
      PRE_CT_a  =  rmMonAnnCycTLL(PRE_CT)
      

    ;;---------get trend pattern-------------------
      PRE_CT_trend = regCoef_n(ispan(1,dimsizes(PRE_CT_a&time),1),PRE_CT_a(:,:,:),0,0)  
      copy_VarCoords(PRE_CT_a(1,:,:), PRE_CT_trend)

      PRE_CT_trend = (/PRE_CT_trend*10*12/)  ;;; convert to decade trend 
    ;;
  return PRE_CT_trend
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
	

    fils1 = "/home/ys17-19/lsh/CESM-data/F/F_1850-PDAY/F_1850-PDAY.cam.h0.197801-201212.pre.nc"
    
    fils2 = "/home/ys17-19/lsh/CESM-data/F/F_1850-PDAY-rmindian/F_1850-PDAY-rmindian.cam.h0.197801-201212.pre.nc"
   
    fils3 = "/home/ys17-19/lsh/CESM-data/F/F_1850-PDAY-rmatlantic/F_1850-PDAY-rmatlan.cam.h0.197801-201212.pre.nc"

    fils4 = "/home/ys17-19/lsh/CESM-data/F/F_1850-PDAY-rmpacific/F_1850-PDAY-rmpacific.cam.h0.197801-201212.nc"

    fils5 = "/home/ys17-19/lsh/CESM-data/F/F_1850-PDAY-rafixed/F_1850-PDAY-rafixed.cam.h0.197801-201212.nc"

    fils6 = "/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc"

    timestart = 197901
    timeend = 201212
    
    ;;;;get trend pattern
    PRE_trend_ct = get_trend_pattern(fils1,timestart,timeend)
    PRE_trend_id = get_trend_pattern(fils2,timestart,timeend)
    PRE_trend_al = get_trend_pattern(fils3,timestart,timeend)
    PRE_trend_pa = get_trend_pattern(fils4,timestart,timeend)
    PRE_trend_ghg = get_trend_pattern(fils5,timestart,timeend)
    PRE_trend_mre = get_trend_pattern(fils6,timestart,timeend)
    ;;;;
      
     PRE_trend_id = (/PRE_trend_ct - PRE_trend_id/)
     PRE_trend_al = (/PRE_trend_ct - PRE_trend_al/)
     PRE_trend_pa = (/PRE_trend_ct - PRE_trend_pa/) 
     PRE_trend_ghg = (/PRE_trend_ct - PRE_trend_ghg/) 


    ;;;;显著性检验
    pro_trend_ct = rc_test(PRE_trend_ct)
    pro_trend_id = rc_test(PRE_trend_id)
    pro_trend_al = rc_test(PRE_trend_al)
    pro_trend_pa = rc_test(PRE_trend_pa)
    ;;;;;



  ;;---------------画趋势图--------------------------------------
    dir = systemfunc("/bin/pwd")
   
    wks=gsn_open_wks("pdf", dir+"/fig/"+"CESM-prec_diff_pattern-F-pday-ssttrend-DIFF")
    ; wks=gsn_open_wks("pdf", dir+"/fig/"+"CESM-prec_diff_pattern-F-pday-rmindian-rmalt-rmallocean_DIFF")
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

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "MPL_BrBG"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.8
    ; res@vpHeightF = 0.2

    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.6
    ; res@vpHeightF = 0.3

    plot=new(6,graphic)
    res@gsnLeftStringFontHeightF = 20.0
    res@gsnLeftString = "Observation"
    plot(0) = gsn_csm_contour_map(wks,PRE_trend_mre,res)
    res@gsnLeftString = "Control run"
    plot(1) = gsn_csm_contour_map(wks,PRE_trend_ct,res)
    ; res@gsnLeftString = "Control run - indian warming removed"
    res@gsnLeftString = "Indian warming"
    plot(2) = gsn_csm_contour_map(wks,PRE_trend_id,res)
    ; res@gsnLeftString = "Control run - atlantic warming removed"
    res@gsnLeftString = "Atlantic warming"
    plot(3) = gsn_csm_contour_map(wks,PRE_trend_al,res)
    ; res@gsnLeftString = "Control run - pacific sst trend removed"
    res@gsnLeftString = "Pacific Lanina-like SST trend"
    plot(4) = gsn_csm_contour_map(wks,PRE_trend_pa,res)
    res@gsnLeftString = "GHG effect"
    plot(5) = gsn_csm_contour_map(wks,PRE_trend_ghg,res)

    ; ;;;打点
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False

    ; iplot = new(6,graphic)

    ; iplot(0) = gsn_csm_contour(wks,pro_trend_ct ,res_pro)
    ; iplot(1) = gsn_csm_contour(wks,pro_trend_id ,res_pro)
    ; iplot(2) = gsn_csm_contour(wks,pro_trend_al ,res_pro)
    ; iplot(3) = gsn_csm_contour(wks,pro_trend_pa ,res_pro)


    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; opt@gsnShadeMid     = 17         ; stipple
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    ; opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

    ; do i=0,3
    ; iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
    ; end do

    ;;;  
    ; do i = 0,3 
    ; overlay(plot(i), iplot(i))
    ; end do 

    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
    
    gsn_panel(wks, plot, (/3,2/), pres)
  end 


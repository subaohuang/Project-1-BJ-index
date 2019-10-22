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
      
      msf_CT  =  f1->msf_time(timeind,{1000:100},:)
      msf_CT = (/msf_CT/10^9/)

      msf_CT_a  =  rmMonAnnCycTLL(msf_CT)
    ;;  
    ;;---------get trend pattern-------------------
      msf_CT_trend = regCoef_n(ispan(1,dimsizes(msf_CT_a&time),1),msf_CT_a(:,:,:),0,0)  
      copy_VarCoords(msf_CT_a(1,:,:), msf_CT_trend)
      printVarSummary(msf_CT_trend)

      msf_CT_trend = (/msf_CT_trend*10*12/)  ;;; convert to decade trend 
    ;;

  return msf_CT_trend
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
	

  fils1  =  "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY_1979-201212.nc"

  fils2 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmindian_197901-201212.nc"

  fils3 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmatlan_197912-201212.nc"

  fils4 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmpacific_197912-201212.nc"

    timestart = 197901
    timeend = 201212
    
    ;;;;get trend pattern
    msf_trend_ct = get_trend_pattern(fils1,timestart,timeend)
    msf_trend_id = get_trend_pattern(fils2,timestart,timeend)
    msf_trend_al = get_trend_pattern(fils3,timestart,timeend)
    msf_trend_pa = get_trend_pattern(fils4,timestart,timeend)
    ;;;;

    ;;;;显著性检验
    pro_trend_ct = rc_test(msf_trend_ct)
    pro_trend_id = rc_test(msf_trend_id)
    pro_trend_al = rc_test(msf_trend_al)
    pro_trend_pa = rc_test(msf_trend_pa)
    ;;;

    ;;;diff from control run
    diff_id = (/msf_trend_id - msf_trend_ct/)
    diff_al = (/msf_trend_al - msf_trend_ct/)
    diff_pa = (/msf_trend_pa - msf_trend_ct/)

    copy_VarCoords(msf_trend_id, diff_id)
    copy_VarCoords(msf_trend_al, diff_al)
    copy_VarCoords(msf_trend_pa, diff_pa)
    ;;;

    
  ;;---------------画趋势图--------------------------------------
    wks=gsn_open_wks("pdf", "CESM-msf_pattern-F-pday-rmindian-rmalt-rmpaci")
    

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLinesOn=True
 
    res@cnFillOn =True
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.5
    res@cnMaxLevelValF = 0.5
    res@cnLevelSpacingF = 0.1
    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False


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
   

    plot=new(4,graphic)
    res@gsnLeftString = "CT"
    plot(0) = gsn_csm_pres_hgt(wks,msf_trend_ct,res)
    res@gsnLeftString = "remove indian"
    plot(1) = gsn_csm_pres_hgt(wks,msf_trend_id,res)
    res@gsnLeftString = "remove atlantic"
    plot(2) = gsn_csm_pres_hgt(wks,msf_trend_al,res)
    res@gsnLeftString = "remove pacific"
    plot(3) = gsn_csm_pres_hgt(wks,msf_trend_pa,res)

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

    iplot(0) = gsn_csm_pres_hgt(wks,pro_trend_ct ,res_pro)
    iplot(1) = gsn_csm_pres_hgt(wks,pro_trend_id ,res_pro)
    iplot(2) = gsn_csm_pres_hgt(wks,pro_trend_al ,res_pro)
    iplot(3) = gsn_csm_pres_hgt(wks,pro_trend_pa ,res_pro)


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

  ;;;; 
    res@cnFillPalette = "BlueRed"
    wks2 = gsn_open_wks("pdf", "CESM-msf_pattern-F-pday-diff")
    plot_diff = new(4,graphic)
    res@gsnLeftString = "CT"
    plot_diff(0) = gsn_csm_pres_hgt(wks2,msf_trend_ct,res)
    res@gsnLeftString = "remove indian"
    plot_diff(1) = gsn_csm_pres_hgt(wks2,diff_id,res)
    res@gsnLeftString = "remove atlantic"
    plot_diff(2) = gsn_csm_pres_hgt(wks2,diff_al,res)
    res@gsnLeftString = "remove pacific"
    plot_diff(3) = gsn_csm_pres_hgt(wks2,diff_pa,res)

 
    
        
    gsn_panel(wks2, plot_diff, (/2,2/), pres)



  end 


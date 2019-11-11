;;;**********************************************************************
;;;计算三个环流的msf强度与海温的回归场
;;;2019.10.29
;;;by shuheng lin
;;;**********************************************************
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
;;;;
begin
;;;读入HadiSST海温数据
  diri="/home/ys17-19/lsh/data/sst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)

  timeind=ind(date.le.201212.and.date.ge.197901)
  
  sst=lonFlip(f1->sst(timeind,{60:-60},:))
  ssta=(/rmMonAnnCycTLL(sst)/)
  copy_VarCoords(sst, ssta)
  ssta=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
  ; ssta_son=month_to_season(ssta_de, "SON")
  ; ssta_djf=month_to_season(ssta_de, "DJF")
  ssta_annual = month_to_annual(ssta, 1)
  ; ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
  ; copy_VarCoords(ssta_son(0:54,:,:), ssta_s)

  
  ssta_annual_id = dim_avg_n_Wrap(ssta_annual(:,{-20:20},{30:100}), (/1,2/))


  ssta_annual_pa =  dim_avg_n_Wrap(ssta_annual(:,{5:-5},{210:270}), (/1,2/))\
    -  (dim_avg_n_Wrap(ssta_annual(:,{45:25},{150:210}), (/1,2/))+dim_avg_n_Wrap(ssta_annual(:,{-25:-45},{150:210}), (/1,2/)))/2

  ssta_annual_ln = lonFlip(ssta_annual)
  ssta_annual_al = dim_avg_n_Wrap(ssta_annual_ln(:,{-30:60},{-80:20}), (/1,2/))



  delete(date)
  delete(timeind)


;;;读入msf的数据

  fils1  =  "/home/ys17-19/lsh/data/rea_msf/msf_time_era_global.nc"
  f1  =  addfile(fils1, "r") 
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201212.and.date.ge.197901)

  msf_time = f1->msf_time(:,:,timeind)

  msf_time = (/msf_time/10^9/)

  msf_time = (/dtrend_msg_n(ispan(1,dimsizes(msf_time&time),1),msf_time,False,True,2)/)

  msf_time =(/rmMonAnnCycLLT(msf_time)/)

  msf_time_annual = month_to_annual(msf_time(time|:,level|:,lon|:), 1)


  delete(date)
  delete(timeind)
;;;

;;;滑动平均，提取年代际变化信号
  ssta_annual_id= (/runave_n(ssta_annual_id, 9, 0, 0)/)
  ssta_annual_pa = (/runave_n(ssta_annual_pa,9, 0, 0)/)
  ssta_annual_al = (/runave_n(ssta_annual_al, 9, 0, 0)/)

  msf_time_annual= (/runave_n(msf_time_annual,9, 0, 0)/)

;;;



; ;;;回归sst场

  msf_id_rc=regCoef_n(ssta_annual_id,msf_time_annual,0,0)  
  copy_VarCoords(msf_time_annual(0,:,:), msf_id_rc)

  msf_pa_rc=regCoef_n(ssta_annual_pa,msf_time_annual,0,0)  
  copy_VarCoords(msf_time_annual(0,:,:), msf_pa_rc)

  msf_al_rc=regCoef_n(ssta_annual_al,msf_time_annual,0,0)  
  copy_VarCoords(msf_time_annual(0,:,:), msf_al_rc)
;;;

  
; ;;;回归u场,v场

;   ssta_rc=regCoef_n(ipindex2,ssta_s,0,0)  
;   copy_VarCoords(ssta_s(0,:,:), ssta_rc)

;   ua_rc=regCoef_n(ipindex2,ua_cold,0,0)  
;   copy_VarCoords(ua(0,:,:), ua_rc)

;   va_rc=regCoef_n(ipindex2,va_cold,0,0)  
;   copy_VarCoords(va(0,:,:), va_rc)
; ;;;显著性检验
;   prob_ssta=rc_test(ssta_rc)
;   prob_ua=rc_test(ua_rc)
;   prob_va=rc_test(va_rc)
; ;;;
; ;;;

; ;;;显著性检验
;   prob_ssta=rtest(ssta_rc,dimsizes(ssta_s&time),0)
;   prob_ua=rtest(ua_rc,dimsizes(ua_cold&time),0)
;   prob_va=rtest(va_rc,dimsizes(va_cold&time),0)
;   prob_ssta=where(prob_ssta.eq.0, prob_ssta@_FillValue, prob_ssta)
;   copy_VarCoords(ssta_rc, prob_ssta)
;   copy_VarCoords(ua_rc, prob_ua)
;   copy_VarCoords(va_rc, prob_va)
; ;;;



; ;;;风场的检验
;   ua_rc=where(prob_ua.le.0.05, ua_rc, ua_rc@_FillValue)
;   va_rc=where(prob_va.le.0.05, va_rc, va_rc@_FillValue)
;   printVarSummary(ua_rc)
;   printVarSummary(va_rc)
; ;;;



;;;画图
wks = gsn_open_wks("pdf", "各海盆sst强度指数回归msf")
  gsn_define_colormap(wks, "BlueWhiteOrangeRed")
  
 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    ; res@cnLevelSelectionMode ="ManualLevels"
    ; res@cnMinLevelValF = -0.5
    ; res@cnMaxLevelValF = 0.5
    ; res@cnLevelSpacingF = 0.1
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
   

    plot=new(3,graphic)
    res@gsnLeftString = "regression to inidan sst "
    plot(0) = gsn_csm_pres_hgt(wks,msf_id_rc,res)
    res@gsnLeftString = "regression to pacific sst"
    plot(1) = gsn_csm_pres_hgt(wks,msf_pa_rc,res)
     res@gsnLeftString = "regression to atlantic sst"
    plot(2) = gsn_csm_pres_hgt(wks,msf_al_rc,res)
   
pres = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/2,2/), pres)
end



;;********************************************************
;;计算ERA的u,v,slp的场到P,  Q-term以及其动力项，还有观测的U,V,SLP的趋势
;;
;;lnshuheng 2020.1.4
;**********************************************************
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

begin
;;;读入GPCP降水场
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.201312)
    p_t=f2->precip(timeind,:,:)
;;;;

;==========================================================
;;;回归SST
;;==========================================================

   ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
   
   
   ;;;
    ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
    ;;;

    ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;


    f1 = addfile(fils1,"r")
    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    sst=lonFlip(f1->sst(timeind,:,:))
    ;;;


    ;;;
    sst = (/rmMonAnnCycTLL(sst)/)
    p_t = (/rmMonAnnCycTLL(p_t)/)



     ;;;choose season

     season = (/"MAM","JJA","SON","DJF"/)


     plot=new(8,graphic)
     iplot =new(8,graphic)
     str = ispan(1, 8, 1) +""
     wks   = gsn_open_wks("pdf", "sst_congruent_on_P-"+"allseason")

 do i=0,3
   ;;;处理成年异常数据
    
     sst_annual = month_to_season(sst, season(i))
     
     p = dim_avg_n_Wrap(p_t(:,{0:15},{110:155}), (/1,2/))
     p_annual = month_to_season(p, season(i))
     p_annual = (/dim_standardize(p_annual, 1)/)
   

     ;;;计算趋势,计算线性一致性之前先把趋势计算了
     
     sst_rc  = regCoef_n(ispan(1, yearnum, 1), sst_annual, 0, 0)
     ; sst_rc = sst_rc*10

     sst_rc = sst_rc*yearnum
     sst_rc = where(sst_rc.eq.0, sst_rc@_FillValue, sst_rc)

     copy_VarCoords(sst_annual(1,:,:), sst_rc)

    
     p_rc  = regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
     p_rc = p_rc*yearnum
  

     sst_annual = (/dtrend_msg_n(ispan(1,dim(0),1),sst_annual,False,True,0)/)
     p_annual = (/dtrend_msg_n(ispan(1,dim(0),1),p_annual,False,True,0)/)


  ;;回归
    ;;1: 回归到P
    sst_p_rc = regCoef_n(p_annual, sst_annual, 0, 0)
    copy_VarCoords(sst_annual(1,:,:), sst_p_rc)
    
   ;;;


   ;;;显著性检验
    
    prob_sst = rc_test(sst_rc)
  
    prob_sst_p = rc_test(sst_p_rc)

    prob_sst_p = prob_sst
   ;;;;

  ;;;对线性一致性进行f检验
    


  ;;;趋势一致性

    sst_p_rc = (/sst_p_rc*p_rc/)
    ; sst_p_rc = (/sst_p_rc*p_rc/sst_rc/)
    ; sst_p_rc = where(prob_sst.le.0.05, sst_p_rc, sst_p_rc@_FillValue)

  ;;;;
    

;;;;
    
print("==========finish preprocessing data =============")

;==============================================================
; Mission 4. draw plot
; ==============================================================
     ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
    ;;;

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.4
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2

  if i.eq.3
   res@lbLabelBarOn=True
   res@lbOrientation = "Horizontal"
   res@pmLabelBarWidthF = 0.5
   res@pmLabelBarHeightF = 0.1
   res@pmLabelBarOrthogonalPosF = 0.3
  end if


  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2 

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0



  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
  
  res@gsnLeftString = "(a)"+season(i)+"-observed"+" SST linear trend"
  plot(i*2) = gsn_csm_contour_map(wks,sst_rc,res)
  res@cnMinLevelValF       =   -0.5
  res@cnMaxLevelValF       =   0.5
  res@cnLevelSpacingF      =  0.1
  ;   res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2

  res@gsnLeftString = "(b)"+season(i)+"-SST linearly congruent with GPCP Pr"
    
  plot(i*2+1)  = gsn_csm_contour_map(wks,sst_p_rc,res)

  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   
   plot@$str(i*2)$ = gsn_add_polyline(wks, plot(i*2), pgx_pool, pgy_pool, pgres)
   plot@$str(i*2+1)$ = gsn_add_polyline(wks, plot(i*2+1), pgx_pool, pgy_pool, pgres)

    ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True


    iplot(i*2)=gsn_csm_contour(wks,prob_sst,res_pro)
    iplot(i*2+1)=gsn_csm_contour(wks,prob_sst_p,res_pro)

    
    overlay(plot(i*2), iplot(i*2))
    overlay(plot(i*2+1), iplot(i*2+1))
    
   delete([/sst_rc,sst_p_rc/])
   delete([/lat1,lat2,lon1,lon2/])
end do

pres = True
pres@gsnPanelLabelBar = False
; pres@lbOrientation  = "Vertical"
gsn_panel(wks, plot, (/4,2/), pres)
  


end
    
  







    
    
;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 601
    endmon = 11012
    startyear = 6
    endyear = 110

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-model/CESM_CTRL/F2000_F19_CAM4_CTRL/pre/"
    fils1=systemfunc("ls "+diri+"CTRL.cam.h1.0101-11012_chazhi.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

  
    lat = f1->lat
    lon = f1->lon

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

;;;;U,V
    u_200  = f1->U(timeind,{200.0},:,:)
    u_850  = f1->U(timeind,{850.0},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

    v_200  = f1->V(timeind,{200},:,:)
    v_850  = f1->V(timeind,{850},:,:)

    v_850 = (/rmMonAnnCycTLL(v_850)/)
    v_200 = (/rmMonAnnCycTLL(v_200)/)
    ; v_200  = short2flt(f2->v(timeind2,{200},:,:))
    ; v_850  = short2flt(f2->v(timeind2,{850},:,:))
    ; v_850 = (/rmMonAnnCycTLL(v_850)/)
    ; v_200 = (/rmMonAnnCycTLL(v_200)/)

    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)
;;;



;;;; 
    precc = f1->PRECC(timeind,:,:)
    precl = f1->PRECL(timeind,:,:)
    
    pre = precc
    pre = (/precc + precl/)
    pre  =   (/pre*1000*3600*24/)  ;; m/s conver to mm/day
    pre = (/rmMonAnnCycTLL(pre)/)
    
    pre = (/dtrend_msg_n(ispan(1,dimsizes(pre&time),1),pre,False,True,0)/)
    
    ssta = pre

;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_850(0:dimsizes(u_850&time)-13,1,1)

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(0:dimsizes(u_850&time)-13,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(0:dimsizes(u_850&time)-13,{5:15},{90:130}), (/1,2/))/)

    ; DMI_ctrl = (/dim_avg_n_Wrap(u_850(0:dimsizes(u_850&time)-13,{5:15},{40:80}), (/1,2/)) - \   ;;;ISM 
    ;               dim_avg_n_Wrap(u_850(0:dimsizes(u_850&time)-13,{20:30},{70:90}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  
;;;


;;;strong and weak monsoon year 
    total_year = ispan(startyear,endyear-1,1)
    std = dim_stddev_Wrap(DMI_ctrl_JJA)

    inds = ind(DMI_ctrl_JJA.ge.std)
    indw = ind(DMI_ctrl_JJA.le.-1*std)

    indnor = ind(DMI_ctrl_JJA.gt.-0.5*std.and.DMI_ctrl_JJA.lt.0.5*std)

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
    
    ; ;;;remove TBO
    ; monstryear2 = monstryear
    ; monwekyear2 = monwekyear

    ; do i = 0,dimsizes(monwekyear)-1
    ;   monstryear = where(monstryear.eq.monwekyear(i)+1,0,monstryear)
    ;   monstryear = where(monstryear.eq.monwekyear(i)-1,0,monstryear)
    ; end do 

    ; ; ind_str_notbo = ind(monstryear.ne.0)
    ; ind_str_notbo = ind(monstryear.eq.0)
    ; monstryear_notbo = monstryear2(ind_str_notbo)

    ; do i = 0,dimsizes(monstryear)-1
    ;   monwekyear = where(monwekyear.eq.monstryear2(i)+1,0,monwekyear)
    ;   monwekyear = where(monwekyear.eq.monstryear2(i)-1,0,monwekyear)
    ; end do 

    ; ; ind_weak_notbo = ind(monwekyear.ne.0)
    ; ind_weak_notbo = ind(monwekyear.eq.0)
    ; monwekyear_notbo = monwekyear2(ind_weak_notbo)

;;
print(monstryear)
print(monwekyear)
; print(monstryear_notbo)
; print(monwekyear_notbo)


; exit()
  
  yearnum = dimsizes(u_850)/12 -1
  dim_u = dimsizes(u_850)



; ;;;split to 12 month nino34_12 is a array size (12,yearnum)


;   ssta_12 = new((/12,yearnum,ssta/),float)
 
;   do i=0,11   
;     if i.le.6 then 
;      nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
;      else
;      nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
;     end if
;   end do 
; ;;;;


; ;;;;split u,v to 12 month 
;   yearnum = dimsizes(nino34)/12 -1
;   nino34_12 = new((/12,yearnum/),float)
 
;   do i=0,11   
;     if i.le.6 then 
;      nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
;      else
;      nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
;     end if
;   end do 


   sst_year = ispan(startyear,endyear,1) 

   ind_ts1 =  get1Dindex(sst_year,monstryear)
   ind_tw1 =  get1Dindex(sst_year,monwekyear)
   ind_tnor1 =  get1Dindex(sst_year,noryear)
  
   u_year = ispan(startyear,endyear,1) 

   ind_us1 =  get1Dindex(u_year,monstryear)
   ind_uw1 =  get1Dindex(u_year,monwekyear)
   ind_unor1 =  get1Dindex(u_year,noryear)
  
   ;  sst_year = ispan(startyear,endyear,1) 

   ; ind_ts1 =  get1Dindex(sst_year,monstryear_notbo)
   ; ind_tw1 =  get1Dindex(sst_year,monwekyear_notbo)
   ; ind_tnor1 =  get1Dindex(sst_year,noryear)
  
   ; u_year = ispan(startyear,endyear,1) 

   ; ind_us1 =  get1Dindex(u_year,monstryear_notbo)
   ; ind_uw1 =  get1Dindex(u_year,monwekyear_notbo)
   ; ind_unor1 =  get1Dindex(u_year,noryear)

  plot = new(12,graphic)
  plot1 = new(12,graphic)
  ; wks = gsn_open_wks("pdf", "WNP-UVSST-notbo_strong-monsoon-after")
  ; wks1 = gsn_open_wks("pdf", "WNP-UVSST-notbo_weak-monsoon-after")

   wks = gsn_open_wks("pdf", "CAM4-WNP-UVSST_strong-monsoon-after")
  wks1 = gsn_open_wks("pdf", "CAM4-WNP-UVSST_weak-monsoon-after")

 month = (/"Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan(+1)","Feb(+1)","Mar(+1)","Apr(+1)","May(+1)"/) 
do i=0,11
;;; strong and weak monsoon year ssta , u ,v 
   if i.le.6 then 
      ssta_month =  ssta(i+5:dimsizes(ssta&time)-1-12:12,:,:)
      u_850_month = u_850(i+5:dimsizes(u_850&time)-1-12:12,:,:)
      v_850_month = v_850(i+5:dimsizes(v_850&time)-1-12:12,:,:)
    else
      ssta_month =  ssta(i+5:dimsizes(ssta&time)-1:12,:,:)
      u_850_month = u_850(i+5:dimsizes(u_850&time)-1:12,:,:)
      v_850_month = v_850(i+5:dimsizes(v_850&time)-1:12,:,:)
    end if

    ssta_month_str = dim_avg_n_Wrap(ssta_month(ind_ts1,:,:),0)

    u_month_str = dim_avg_n_Wrap(u_850_month(ind_us1,:,:),0)

    v_month_str = dim_avg_n_Wrap(v_850_month(ind_us1,:,:),0)

    
    ssta_month_weak = dim_avg_n_Wrap(ssta_month(ind_tw1,:,:),0)

    u_month_weak = dim_avg_n_Wrap(u_850_month(ind_uw1,:,:),0)

    v_month_weak = dim_avg_n_Wrap(v_850_month(ind_uw1,:,:),0)


;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@cnFillPalette = "CBR_drywet"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@mpFillOn = False

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; ; res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False
  
  res@gsnLeftString = "Pre&850UV " + month(i)
  res@gsnLeftStringFontHeightF = 0.02

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

   res@gsnScalarContour        = True

    ;;;;;;;LineArrow
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.4
    res@vcRefMagnitudeF          = 4.0
    res@vcRefLengthF             = 0.1     ; define length of v
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF = 1.0
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.006
  
 
   plot(i) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, ssta_month_str, res)
   plot1(i) = gsn_csm_vector_scalar_map(wks1,u_month_weak, v_month_weak, ssta_month_weak, res)
  

   ; str1 = unique_string("dum1")  ; "unique_string" will return a unique    

   ; ;;;打点
   ;  res_pro=True
   ;  res_pro@gsnAddCyclic = False
   ;  res_pro@gsnDraw=False
   ;  res_pro@gsnFrame=False
   ;  res_pro@cnFillOn = True

   ;  res_pro@cnLevelSelectionMode = "ExplicitLevels"
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  res_pro@lbLabelBarOn = False
   ;  res_pro@cnLinesOn=False
   ;  ; res_pro@cnMonoLineColor =False
   ;  ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

   ;  res_pro@cnMonoFillPattern = True
   ;  res_pro@cnFillPattern = 17
   ;  res_pro@cnFillDotSizeF = 0.001
   ;  res_pro@cnMonoFillColor = False
   ;  res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
   ;  ; res_pro@cnLevels = ispan(-6,6,1)
   ;  res_pro@cnLineLabelsOn=False
   ;  res_pro@tmYRMode = "Automatic"
   ;  res_pro@cnInfoLabelOn=False

   ;   iplot = new(4,graphic)
   ;   iplot(0) = gsn_csm_contour(wks,prob_mam  ,res_pro)
   ;   iplot(1) = gsn_csm_contour(wks,prob_jja  ,res_pro)
   ;   iplot(2) = gsn_csm_contour(wks,prob_son  ,res_pro)
   ;   iplot(3) = gsn_csm_contour(wks,prob_djf  ,res_pro)
   ;  ;  
   ;  do i=0,3
   ;   overlay(plot(i), iplot(i))
   ;  end do
end do
pres = True
pres@gsnPanelLabelBar  = True

; pres@pmLabelBarOrthogonalPosF = -0.008

; pres@lbLabelStride=1
gsn_panel(wks, plot, (/4,3/), pres)
gsn_panel(wks1, plot1, (/4,3/), pres)

end
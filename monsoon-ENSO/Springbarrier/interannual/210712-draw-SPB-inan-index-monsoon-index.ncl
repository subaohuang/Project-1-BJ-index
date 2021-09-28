begin 

  startmon = 197901
  endmon = 201912
  startyear = 1979
  endyear = 2019

;;;读入ERA wind data
    fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    ; fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1900-2012
    startmon = startmon 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  

 
;;;;U,V

    u_200  = short2flt(f1->u(timeind1,{200},:,:))
    u_850  = short2flt(f1->u(timeind1,{1000},:,:))

    ; u_200  = f1->uwnd(timeind1,{200},:,:)
    ; u_850  = f1->uwnd(timeind1,{850},:,:)

    ; u_850 = (/rmMonAnnCycTLL(u_850)/)
    ; u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

   ; fils2 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
   ; f2 = addfile(fils2, "r")
   ; time2 = f2->time
   ; date2 = cd_calendar(time2, 1)
   ; timeind2 = ind(date2.ge.197901.and.date2.le.201912)
;;;;

;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_850(:,1,1)

    ; DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(:,{22.5:32.5},{110:140}), (/1,2/)) - \
    ;               dim_avg_n_Wrap(u_850(:,{5:15},{90:130}), (/1,2/))/)

    ; DMI_ctrl = (/dim_avg_n_Wrap(u_850(:,{5:15},{40:80}), (/1,2/)) - \   ;;;ISM 
    ;               dim_avg_n_Wrap(u_850(:,{20:30},{70:90}), (/1,2/))/)

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(:,{-10:10},{130:180}), (/1,2/))/)        ;;;; equtorial wind 

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;1979-2019

    DMI_ctrl_JJA  =dim_standardize(DMI_ctrl_JJA, 1)



;;;

;;;读入SPB 指数数据

   opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual_1980-2019_method1.ascii" 
  ;opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual_1980-2019_method2.ascii" 
   spb_year = asciiread(opt1, -1, "float") ;;;1980
  ; opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual_1901-2019.ascii"     
  ; spb_year = asciiread(opt1, 112, "float") 
;;;
 
  yearnum  = endyear - startyear 


r = escorc_n(DMI_ctrl_JJA(0:dimsizes(DMI_ctrl_JJA)-2), spb_year, 0, 0)

r2 = escorc_n(DMI_ctrl_JJA(1:dimsizes(DMI_ctrl_JJA)-1), spb_year, 0, 0)

print("r="+r)
print("r2="+r2)
; exit()
;;; plotr 
;;------------------------
;;------------------------
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.9
  res@vpHeightF = 0.4
  res@vpXF      = 0.06
  res@vpYF      = 0.75


  res@tiXAxisOn      = False
  res@tiYAxisOn      = False
  res@tmXTLabelsOn   = False
  res@tmXTOn         = False

  res@trXMinF        = 0
  res@trXMaxF        = yearnum

  res@trYMinF        = -4.0
  res@trYMaxF        = 5.0


  res@tmXBLabelFontHeightF = 0.017
  res@tmYLLabelFontHeightF = 0.017
  res@gsnLeftStringFontHeightF   = 0.023
  res@gsnRightStringFontHeightF  = 0.023

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(1, yearnum, 5)
  res@tmXBLabels = ispan(startyear+1, endyear, 5)

  res@tmXBMinorOn = True
  res@tmXBMinorValues = ispan(1, yearnum, 1)
  res@tmYLMode = "Explicit"
  res@tmYLValues = fspan(-5.0, 5.0, 11)
  res@tmYLLabels = fspan(-5.0, 5.0, 11)
  res@tmYLMinorOn = True
  res@tmYLMinorValues = fspan(-5.0, 5.0, 21)
  bres = res


  res@xyMarkLineModes = (/"MarkLines","Lines"/)
  res@xyMarkers = (/1,1/)
  res@xyMarkerSizes = (/0.03,0.001/)
  res@xyMarkerColors= (/"black","black"/)
  res@xyMarkerThicknesses = (/4.0,2.0/)
  res@xyLineThicknesses  = (/4.0,4.0/)
  res@xyLineColors       = (/"black","black"/)
  res@xyDashPatterns     = (/0,4/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = "(a) Interannual SPB index Time Series"


  ;;------------------------ FIG.2 ----------------------------
  pltType = "pdf"
  ; pltName = "/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB-interannual-index&monsoon-before-1901-2012"
  pltName = "/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB-interannual-index&monsoon-before"
  wks = gsn_open_wks(pltType, pltName)
  res@xyLineColors  = (/"black","black"/)
  res@xyMarkerColors= (/"black","black"/)

  plot  = gsn_csm_xy (wks, ispan(1, yearnum, 1), spb_year, res)

 

    ;;---------------- 参考线 --------------
  x=ispan(1, yearnum, 1)
    y= x*0.1
    y = 0.0
   lnres   =    True
   lnres@gsLineColor="gray30"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=0
   lnres@gsnFrame = False
   dum000 = gsn_add_polyline(wks,plot,x,y,lnres)
   y= 0.5
   lnres@gsLineColor="red"
   lnres@gsLineDashPattern=11
   dum111 = gsn_add_polyline(wks,plot,x,y,lnres)
   y=-0.5
   lnres@gsLineColor="blue"
   dum222 = gsn_add_polyline(wks,plot,x,y,lnres)

   ; lnres@gsLineDashPattern=0
   ; lnres@gsLineColor="blue"
   ; ; dum223 = gsn_add_polyline(wks,plot,x,DMI_ctrl_JJA(1:),lnres)
   ; dum223 = gsn_add_polyline(wks,plot,ispan(0, yearnum-1, 1),DMI_ctrl_JJA(0:dimsizes(DMI_ctrl_JJA)-2),lnres)




  ;;---------------- 标题 --------------
  txres               = True
  txres@txFontHeightF = 0.019
  txres@txFont = "helvetica-bold"
  txres@txFontColor   = "black"
  ; gsn_text_ndc(wks, "Spring prediction barrier index",  0.5, 0.775, txres)
  
  draw(plot)
  frame(wks)

print("****done******")
;;;;

end 






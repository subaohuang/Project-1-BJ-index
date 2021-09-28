;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 187101
    endmon = 201912
    startyear = 1871
    endyear = 2019

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 
    
    std_nino = dim_stddev(nino34_djf)
    sst_total_year = ispan(startyear, endyear, 1)
    nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    ssta = (/center_finite_diff_n(ssta, 1, False, 0, 0)/)

;;;


;;;读入ERA wind data

    ; fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1900-2019
    startmon = startmon 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  

 
;;;;U,V

    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_200  = f1->uwnd(timeind1,{200},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

   ; fils2 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
   ; f2 = addfile(fils2, "r")
   ; time2 = f2->time
   ; date2 = cd_calendar(time2, 1)
   ; timeind2 = ind(date2.ge.197901.and.date2.le.201912)
;;;;


;;;
    fils2 = "/home/yangsong3/data-observation/20thCRV2/vwnd.mon.mean.nc"
    f2 = addfile(fils2, "r")
    time2 = f2->time
    date2 = cd_calendar(time2, 1)

    timeind2 = ind(date2.ge.startmon.and.date2.le.endmon)  ;;;风场数据从1979年开始

;;;;U,V

    v_200  = f2->vwnd(timeind1,{200},:,:)
    v_850  = f2->vwnd(timeind1,{850},:,:)

    v_850 = (/rmMonAnnCycTLL(v_850)/)
    v_200 = (/rmMonAnnCycTLL(v_200)/)
    ; v_200  = short2flt(f2->v(timeind2,{200},:,:))
    ; v_850  = short2flt(f2->v(timeind2,{850},:,:))
    ; v_850 = (/rmMonAnnCycTLL(v_850)/)
    ; v_200 = (/rmMonAnnCycTLL(v_200)/)

    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)
;;;


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

    ; ind_str_notbo = ind(monstryear.ne.0)
    ; ; ind_str_notbo = ind(monstryear.eq.0)
    ; monstryear_notbo = monstryear2(ind_str_notbo)

    ; do i = 0,dimsizes(monstryear)-1
    ;   monwekyear = where(monwekyear.eq.monstryear2(i)+1,0,monwekyear)
    ;   monwekyear = where(monwekyear.eq.monstryear2(i)-1,0,monwekyear)
    ; end do 

    ; ind_weak_notbo = ind(monwekyear.ne.0)
    ; ; ind_weak_notbo = ind(monwekyear.eq.0)
    ; monwekyear_notbo = monwekyear2(ind_weak_notbo)

  ;;;select no ENSO effect year 
    monstryear2 = monstryear
    monwekyear2 = monwekyear

    do i = 0,dimsizes(nino_year)-1
      monstryear = where(monstryear.eq.nino_year(i),0,monstryear)
    end do 
    
    do i = 0,dimsizes(lanina_year)-1
      monstryear = where(monstryear.eq.lanina_year(i),0,monstryear)
    end do 

    ind_str_notbo = ind(monstryear.ne.0)
    ; ind_str_notbo = ind(monstryear.eq.0)
    monstryear_notbo = monstryear2(ind_str_notbo)

    do i = 0,dimsizes(lanina_year)-1
      monwekyear = where(monwekyear.eq.lanina_year(i),0,monwekyear)
    end do 
     
    do i = 0,dimsizes(nino_year)-1
      monwekyear = where(monwekyear.eq.nino_year(i),0,monwekyear)
    end do 


    ind_weak_notbo = ind(monwekyear.ne.0)
    ; ind_weak_notbo = ind(monwekyear.eq.0)
    monwekyear_notbo = monwekyear2(ind_weak_notbo)


;;
; print(monstryear)
; print(monwekyear)
; exit()

  yearnum = dimsizes(nino34)/12 -1
  dim_sst = dimsizes(ssta)
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


   ; sst_year = ispan(startyear,endyear,1) 

   ; ind_ts1 =  get1Dindex(sst_year,monstryear2)
   ; ind_tw1 =  get1Dindex(sst_year,monwekyear2)
   ; ind_tnor1 =  get1Dindex(sst_year,noryear)
  
   ; u_year = ispan(startyear,endyear,1) 

   ; ind_us1 =  get1Dindex(u_year,monstryear2)
   ; ind_uw1 =  get1Dindex(u_year,monwekyear2)
   ; ind_unor1 =  get1Dindex(u_year,noryear)
  
    sst_year = ispan(startyear,endyear,1) 

   ind_ts1 =  get1Dindex(sst_year,monstryear_notbo)
   ind_tw1 =  get1Dindex(sst_year,monwekyear_notbo)
   ind_tnor1 =  get1Dindex(sst_year,noryear)
  
   u_year = ispan(startyear,endyear,1) 

   ind_us1 =  get1Dindex(u_year,monstryear_notbo)
   ind_uw1 =  get1Dindex(u_year,monwekyear_notbo)
   ind_unor1 =  get1Dindex(u_year,noryear)

  plot = new(2,graphic)
  
  wks = gsn_open_wks("pdf", "SST-tendency-nocoenso_strong-monsoon-after")
  ; wks = gsn_open_wks("pdf", "SST-tendency-strong-monsoon-after")

 month = (/"J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/) 

 ssta_month_str = new(15,typeof(ssta))
 ssta_month_weak = ssta_month_str
do i=0,14
;;; strong and weak monsoon year ssta , u ,v 
   if i.le.6 then 
      ssta_month =  ssta(i+5:dimsizes(ssta&time)-1-12:12,:,:)
    else
      ssta_month =  ssta(i+5:dimsizes(ssta&time)-1:12,:,:)
    end if

    ssta_month_str(i) = avg(dim_avg_n_Wrap(ssta_month(ind_ts1,{-5:5},{190:240}),0))

    ssta_month_weak(i) = avg(dim_avg_n_Wrap(ssta_month(ind_tw1,{-5:5},{190:240}),0))
end do
;;**************************************************************
;画图
;;***************************************************************

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    ; res@tiMainString  = "15-years window running std for JFMAM MSE "   ; title
    ; res@tiMainString  = " SSTA tendency evolution"   ; title

    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02

    ; res@tiYAxisString = "Std"           ; y axis title
    res@tiYAxisString = "SSTA tendency"           ; y axis title
    res@tiXAxisString = "Month"

    
    res@xyLineColor      = (/"black"/)  ; line colors
   

    res@xyDashPattern    = (/0.0/)        ; line patterns
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    res@gsnYRefLine = 0.0

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 14
    res@trXMinF = 0
    ; res@trYMinF = 0.0
    ; res@trYMaxF = 2.0

    res@trYMinF = -0.4
    res@trYMaxF = 0.4

    res@tmXTOn = False
    res@tmYROn = False

     
    res@vpHeightF = 0.4
    res@vpWidthF = 0.7

   

    res@tmXBValues=ispan(0,14,1)
    res@tmXBLabels=month
 
    res@tmYRMode="Automatic"

    res@gsnLeftString = "SSTA tendency after JJA Strong monsoon"
    
    plot(0) = gsn_csm_xy(wks, ispan(0, 14, 1), (/ssta_month_str/), res)

    res@gsnLeftString = "SSTA tendency after JJA Weak monsoon"

    plot(1) = gsn_csm_xy(wks, ispan(0, 14, 1), (/ssta_month_weak/), res)

  
    ; ;;添加legend 
    ; gres = True
    ; gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;       ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ; gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ; lineres = True
    ; ; lineres@lgLineColors = (/"black"，/) ; line colors
    ; lineres@lgLineColors       = res@xyLineColors ; colors chosen
    ; lineres@lgLineThicknesses = res@xyLineThicknessF                       ; line thicknesses
    ; lineres@LineLengthPercent = 9.
    ; lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    ; textres = True

    ; textres@lgLabels = (/""/)
    ; textres@lgLabelFontHeightF = 0.01
    ; ; textres@lgLabelOpacity = (/0/) 
   
    ; plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   

  pres = True 
  gsn_panel(wks, plot, (/2,1/), pres)
end
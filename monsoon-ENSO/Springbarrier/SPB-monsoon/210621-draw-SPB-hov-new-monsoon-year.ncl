;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
function spb_cal(data,data2)  ;;  data:two dimension 12×yearnum array i.e.12×10
;;data1 : the monsoon anomous year  data2: next year
begin 
  dim = dimsizes(data)
  yearnum = dim(1)
  ;;; cal ACF, want to get a two dimensional array 12×13,initial month start June 
     acr = new((/12,13/),float)
     do i=0,11
         do j=0,12
            k=i+j 
            if k.le.11 then 
              acr(i,j) = escorc(data(i,0:yearnum-1),data(i+j,0:yearnum-1))
            else 
              acr(i,j) = escorc(data(i,0:yearnum-1),data2(i+j-12,0:yearnum-1))
            end if 
         end do 
     end do
   ;;;
 ;;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 2, False, 1, 1)

   acr_gra = -1*acr_gra  ;;;change to positive
     print(acr_gra(0,:))
   exit()
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"
   maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   SPB = dim_sum_n_Wrap(acr_gra_max(6:), 0)
   acr@maind = maind
   acr@SPB  = SPB
return acr
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 190001
    endmon = 201912
    startyear = 1900
    endyear = 2019
;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

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
    ; nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{210:270}), (/1,2/))
;;;



;;;读入ERA wind data

    ; fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon - 100
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

;;;;U,V
    ; v_200  = short2flt(f2->v(timeind2,{200},:,:))
    ; v_850  = short2flt(f2->v(timeind2,{850},:,:))
    ; v_850 = (/rmMonAnnCycTLL(v_850)/)
    ; v_200 = (/rmMonAnnCycTLL(v_200)/)

     u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
;;;


;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_850(:,1,1)

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(:,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(:,{5:15},{90:130}), (/1,2/))/)


    ; DMI_ctrl = (/dim_avg_n_Wrap(u_850(:,{5:15},{40:80}), (/1,2/)) - \
    ;               dim_avg_n_Wrap(u_850(:,{20:30},{70:90}), (/1,2/))/)


    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;;jja 少一年，从25601开始
;;;

;;;strong and weak monsoon year 
   total_year = ispan(startyear,endyear-1,1)
   std = dim_stddev_Wrap(DMI_ctrl_JJA)

   inds = ind(DMI_ctrl_JJA.ge.std)
   indw = ind(DMI_ctrl_JJA.le.-1*std)

   indnor = ind(DMI_ctrl_JJA.gt.-0.5*std.and.DMI_ctrl_JJA.lt.0.5*std)

   monstryear = total_year(inds)
   monwekyear = total_year(indw)

; ;;;ISM index  850hPA ZONAL wind diff 5-15N 40-80E - 20-30N,70-90E
;     ISM = u_850(:,1,1)

;     ISM = (/dim_avg_n_Wrap(u_850(:,{5:15},{40:80}), (/1,2/)) - \
;                   dim_avg_n_Wrap(u_850(:,{20:30},{70:90}), (/1,2/))/)

;     ISM_JJA = month_to_season(ISM, "JJA")  ;;;jja 少一年，从25601开始
; ;;;;


; ;;;strong and weak monsoon year 
;    total_year = ispan(startyear,endyear-1,1)
;    std = dim_stddev_Wrap(ISM_JJA)

;    inds = ind(ISM_JJA.ge.std)
;    indw = ind(ISM_JJA.le.-1*std)

;    indnor = ind(ISM_JJA.gt.-0.5*std.and.ISM_JJA.lt.0.5*std)

;    monstryear = total_year(inds)
;    monwekyear = total_year(indw)
; ; print(monstryear)
; ; print(monwekyear)
; ; exit()

;;;;;  
  
;;;split to 12 month nino34_12 is a array size (12,yearnum)
  yearnum = dimsizes(nino34)/12 -1
  nino34_12 = new((/12,yearnum/),float)

  do i=0,11   
    if i.le.9 then 
     nino34_12(i,:) = nino34(i+2:dimsizes(nino34)-1-12:12)
     else
     nino34_12(i,:) = nino34(i+2:dimsizes(nino34)-1:12)
    end if
  end do 

 ; printVarSummary(nino34_12)


;;; strong and weak monsoon year ssta 

  nino34_12_strmon = nino34_12(:,inds)
  nino34_12_weakmon = nino34_12(:,indw)
  nino34_12_normon = nino34_12(:,indnor)

  nino34_12_strmon2 = nino34_12(:,inds+1)
  nino34_12_weakmon2 = nino34_12(:,indw+1)
  nino34_12_normon2 = nino34_12(:,indnor+1)

  ; spb_strmon = spb_cal(nino34_12_strmon,nino34_12_strmon2)
  spb_weakmon = spb_cal(nino34_12_weakmon,nino34_12_weakmon2)
  spb_normon = spb_cal(nino34_12_normon,nino34_12_normon2)
 
  ; print("spb_strmon="+spb_strmon@SPB)
  print("spb_weakmon="+spb_weakmon@SPB)
  print("spb_normon="+spb_normon@SPB)

  plot = new(3,graphic)

;;;;


;;;draw root mean square R(lagt) =  sqrt(((1/N)*sum(ACF(lagt)²))

   ;;分成3-8月 ：before and during monsoon   9-2月： after monsoon
   N = 8-3+1

   acf_square_strmon = spb_strmon^2
   acf_square_weakmon = spb_weakmon^2
   acf_square_normon = spb_normon^2

   R38_strmon =  sqrt(dim_sum_n(acf_square_strmon(0:5,:),0)/N)

   R38_weakmon =  sqrt(dim_sum_n(acf_square_weakmon(0:5,:),0)/N)

   R38_normon =  sqrt(dim_sum_n(acf_square_normon(0:5,:),0)/N)

   R92_strmon =  sqrt(dim_sum_n(acf_square_strmon(6:,:),0)/N)

   R92_weakmon =  sqrt(dim_sum_n(acf_square_weakmon(6:,:),0)/N)

   R92_normon =  sqrt(dim_sum_n(acf_square_normon(6:,:),0)/N)

   
   R38_strmon = center_finite_diff(R38_strmon, 2, False, 1)*-10
   R38_weakmon = center_finite_diff(R38_weakmon, 2, False, 1)*-10
   R38_normon = center_finite_diff(R38_normon, 2, False, 1)*-10

   R92_strmon = center_finite_diff(R92_strmon, 2, False, 1)*-10
   R92_weakmon = center_finite_diff(R92_weakmon, 2, False, 1)*-10
   R92_normon = center_finite_diff(R92_normon, 2, False, 1)*-10
   
 
   





   ; wks = gsn_open_wks("pdf", "ISM-acr_hov-newmonsoon-strong-weak-monsoon")
   ; wks1 = gsn_open_wks("pdf", "ISM-root_mean_square-newmonsoon")


   wks = gsn_open_wks("pdf", "WNP-acr_hov-newmonsoon-strong-weak-monsoon")
   ; wks1 = gsn_open_wks("pdf", "WNP-root_mean_square-newmonsoon")

   wks1 = gsn_open_wks("pdf", "central-diff-WNP-root_mean_square-newmonsoon")

;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=True
    res@gsnAddCyclic = False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn= True
    ; res@pmLabelBarOrthogonalPosF = 0.2
    ; res@lbLabelStride=1
    res@lbOrientation="Vertical"
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmBorderThicknessF=1.0
    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF = 1.5

    res@tmXBValues = ispan(0,12,1)
    res@tmXBLabels = ispan(0,12,1) +""

    res@tmYLValues = ispan(0,11,1)
    res@tmYLLabels = (/"Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb"/) 

    ; res@cnFillPalette = "BlueDarkOrange18"
    gsn_define_colormap(wks,"BlueDarkRed18")


    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels

    res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 
  ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

   res@gsnLeftString = "strongmonsoon"
   plot(0) = gsn_csm_contour(wks,spb_strmon,res)
   res@gsnLeftString = "weakmonsoon"
   plot(1) = gsn_csm_contour(wks,spb_weakmon,res)
   res@gsnLeftString = "normalmonsoon"
   plot(2) = gsn_csm_contour(wks,spb_normon,res)

   resp = True
   resp@gsMarkerColor = "black"
   resp@gsMarkerSizeF = 0.1
   resp@gsMarkerIndex = 1

   plot_gs = gsn_add_polymarker(wks, plot(0), spb_strmon@maind, ispan(0,11,1), resp)
   plot_gs1 = gsn_add_polymarker(wks, plot(1), spb_weakmon@maind, ispan(0,11,1), resp)
   plot_gs2 = gsn_add_polymarker(wks, plot(2), spb_normon@maind, ispan(0,11,1), resp)


   pres = True

   gsn_panel(wks, plot, (/2,2/), pres)



   ;;;;;;; draw root_mean_square

      resxy = True
      resxy@gsnDraw = False
      resxy@gsnFrame = False
      resxy@trXMinF = 0
      resxy@trXMaxF = 12

      ; resxy@trYMaxF = 1.0
      ; resxy@trYMinF = 0.0
      
      resxy@xyMonoLineColor = False
      resxy@xyLineColors = (/"red","blue","black"/)
      resxy@xyDashPatterns= (/0,0,0/)

      resxy@gsnLeftString = "Rb"
      resxy@vpHeightF = 0.4
      resxy@vpWidthF = 0.6
      plotxy = new(2,graphic)
      plotxy(0) = gsn_csm_xy(wks1, ispan(0,12,1), (/R38_strmon,R38_weakmon,R38_normon/) , resxy) 
      resxy@gsnLeftString = "Ra"
      plotxy(1) = gsn_csm_xy(wks1, ispan(0,12,1), (/R92_strmon,R92_weakmon,R92_normon/) , resxy) 

        ;;添加legend 
        gres = True
        gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
        gres@XPosPercent = 68    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

        lineres = True
        ; lineres@lgLineColors = (/"black"，/) ; line colors
        lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
        lineres@lgLineThicknesses = 1.0                       ; line thicknesses
        lineres@LineLengthPercent = 9.
        lineres@lgDashIndexes   = resxy@xyDashPatterns                     ; expressed as %, 0->100, length of line

        textres = True

        textres@lgLabels = (/"Strong Monsoon","Weak Monsoon","Normal Monsoon"/)
        textres@lgLabelFontHeightF = 0.02
        ; textres@lgLabelOpacity = (/0/) 

        plotlg = simple_legend(wks1,plotxy(0),gres,lineres,textres)   
        plotlg1 = simple_legend(wks1,plotxy(1),gres,lineres,textres)   

  
   pres = True
   gsn_panel(wks1, plotxy, (/2,1/), pres)



   ; draw(plot)
   ; frame(wks)

   end 

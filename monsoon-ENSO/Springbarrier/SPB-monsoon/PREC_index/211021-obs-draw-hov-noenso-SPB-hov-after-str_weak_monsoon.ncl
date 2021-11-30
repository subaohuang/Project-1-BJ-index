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
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"
   maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   SPB = dim_sum_n_Wrap(acr_gra_max, 0)
   acr@maind = maind
   acr@SPB = SPB
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

    nino3_4 = month_to_season(nino34, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)

    nino3_4(0:dimsizes(nino3_4)-2) = (/nino3_4(1:)/)
    nino3_4(dimsizes(nino3_4)-1) =nino34(dimsizes(nino3_4)*12-1) 

   
    std_sst = dim_stddev_n(nino3_4, 0)
    

    year_sst = ispan(startyear,endyear,1)

    enino_year = year_sst(ind(nino3_4.ge.1*std_sst))
    lanina_year = year_sst(ind(nino3_4.le.-1*std_sst))

    

;;;


;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:2016}) 

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2016-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;
;;;

;;;strong and weak monsoon year 
   total_year = ispan(startyear,endyear-2,1)
   std = dim_stddev_Wrap(DMI_ctrl_JJA)

   ; inds = ind(DMI_ctrl_JJA.ge.std)
   ; indw = ind(DMI_ctrl_JJA.le.-1*std)
   ; indnor = ind(DMI_ctrl_JJA.gt.-0.5*std.and.DMI_ctrl_JJA.lt.0.5*std)

    inds = ind(DMI_ctrl_JJA.ge.0.43*std)
    indw = ind(DMI_ctrl_JJA.le.-0.43*std)

    indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)


   monstryear = total_year(inds)
   monwekyear = total_year(indw)
   noryear = total_year(indnor)
; print(monstryear)
; print(monwekyear)
; exit()

  ;;;select no ENSO effect year : revised version 
    monstryear2 = monstryear
    monwekyear2 = monwekyear

    do i = 0,dimsizes(enino_year)-1
      monstryear = where(monstryear.eq.enino_year(i),0,monstryear)
    end do 
    
    do i = 0,dimsizes(lanina_year)-1
      monstryear = where(monstryear.eq.lanina_year(i),0,monstryear)
    end do 

    ind_str_notbo = ind(monstryear.ne.0)
    ; ind_str_notbo = ind(monstryear.eq.0)
    monstryear_notbo1 = monstryear2(ind_str_notbo)  ;;;;让样本数量一致

    monstryear_notbo = monstryear_notbo1(9:)

    do i = 0,dimsizes(lanina_year)-1
      monwekyear = where(monwekyear.eq.lanina_year(i),0,monwekyear)
    end do 
     
    do i = 0,dimsizes(enino_year)-1
      monwekyear = where(monwekyear.eq.enino_year(i),0,monwekyear)
    end do 


    ind_weak_notbo = ind(monwekyear.ne.0)
    ; ind_weak_notbo = ind(monwekyear.eq.0)
    monwekyear_notbo = monwekyear2(ind_weak_notbo)

    

    print(monstryear_notbo)
    print(monwekyear_notbo)
  ;;;;
   

  
;;;split to 12 month nino34_12 is a array size (12,yearnum)
  yearnum = dimsizes(nino34)/12 -1
  nino34_12 = new((/12,yearnum/),float)

  do i=0,11   
    if i.le.6 then 
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
     else
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
    end if
  end do 


   sst_year = ispan(startyear,endyear,1) ;; nino34_12 start from startyear-1 because nino34 start from startyear-2

   ind_s1 =  get1Dindex(sst_year,monstryear_notbo)
   ind_w1 =  get1Dindex(sst_year,monwekyear_notbo)
   ind_nor1 =  get1Dindex(sst_year,noryear)
 ; printVarSummary(nino34_12)


;;; strong and weak monsoon year ssta 

  nino34_12_strmon = nino34_12(:,ind_s1)
  nino34_12_weakmon = nino34_12(:,ind_w1)
  nino34_12_normon = nino34_12(:,ind_nor1)

  nino34_12_strmon2 = nino34_12(:,ind_s1+1)
  nino34_12_weakmon2 = nino34_12(:,ind_w1+1)
  nino34_12_normon2 = nino34_12(:,ind_nor1+1)

  spb_strmon = spb_cal(nino34_12_strmon,nino34_12_strmon2)
  spb_weakmon = spb_cal(nino34_12_weakmon,nino34_12_weakmon2)
  spb_normon = spb_cal(nino34_12_normon,nino34_12_normon2)
  
  ; system("rm WNP_notbo_acr_hov_after_monsoon.nc")
  ; fi_out = addfile("WNP_notbo_acr_hov_after_monsoon.nc", "c")
  ;  system("rm WNP_tbo_acr_hov_after_monsoon.nc")
  ; fi_out = addfile("WNP_tbo_acr_hov_after_monsoon.nc", "c")
  ; fi_out->spb_strmon = spb_strmon 
  ; fi_out->spb_weakmon = spb_weakmon
  ; fi_out->spb_normon = spb_normon

print("spb_strmon="+spb_strmon@SPB)
print("spb_weakmon="+spb_weakmon@SPB)




  plot = new(3,graphic)

;;;;
;;;draw root mean square R(lagt) =  sqrt(((1/N)*sum(ACF(lagt)²))

   ;;分成3-8月 ：before and during monsoon   9-2月： after monsoon
   N = 12

   acf_square_strmon = spb_strmon^2
   acf_square_weakmon = spb_weakmon^2
   acf_square_normon = spb_normon^2

   R_strmon =  sqrt(dim_sum_n(acf_square_strmon(:,:),0)/N)

   R_weakmon =  sqrt(dim_sum_n(acf_square_weakmon(:,:),0)/N)

   R_normon =  sqrt(dim_sum_n(acf_square_normon(0,:),0)/N)


; wks = gsn_open_wks("pdf", "WNP-after-acr_hov-strong-weak-monsoon")
; wks1 = gsn_open_wks("pdf", "WNP-after-root_mean_square-newmonsoon.pdf")

; wks = gsn_open_wks("pdf", "NOsimutaousenso-WNP-noenso-after-acr_hov-strong-weak-monsoon")
; wks1 = gsn_open_wks("pdf", "NOsimutaousenso-WNP-noenso-after-root_mean_square.pdf")

wks = gsn_open_wks("pdf", "NOsimutaousenso-AISMR-noenso-after-acr_hov-strong-weak-monsoon")
wks1 = gsn_open_wks("pdf", "NOsimutaousenso-AISMR-noenso-after-root_mean_square.pdf")


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
    res@tmYLLabels = (/"Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May"/) 

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
   ; draw(plot)
   ; frame(wks)
   ;;;;;;; draw root_mean_square

      resxy = True
      resxy@gsnDraw = False
      resxy@gsnFrame = False
      resxy@trXMinF = 0
      resxy@trXMaxF = 12

      resxy@trYMaxF = 1.0
      resxy@trYMinF = 0.0
      
      resxy@xyMonoLineColor = False
      resxy@xyLineColors = (/"red","blue"/)
      resxy@xyDashPatterns= (/0,1/)

      resxy@gsnLeftString = "R"
      resxy@vpHeightF = 0.4
      resxy@vpWidthF = 0.6
      plotxy = new(1,graphic)
      plotxy(0) = gsn_csm_xy(wks1, ispan(0,12,1), (/R_strmon,R_weakmon/) , resxy) 
     
    draw(plotxy)
    frame(wks)
   end 
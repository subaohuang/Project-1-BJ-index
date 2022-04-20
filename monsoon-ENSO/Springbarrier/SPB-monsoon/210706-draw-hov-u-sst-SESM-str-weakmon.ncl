function one_test(pv_ano,pv_clm,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east- pv_clm)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
 end

function corr(data,inds)  ;;  
;;data1 : the monsoon index 1901-2011 data2: 1900-2012
begin 
  
  dim1 =dimsizes(data)
  yearnum1 = dim1(1)
  acr = new((/36,dim1(2)/),float)
  r_test = new((/36,dim1(2)/),float)
  n = dimsizes(inds)



 do i=0,35
    if i.le.11 then 
         acr(i,:) = dim_avg_n_Wrap(data(i,inds-1,:),0)
         var =  dim_variance_n(data(i,inds-1,:), 0)    ; 13.477
         var = sqrt(var/n) 
         t = acr(i,:)/var
         r_test(i,:) = student_t(t, n-1)
    else 
      if i.le.23 then 
         acr(i,:) = dim_avg_n_Wrap(data(i-12,inds,:),0)
         var =  dim_variance_n(data(i-12,inds,:), 0)    ; 13.477
         var = sqrt(var/n) 
         t = acr(i,:)/var
         r_test(i,:)  = student_t(t, n-1)
      else
         acr(i,:) = dim_avg_n_Wrap(data(i-24,inds+1,:),0)
         var =  dim_variance_n(data(i-24,inds+1,:), 0)    ; 13.477
         var = sqrt(var/n) 
         t = acr(i,:)/var
         r_test(i,:)  = student_t(t, n-1)

      end if
    end if
  end do 

  acr!0 = "time"
  acr!1 = "lon"
  acr&lon = data&lon
  acr&lon@units = "degree_east"
  acr@prob = r_test
return acr
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 190001
    endmon = 201212
    startyear = 1900
    endyear = 2012
;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},:), 1)

    nino = dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
    nino3_4 = month_to_season(nino, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)

    nino3_4(0:dimsizes(nino3_4)-2) = (/nino3_4(1:)/)
    nino3_4(dimsizes(nino3_4)-1) =nino(dimsizes(nino&time)-1) 


    std_sst = dim_stddev_n(nino3_4, 0)


    year_sst = ispan(startyear,endyear,1)

    enino_year = year_sst(ind(nino3_4.ge.std_sst))
    lanina_year = year_sst(ind(nino3_4.le.-1*std_sst))

;;;


;;;读入ERA wind data

   ; fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
   fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
   f1 = addfile(fils1, "r")
   time1 = f1->time
   date1=cd_calendar(time1, 1)

   endmon = endmon 
   timeind1=ind(date1.ge.startmon.and.date1.le.endmon) ;;1900-2012
 
  
;;;;U,V
    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_1000  = f1->uwnd(timeind1,{1000},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_1000 = (/rmMonAnnCycTLL(u_1000)/)
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
     u_1000 = (/dtrend_msg_n(ispan(1,dimsizes(u_1000&time),1),u_1000,False,True,0)/)
;;;


;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_850(12:dimsizes(u_850&time)-13,1,1)

    ; DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(12:dimsizes(u_850&time)-13,{22.5:32.5},{110:140}), (/1,2/)) - \
    ;               dim_avg_n_Wrap(u_850(12:dimsizes(u_850&time)-13,{5:15},{90:130}), (/1,2/))/)
    DMI_ctrl = (/dim_avg_n_Wrap(u_850(12:dimsizes(u_850&time)-13,{5:15},{40:80}), (/1,2/)) - \   ;;;ISM 
                  dim_avg_n_Wrap(u_850(12:dimsizes(u_850&time)-13,{20:30},{70:90}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;;1901-2011

   ;; surface wind

    u_zonal = dim_avg_n_Wrap(u_1000(:,{-5:5},:),1)  ;;1900-2012
;;;

;;;strong and weak monsoon year 
   total_year = ispan(startyear+1,endyear-1,1)
   std = dim_stddev_Wrap(DMI_ctrl_JJA)

   inds = ind(DMI_ctrl_JJA.ge.std)
   indw = ind(DMI_ctrl_JJA.le.-1*std)

   indnor = ind(DMI_ctrl_JJA.gt.-0.5*std.and.DMI_ctrl_JJA.lt.0.5*std)

   monstryear = total_year(inds)
   monwekyear = total_year(indw)
   noryear = total_year(indnor)
;;;;;  

; ;;; tbo
;     monstryear2 = monstryear
;     monwekyear2 = monwekyear

;     do i = 0,dimsizes(monwekyear)-1
;       monstryear = where(monstryear.eq.monwekyear(i)+1,0,monstryear)
;       monstryear = where(monstryear.eq.monwekyear(i)-1,0,monstryear)
;     end do 

;     ; ind_str_notbo = ind(monstryear.ne.0)
;     ind_str_notbo = ind(monstryear.eq.0)
;     monstryear_notbo = monstryear2(ind_str_notbo)

;     do i = 0,dimsizes(monstryear)-1
;       monwekyear = where(monwekyear.eq.monstryear2(i)+1,0,monwekyear)
;       monwekyear = where(monwekyear.eq.monstryear2(i)-1,0,monwekyear)
;     end do 

;     ; ind_weak_notbo = ind(monwekyear.ne.0)
;     ind_weak_notbo = ind(monwekyear.eq.0)
;     monwekyear_notbo = monwekyear2(ind_weak_notbo)
; ;;;;

;;;select no ENSO effect year 
    monstryear2 = monstryear
    monwekyear2 = monwekyear

    do i = 0,dimsizes(enino_year)-1
      monstryear = where(monstryear.eq.enino_year(i)+1,0,monstryear)
    end do 
    
    do i = 0,dimsizes(lanina_year)-1
      monstryear = where(monstryear.eq.lanina_year(i),0,monstryear)
    end do 

    ind_str_notbo = ind(monstryear.ne.0)
    ; ind_str_notbo = ind(monstryear.eq.0)
    monstryear_notbo = monstryear2(ind_str_notbo)

    do i = 0,dimsizes(lanina_year)-1
      monwekyear = where(monwekyear.eq.lanina_year(i)+1,0,monwekyear)
    end do 
     
    do i = 0,dimsizes(enino_year)-1
      monwekyear = where(monwekyear.eq.enino_year(i),0,monwekyear)
    end do 


    ind_weak_notbo = ind(monwekyear.ne.0)
    ; ind_weak_notbo = ind(monwekyear.eq.0)
    monwekyear_notbo = monwekyear2(ind_weak_notbo)
;;;


;;;;
   sst_year = ispan(startyear,endyear,1) 

   ind_ts1 =  get1Dindex(sst_year,monstryear2)
   ind_tw1 =  get1Dindex(sst_year,monwekyear2)
   ind_tnor1 =  get1Dindex(sst_year,noryear)

   ; ind_ts1 =  get1Dindex(sst_year,monstryear_notbo)
   ; ind_tw1 =  get1Dindex(sst_year,monwekyear_notbo)
   ; ind_tnor1 =  get1Dindex(sst_year,noryear)
  
   ; u_year = ispan(startyear,endyear,1) 
   ; ind_us1 =  get1Dindex(u_year,monstryear_notbo)
   ; ind_uw1 =  get1Dindex(u_year,monwekyear_notbo)
   ; ind_unor1 =  get1Dindex(u_year,noryear)


;;;split to 12 month nino34_12 is a array size (12,yearnum) monsoon year 
  dim_sst = dimsizes(nino34)
  yearnum = dimsizes(nino34&time)/12
  nino34_12 = new((/12,yearnum,dim_sst(1)/),float)

  dim_u = dimsizes(u_zonal)

  u_12 = new((/12,yearnum,dim_u(1)/),float)

  do i=0,11   
     nino34_12(i,:,:) = nino34(i:dim_sst(0)-1:12,:)
     u_12(i,:,:) = u_zonal(i:dim_u(0)-1:12,:)
  end do 


  nino34_12!2 = "lon"
  nino34_12&lon = nino34&longitude
  nino34_12&lon@units = "degree_east"

  u_12!2 = "lon"
  u_12&lon = u_zonal&lon
  u_12&lon@units = "degree_east"


  printVarSummary(nino34_12)
  printVarSummary(u_12)
;;;

;;;相关系数

 corr_sst_str = corr(nino34_12(:,:,{59.5:280.5}),ind_ts1)
 corr_u1000_str = corr(u_12(:,:,{60:280}),ind_ts1)

 corr_sst_weak = corr(nino34_12(:,:,{59.5:280.5}),ind_tw1)
 corr_u1000_weak = corr(u_12(:,:,{60:280}),ind_tw1)


 prob_sst_str = corr_sst_str
 prob_sst_str = (/corr_sst_str@prob/)

 prob_u_str = corr_u1000_str
 prob_u_str = (/corr_u1000_str@prob/)
 
 prob_sst_weak = corr_sst_weak
 prob_sst_weak = (/corr_sst_weak@prob/)

 prob_u_weak = corr_u1000_weak
 prob_u_weak = (/corr_u1000_weak@prob/)

 ; corr_sst = where(prob_sst.le.0.05, corr_sst, corr_sst@_FillValue)
 ; corr_u1000 = where(prob_u.le.0.05, corr_u1000, corr_u1000@_FillValue)

;;;;

; wks = gsn_open_wks("pdf", "nocoenso-Hov-U-SST-str-weak-mon")
; wks = gsn_open_wks("pdf", "noallenso-Hov-U-SST-str-weak-mon")
wks = gsn_open_wks("pdf", "ISM-Hov-U-SST-str-weak-mon")



;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=True
    res@gsnAddCyclic = False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn= True
    ; res@pmLabelBarOrthogonalPosF = -0.05
    ; res@lbLabelStride=1
    ; res@lbOrientation="Horizontal"
    res@lbOrientation="Vertical"
    ; res@pmLabelBarHeightF = 0.05
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmBorderThicknessF=1.0
    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF = 1.5
    res@cnInfoLabelOn = False
    

    res@tmXBValues = ispan(60,280,40)
    ; res@tmXBLabels = (/"60E","80E","100E","120E","140E","160E","180","160W","140W","120W","100W","80W"/)
    res@tmXBLabels = (/"60E","100E","140E","180","140W","100W","80W"/)


      res@tmXBMinorValues = ispan(60,280,10)

    res@tmYLValues = ispan(0,35,3)
    res@tmYLMinorValues = ispan(0,35,1)
    


    ; res@tmYLLabels = (/"Jun(-1)","Sep(-1)","Nov(-1)","Mar","Jun","Sep","Nov","Mar(1)","Jun(+1)"/) 
    res@tmYLLabels = (/"Jan(-1)","Apr(-1)","Jul(-1)","Oct(-1)","Jan","Apr","Jul","Oct","Jan(+1)","Apr(+1)","Jul(+1)","Oct(+1)"/) 

    ; res@cnFillPalette = "BlueDarkOrange18"
    gsn_define_colormap(wks,"BlueDarkRed18")


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1

    ; res@cnMinLevelValF       =   -0.8
    ; res@cnMaxLevelValF       =   0.8
    ; res@cnLevelSpacingF      =  0.2

    
    ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    
    ; res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    ; res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 
    ; gsn_define_colormap(wks, "precip_diff_12lev")
    ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
    ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

    plot = new(4,graphic)
    res@gsnRightString = ""

    res@gsnLeftString = "(a) SST anomaly during Str monsoon"
    plot(0) = gsn_csm_contour(wks,corr_sst_str,res)
    res@gsnLeftString = "(b) Surface wind during Str monsoon"

    res@cnMinLevelValF       =   -2.0
    res@cnMaxLevelValF       =   2.0
    res@cnLevelSpacingF      =  0.4
   

    plot(1) = gsn_csm_contour(wks,corr_u1000_str,res)

    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1
    ;    res@cnMinLevelValF       =   -0.8
    ; res@cnMaxLevelValF       =   0.8
    ; res@cnLevelSpacingF      =  0.2

    ; res@lbLabelBarOn= True

    res@gsnLeftString = "(c) SST anomaly during Weak monsoon"
    plot(2) = gsn_csm_contour(wks,corr_sst_weak,res)
    res@gsnLeftString = "(d) Surface wind during Weak monsoon"

    res@cnMinLevelValF       =   -2.0
    res@cnMaxLevelValF       =   2.0
    res@cnLevelSpacingF      =  0.4
    
    plot(3) = gsn_csm_contour(wks,corr_u1000_weak,res)


;;;

;;;打点

    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False

    res_pro@cnLinesOn=True
    res_pro@cnMonoLineColor =False
    res_pro@cnLineColors = ((/"Transparent","green","Transparent"/))
    res_pro@cnMonoLineThickness = True
    res_pro@cnLineThicknessF = 2.

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","green","Transparent"/)
  
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    res_pro@gsnLeftString = ""
    res_pro@gsnRightString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(4,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_sst_str,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob_u_str ,res_pro)

    iplot(2)=gsn_csm_contour(wks,prob_sst_weak,res_pro)
    iplot(3)=gsn_csm_contour(wks,prob_u_weak ,res_pro)



   do i=0,3
   overlay( plot(i),iplot(i))
end do

   pres = True

   gsn_panel(wks, plot, (/2,2/), pres)


end 


   


function corr(data,data2)  ;;  
;;data1 : the monsoon index 1901-2011 data2: 1900-2012
begin 
  
  dim1 =dimsizes(data)
  dim2 = dimsizes(data2)
  yearnum1 = dim1
  yearnum2 = dim2(1)

  
  acr = new((/36,dim2(2)/),float)
  r_test = new((/36,dim2(2)/),float)

 do i=0,35
    if i.le.11 then 
      acr(i,:) = escorc_n(data(0:yearnum1-1),data2(i,0:yearnum2-3,:),0,0)
      r_test(i,:) = rtest(acr(i,:),yearnum1-1,0)
    else 
      if i.le.23 then 
         acr(i,:) = escorc_n(data(0:yearnum1-1),data2(i-12,1:yearnum2-2,:),0,0)
         r_test(i,:) = rtest(acr(i,:),yearnum1-1,0)
      else
         acr(i,:) = escorc_n(data(0:yearnum1-1),data2(i-24,2:yearnum2-1,:),0,0)
         r_test(i,:) = rtest(acr(i,:),yearnum1-1,0)
      end if
    end if
  end do 
  acr!0 = "time"
  acr!1 = "lon"
  acr&lon = data2&lon
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

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(12:dimsizes(u_850&time)-13,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(12:dimsizes(u_850&time)-13,{5:15},{90:130}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;;1901-2011

   ;; surface wind

    u_zonal = dim_avg_n_Wrap(u_1000(:,{-5:5},:),1)  ;;1900-2012
;;;

;;;strong and weak monsoon year 
   total_year = ispan(startyear,endyear-1,1)
   std = dim_stddev_Wrap(DMI_ctrl_JJA)

   inds = ind(DMI_ctrl_JJA.ge.std)
   indw = ind(DMI_ctrl_JJA.le.-1*std)

   indnor = ind(DMI_ctrl_JJA.gt.-0.5*std.and.DMI_ctrl_JJA.lt.0.5*std)

   monstryear = total_year(inds)
   monwekyear = total_year(indw)
;;;;;  
  


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

 corr_sst = corr(DMI_ctrl_JJA,nino34_12(:,:,{59.5:280.5}))

 corr_u1000 = corr(DMI_ctrl_JJA,u_12(:,:,{60:280}))

 printVarSummary(corr_sst)
 printVarSummary(corr_u1000)

 prob_sst = corr_sst
 prob_sst = (/corr_sst@prob/)

 prob_u = corr_u1000
 prob_u = (/corr_u1000@prob/)
 

 ; corr_sst = where(prob_sst.le.0.05, corr_sst, corr_sst@_FillValue)
 ; corr_u1000 = where(prob_u.le.0.05, corr_u1000, corr_u1000@_FillValue)

;;;;

wks = gsn_open_wks("pdf", "SESM-NINO34-corr")

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


    
    ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    
    ; res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    ; res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 
    ; gsn_define_colormap(wks, "precip_diff_12lev")
    ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
    ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

    plot = new(2,graphic)

    res@gsnLeftString = "Corr SEASM&SST"
    plot(0) = gsn_csm_contour(wks,corr_sst,res)
    res@gsnLeftString = "Corr SEASM&surface wind"
    plot(1) = gsn_csm_contour(wks,corr_u1000,res)
;;;

;;;打点

    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.05/)
    res_pro@lbLabelBarOn = False

    res_pro@cnLinesOn=True
    res_pro@cnMonoLineColor =False
    res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
  
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(2,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_sst,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob_u ,res_pro)
   do i=0,1
   overlay( plot(i),iplot(i))
end do

   pres = True

   gsn_panel(wks, plot, (/2,2/), pres)


end 


   

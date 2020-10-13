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
    diri="/mnt/f/share/era/data/data_series/"
    fils1=systemfunc("ls "+diri+"era_pre_global.nc")
    ; fils1=systemfunc("ls "+diri+"w_1.5/era_interim_moda_2000.nc")
    f1 = addfile(fils1, "r")
    time = f1->time
    date = cd_calendar(time, 1)
    ; fils = systemfunc("ls /mnt/f/share/era/data/data_series/w_1.5/era_interim_moda_"+"*"+".nc")
    ; f1 = addfiles(fils, "r")
    ; time=f1[:]->time
    ; date=cd_calendar(f1[:]->time, 1)
    tind=ind(date.le.201612.and.date.ge.197901)
    varname = getfilevarnames(f1)
    timeind=ind(date.le.201612.and.date.ge.197901)
    w_r = f1->w(timeind,{1000:100},{10:-10},{90:290}) 
    
    if isatt(w_r, "add_offset") then 
     w = w_r@scale_factor*w_r + w_r@add_offset
     copy_VarCoords(w_r, w)
  else 
     w = w_r
  end if

  if any(getVarDimNames(w).eq."longitude") then 
     w!2 = "lat"
     w!3 = "lon"
  end if
  printVarSummary(w)


   
;;
  w_avg = dim_avg_n_Wrap(w, 2)
  w_clm = dim_avg_n_Wrap(w_avg, 0)
  
  printMinMax(w_clm, 1)
  exit()
 
  w_clm_pi = lonPivot(w_clm, 30)
  w_clm_pi = (/-1*w_clm_pi/)

;;;画图 
  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/垂直速度的plane")

;;  ************************************************
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="vectical vectocal "
   
    res@gsnMaximize=False

    res@cnFillPalette="hotcold_18lev"
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnLineThicknessF=2.0
    res@cnMinLevelValF = -0.06
    res@cnMaxLevelValF = 0.06
    res@cnLevelSpacingF = 0.01
    res@cnFillOn = True 

    res@cnLineLabelsOn=False
    res@cnLineLabelInterval=2
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False


    res@lbLabelBarOn=True
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"

    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)
    
    res@tmXBLabelFontHeightF=0.05
    res@tmYLLabelFontHeightF=0.02
    res@tmYLMajorThicknessF=2.0
    res@tmXBMajorThicknessF=2.0

    res@tmBorderThicknessF=4.0

    res@tiYAxisString=""
    res@vpHeightF=0.3
    res@vpWidthF=0.7
  ;;; *************************************************
  plot=new(1,graphic)

  plot(0)= gsn_csm_pres_hgt(wks,w_clm_pi,res)

 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 30, 25)+360.0
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon6=fspan(30, 39, 25)


      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel6=conform_dims(dimsizes(ploylon6), ploylevel, 0)
    

      str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
      str2 = unique_string("dum2")
      str3 = unique_string("dum3")
      str31 = unique_string("dum31")
      ployres@gsLineDashPattern=0
       
      plot@$str1$= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)

      ployres@gsLineDashPattern=16

      plot@$str2$= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)

      ployres@gsLineDashPattern=0
      plot@$str3$= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)
      ployres@gsLineDashPattern=0
      plot@$str31$= gsn_add_polyline(wks,plot,  ploylon6,ploylevel6,  ployres)
   draw(plot)
   frame(wks)

end 


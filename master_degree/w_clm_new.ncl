;;***************************************
;;w的气候态
;*******************************************
begin

    fils = systemfunc("ls /mnt/f/share/era/data/data_series/w_1.5/era_interim_moda_"+"*"+".nc")
    
    f1 = addfiles(fils, "r")

    time=f1[:]->time
    date=cd_calendar(time, 1)
   
    tind=ind(date.le.201612.and.date.ge.197901)
    ; w = short2flt(f3[:]->w(tind,{1000:100},{10:-10},{90:290}))
    w_r = f1[:]->w(tind,:,:,:)
    printVarSummary(w_r)

    w = new((/456,dimsizes(w_r&level),dimsizes(w_r&latitude),dimsizes(w_r&longitude)/), float)
   
    ;;; 不同的decade文件有不用的scale_factor
    f2 = addfile(fils(0), "r")
    w(0:11,:,:,:) = short2flt(f2->w)
    
    do i = 0,2
        f3 = addfile(fils(i+1), "r")
        w(12+120*i:12+120*(i+1)-1,:,:,:) = short2flt(f3->w)    ;;;1980-2010
    end do 
    f3 = addfile(fils(4), "r")
    w(12+120*(2+1):,:,:,:) = short2flt(f3->w) 
    
    copy_VarMeta(w_r, w)
    printVarSummary(w)

    ;;
    w_avg = dim_avg_n_Wrap(w(:,{1000:100},{10:-10},{90:290}), 2)
    w_clm = dim_avg_n_Wrap(w_avg, 0)

    printMinMax(w_clm, 1)

    w_clm = (/-1*w_clm/)
    
    w_sea = w_clm({500},{120:270})
    location_max = w_sea&longitude(maxind(w_sea))


;;----------------------------------------------
;;draw gsn_contour plot
;;-----------------------------------------------
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/垂直速度的气候态")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="vectical vectocal "
   
    res@gsnMaximize=False
    res@lbOrientation = "Vertical"
    res@pmLabelBarWidthF = 0.05
    res@lbLabelFontHeightF = 0.01
    res@gsnRightString = ""

  

    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnLineThicknessF=2.0
    res@cnMinLevelValF = -0.02
    res@cnMaxLevelValF = 0.06
    res@cnLevelSpacingF = 0.01
 

    res@cnFillOn = True 
    res@cnFillPalette="BlueYellowRed"
    res@cnFillColors = array_append_record((/49,79,110/), toint(fspan(129,250,7)), 0) 
    res@cnLineLabelsOn=False
    res@cnLineLabelInterval=2
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False
   

    res@lbLabelBarOn=True
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"

    res@tmXBValues=ispan(90,270,30)
    res@tmXBLabels=(/"90~S~o~N~E","120~S~o~N~E","150~S~o~N~E","180~S~o","150~S~o~N~W"\
                    ,"120~S~o~N~W","90~S~o~N~W"/)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(90,270,10)
    
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

    plot(0)= gsn_csm_pres_hgt(wks,w_clm,res)

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
      

      pres = True
      pres@gsMarkerColor = "green"
      pres@gsMarkerIndex = 12
      pres@gsMarkerSizeF = 10
      dum = gsn_add_polymarker(wks, plot, location_max, 500, pres)



   draw(plot)
   frame(wks)



end 

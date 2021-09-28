;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
    fils1 = "/home/yangsong3/data-observation/linshh/data/rea_msf/msf_time_era5_global_2020.nc"  

    f1=addfile(fils1, "r")

    date=cd_calendar(f1->time, 1)

    timestart = 197901
    timeend = 202012
    timeind=ind(date.ge.timestart.and.date.le.timeend)

    msf_time_ctrl=lonPivot(f1->msf_time(timeind,{1000:100},:),30)
    
    msf_time_ctrl = (/msf_time_ctrl/10^9/)


    msf_annual = month_to_annual(msf_time_ctrl, 1)

    
    msf_clm_ctrl = dim_avg_n_Wrap(msf_annual, 0)

    


    size = dimsizes(msf_annual)
    yearnum = size(0)

    msf_rc = msf_clm_ctrl
    msf_rc = (/regCoef_n(ispan(1,yearnum,1), msf_annual , 0, 0)/)
    ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:30},{290:360}), area({0:30}), 1.0, 0) , 0, 0)

    msf_rc  = msf_rc*yearnum



     

 
    ;;;;画图
     wks = gsn_open_wks("pdf", "F_2000-msf-change-ERA5-"+timestart + "-"+ timeend)
      res2 = True
      res2@gsnDraw=False
      res2@gsnFrame=False
      res2@cnFillPalette = "BlueWhiteOrangeRed"
      res2@gsnMaximize=False

      res2@cnLinesOn=True
      res2@cnInfoLabelOn=False
      res2@lbLabelBarOn=True
      res2@cnLineLabelsOn=False
      
      res2@gsnLeftStringFontHeightF=18

      ; res2@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
      res2@tmYRMode="Automatic"
      res2@tmXTOn=False
      res2@tmYROn=False

     res2@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

     res2@cnMinLevelValF       =   -1.0
     res2@cnMaxLevelValF       =   1.0
     res2@cnLevelSindingF      =  0.2

      
       
      res2@tiYAxisString=""
   
      res2@cnFillOn=True
      res2@cnLineThicknessF =4.0
      
      res2@tmXBMode="Explicit"
      res2@tmXBValues=(/30,90,150,210,270,330,387.5/)
      res2@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
              ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
      res2@tmXBMinorOn=True
      res2@tmXBMinorValues=ispan(30,390,10)

      res2@vpHeightF=0.3                      
      res2@vpWidthF=0.7
      res2@cnLinesOn=False

  
       plot=gsn_csm_pres_hgt(wks, msf_rc, res2)
      
        res2@gsnRightString=""
        res2@gsnLeftString=""
        res2@cnLevelSelectionMode = "ExplicitLevels"
        res2@cnLinesOn=True
        res2@cnLevels = ispan(-6,6,1)
        res2@cnLineLabelsOn=True
        res2@cnLineLabelBackgroundColor="Transparent"
        res2@cnLineThicknessF = 2.
        res2@cnFillOn = False
        res2@gsnContourNegLineDashPattern=1
        res2@gsnContourZeroLineThicknessF=4.0
        

        plot_clm = gsn_csm_pres_hgt(wks, msf_clm_ctrl, res2)

      overlay(plot, plot_clm)
  
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
      ploylon4=fspan(30, 39, 25)

      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel4=conform_dims(dimsizes(ploylon4), ploylevel, 0)
      
      str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
      str2 = unique_string("dum2")
      str3 = unique_string("dum3")
      str4 = unique_string("dum4")
      ployres@gsLineDashPattern=0

      plot@$str1$= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)

      ployres@gsLineDashPattern=16

      plot@$str2$= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)

      ployres@gsLineDashPattern=0
      plot@$str3$= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)

      ployres@gsLineDashPattern=0
      plot@$str4$= gsn_add_polyline(wks,plot,  ploylon4,ploylevel4,  ployres)

      text=(/"IWC","PWC","AWC"/)
      x=(/90,210,330/)
      y=(/400,400,400/)
      text1=new(3,graphic)
      textres=True
      textres@txFontThicknessF=0.1
      textres@txFontHeightF=0.02
      do i=0,2
      text1(i)=gsn_add_text(wks, plot, text(i), x(i), y(i), textres)
      end do

  draw(plot)
  frame(wks)
  end 


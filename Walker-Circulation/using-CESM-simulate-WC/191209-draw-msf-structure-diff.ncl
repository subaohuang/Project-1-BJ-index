;;;************************************************************
;; 查看F-pday控制试验以及去掉印度洋海温增温和大西洋增温,太平洋增温的敏感性试验的
;; MSF空间结构的变化。
;; linshuheng 2019.12.06
;;;*********************************************************

;;;;funtion of read data 
function read_msf(fils1)
begin
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.ge.197901.and.date.le.201212)
   
    msf_time_20cr_r=f1->msf_time(timeind,{1000:100},:)


    msf_time_20cr_pi=lonPivot(msf_time_20cr_r, 30)
    msf_time_20cr=msf_time_20cr_pi

    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    printVarSummary(msf_time_20cr)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
      msf_year=month_to_annual(msf_time_20cr, 1)
      msf_year&year=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
      msf_annual_1=dim_avg_n_Wrap(msf_year({1979:1994},:,:), 0)
      msf_annual_2=dim_avg_n_Wrap(msf_year({1995:2012},:,:),0)
    ;;;

    ;;;clmatology 
     msf_clm = dim_avg_n_Wrap(msf_year, 0)
    ;;;
     msf_trend = regCoef_n(msf_year&year, msf_year, 0, 0)
     
    copy_VarCoords(msf_annual_1, msf_clm)
    copy_VarCoords(msf_annual_1, msf_trend)
    ;;;
  return [/msf_clm,msf_trend/]
end 
;;;function of draw plot 
function draw_msf(msf_clm,msf_trend, wks,LeftString)
begin
     numyear = msf_trend
     numyear = 34
     
     msf_add = msf_clm + msf_trend*numyear
     copy_VarCoords(msf_trend, msf_add)
     printVarSummary(msf_add)
     ;;;;
     msf_annual_1 = msf_clm
     msf_annual_2 = msf_add
     ;;;

    ;;;;画图
      res2 = True
      res2@gsnDraw=False
      res2@gsnFrame=False

      res2@gsnMaximize=False

      res2@cnLinesOn=True
      res2@cnInfoLabelOn=False
      res2@lbLabelBarOn=False
      res2@cnLineLabelsOn=False
      res2@gsnLeftString= LeftString
      res2@gsnLeftStringFontHeightF=18

      ; res2@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
      res2@tmYRMode="Automatic"
      res2@tmXTOn=False
      res2@tmYROn=False
      
       
      res2@tiYAxisString=""
      res2@cnLevelSelectionMode ="ExplicitLevels"
      res2@cnLevels=(/0.5/)

      res2@cnFillOn=False
      res2@cnLineColor="blue"
      res2@cnMonoLineColor=True
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

  
     plot=gsn_csm_pres_hgt(wks, msf_annual_1(:,:), res2)
     res2@cnLinesOn=True
     delete(res2@gsnLeftString)
     res2@cnLevels=(/-1/)
     plot_blue1=gsn_csm_pres_hgt(wks, msf_annual_1(:,{0:180}), res2)
     res2@cnLevels=(/0/)

     plot_blue2=gsn_csm_pres_hgt(wks, msf_annual_1(:,{120:280}), res2)
     ; plot_blue2=gsn_csm_pres2_hgt_vector(wks, msf_annual_1(:,{120:270}),ud_djf_1(:,{160:161}),w_djf_1(:,{160:161}), res2)
     delete(res2@cnLevels)
     res2@cnLevels=(/0.5/)
     plot_blue3=gsn_csm_pres_hgt(wks, msf_annual_1(:,{280:}), res2)
     delete(res2@cnLevels)
     res2@cnLevels=(/0/)
     plot_blue4=gsn_csm_pres_hgt(wks, msf_annual_1(:,{180:}), res2)

     overlay(plot, plot_blue1)
     overlay(plot, plot_blue2)
     overlay(plot, plot_blue3)
     overlay(plot, plot_blue4)

    res2@cnLineColor="red"
    res2@cnLinesOn=False
    ; plot1=gsn_csm_pres_hgt(wks, msf_annual_2(:,:), res)
     res2@cnLinesOn=True
     delete(res2@cnLevels)
     res2@cnLevels=(/-1/)
     plot_red1=gsn_csm_pres_hgt(wks, msf_annual_2(:,{0:180}), res2)
     delete(res2@cnLevels)
     res2@cnLevels=(/0/)
   
     plot_red2=gsn_csm_pres_hgt(wks, msf_annual_2(:,{120:280}), res2)
     ; plot_red2=gsn_csm_pres_hgt_vector(wks, msf_annual_2(:,{120:270}),ud_djf_2(:,{140:170}),w_djf_2(:,{140:170}), res)
     delete(res2@cnLevels)
     res2@cnLevels=(/0.5/)
     plot_red3=gsn_csm_pres_hgt(wks, msf_annual_2(:,{280:}), res2)
     delete(res2@cnLevels)
     res2@cnLevels=(/0/)
     plot_red4=gsn_csm_pres_hgt(wks, msf_annual_2(:,{180:}), res2)
     
     overlay(plot, plot_red1)
     overlay(plot, plot_red2)
     overlay(plot, plot_red3)
     overlay(plot, plot_red4)
     
  
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

      ;***********************************************
      ; legend resources
      ;***********************************************
      lgres                    = True
      lgres@lgLineColors     = (/"blue","red"/)
      lgres@lgLineThicknessF   = res2@cnLineThicknessF  ; legend line thickness = contour line thickness
      lgres@lgMonoDashIndex=True
      lgres@lgDashIndex=0
      lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
      lgres@vpWidthF           = 0.15           ; width of legend (NDC)
      lgres@vpHeightF          = 0.1            ; height of legend (NDC)
      lgres@lgPerimColor       = "orange"       ; draw the box perimeter in orange
      lgres@lgPerimOn=False
      lgres@lgPerimThicknessF  = 5.0            ; thicken the box perimeter

      lbid = gsn_create_legend(wks,2,(/" 1979-1994","1995-2012"/),lgres)         ; create legend

      amres = True
      amres@amParallelPosF   = -0.38         ; move legend to the right
      amres@amOrthogonalPosF = -0.38             ; move the legend down
      annoid1 = gsn_add_annotation(plot,lbid,amres)   ; add legend to
return plot
end 

begin 
 ;;------------------read data---------------------------------
	

  fils1  =  "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY_1979-201212.nc"

  fils2 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmindian_197901-201212.nc"

  fils3 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmatlan_197912-201212.nc"

  fils4 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmpacific_197912-201212.nc"
  ; fils4 = "/home/ys17-19/lsh/data/wc-result/msf_F_1850-PDAY-rmallocean_197912-201212.nc"

   
    data_ct = read_msf(fils1)
    data_id = read_msf(fils2)
    data_pa = read_msf(fils3)
    data_al = read_msf(fils4)
    
    msf_trend_ct = data_ct[1]

    msf_trend_id = (/data_ct[1] - data_id[1]/)
    msf_trend_pa = (/data_ct[1] - data_pa[1]/)
    msf_trend_al = (/data_ct[1] - data_al[1]/)

    copy_VarCoords(msf_trend_ct, msf_trend_id)
    copy_VarCoords(msf_trend_ct, msf_trend_pa)
    copy_VarCoords(msf_trend_ct, msf_trend_al) 

    wks = gsn_open_wks("pdf", "cesm-msf-structure-diff")
    plot = new(4,graphic)
    ;;;;get trend pattern
    plot(0) = draw_msf(data_ct[0],data_ct[1],wks,"CT")
    plot(1) = draw_msf(data_ct[0],msf_trend_id,wks,"rmindian")
    plot(2) = draw_msf(data_ct[0],msf_trend_pa, wks,"rmatlan")
    plot(3) = draw_msf(data_ct[0],msf_trend_al,wks,"rmpacific")
    ;;;;

    pres = True
    gsn_panel(wks, plot, (/2,2/), pres)


  end 


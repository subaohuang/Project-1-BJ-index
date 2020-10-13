load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/arrow.ncl"

;************************************************
; 画20cr,ncep,era的msf气候态
;************************************************
begin
    diri="/mnt/f/share/rea_msf/"
    ; diri="/mnt/f/share/CMIP5/msf_AMIP/"
    LeftString="mre"
    ; LeftString="amip"
     ; fils1=systemfunc("ls "+diri+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")
    fils1=systemfunc("ls "+diri+"msf_time_"+LeftString+"_global.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind1=ind(f1->level.eq.400)
    levind2=ind(f1->level.eq.600)
    msf_time_20cr_r=f1->msf_time(:,:,timeind)

    ud_time_20cr_r=f1->ud_time(:,:,timeind)
    w_time_20cr_r=f1->w_time(:,:,timeind)

    msf_time_20cr_pi=lonPivot(msf_time_20cr_r(time|:,level|:,lon|:), 30)
    msf_time_20cr=msf_time_20cr_pi(level|:,lon|:,time|:)

    ud_time_20cr_pi=lonPivot(ud_time_20cr_r(time|:,level|:,lon|:), 30)
    ud_time_20cr=ud_time_20cr_pi(level|:,lon|:,time|:)
    
    w_time_20cr_pi=lonPivot(w_time_20cr_r(time|:,level|:,lon|:), 30)
    w_time_20cr=w_time_20cr_pi(level|:,lon|:,time|:)
 
    msf_time_20cr=(/msf_time_20cr/10^9/)
    
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    
     
    delete(ud_time_20cr&time)
    ud_time_20cr&time=date(timeind)

    delete(w_time_20cr&time)
    w_time_20cr&time=date(timeind)


    printVarSummary(ud_time_20cr)
    printVarSummary(msf_time_20cr)
   
    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)


;;; 
      ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
      msf_son=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "SON")
      msf_jja=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "JJA")
      msf_djf=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "DJF")
      msf_year=month_to_annual(msf_time_20cr(time|:,level|:,lon|:), 1)
      
      msf_son&time=ispan(1979,2012,1)
      msf_year&year=ispan(1979,2012,1)
      msf_jja&time=ispan(1979,2012,1)
      msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
      ;;注意点
      msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
      msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr(:,:,{201212})/)
      ;;;
      printVarSummary(msf_djf)
      msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
      msf_djf_1=dim_avg_n_Wrap(msf_djf({1979:1994},:,:), 0)
      msf_djf_2=dim_avg_n_Wrap(msf_djf({1995:2012},:,:),0)

      msf_jja_1=dim_avg_n_Wrap(msf_jja({1979:1994},:,:), 0)
      msf_jja_2=dim_avg_n_Wrap(msf_jja({1995:2012},:,:),0)

      msf_year_1=dim_avg_n_Wrap(msf_year({1979:1994},:,:), 0)
      msf_year_2=dim_avg_n_Wrap(msf_year({1995:2012},:,:),0)

      ud_djf=month_to_season(ud_time_20cr(time|:,level|:,lon|:), "DJF")
      ud_djf&time=ispan(1979,2012,1)  ;
      ;;注意点
      ud_djf(0:dimsizes(ud_djf&time)-2,:,:)=(/ud_djf(1:,:,:)/)
      ud_djf(dimsizes(ud_djf&time)-1,:,:)=(/ud_time_20cr(:,:,{201212})/)
     
      w_djf=month_to_season(w_time_20cr(time|:,level|:,lon|:), "DJF")
      w_djf&time=ispan(1979,2012,1)  ;
      ;;注意点
      w_djf(0:dimsizes(w_djf&time)-2,:,:)=(/w_djf(1:,:,:)/)
      w_djf(dimsizes(w_djf&time)-1,:,:)=(/w_time_20cr(:,:,{201212})/)



    ud_clm_20cr=dim_avg_n_Wrap(month_to_season(ud_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
    w_clm_20cr=dim_avg_n_Wrap(month_to_season(w_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
    
     ud_djf_1=dim_avg_n_Wrap(ud_djf({1979:1994},:,:), 0)
      ud_djf_2=dim_avg_n_Wrap(ud_djf({1995:2012},:,:),0)

      w_djf_1=dim_avg_n_Wrap(w_djf({1979:1994},:,:), 0)
      w_djf_2=dim_avg_n_Wrap(w_djf({1995:2012},:,:),0)

    
;************************************************
; 画图气候态
;************************************************
 wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料DJF平均_全球2个时间段的比较_0")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    
   
    res@gsnMaximize=False

    res@cnLinesOn=True
    res@cnInfoLabelOn=False
    res@lbLabelBarOn=False
    res@cnLineLabelsOn=False
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
 

    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)
     
    res@tiYAxisString=""
    res@cnLevelSelectionMode ="ExplicitLevels"
    res@cnLevels=(/0.5/)

    res@cnFillOn=False
    res@cnLineColor="blue"
    res@cnMonoLineColor=True
    res@cnLineThicknessF =4.0
   
   
    res@vpHeightF=0.3                      
    res@vpWidthF=0.7
    res@cnLinesOn=False

  
     plot=gsn_csm_pres_hgt(wks, msf_djf_1(:,:), res)
     res@cnLinesOn=True
     
     res@cnLevels=(/-1/)
     plot_blue1=gsn_csm_pres_hgt(wks, msf_djf_1(:,{0:180}), res)
     res@cnLevels=(/0/)

     plot_blue2=gsn_csm_pres_hgt(wks, msf_djf_1(:,{120:270}), res)
     ; plot_blue2=gsn_csm_pres_hgt_vector(wks, msf_djf_1(:,{120:270}),ud_djf_1(:,{160:161}),w_djf_1(:,{160:161}), res)
     delete(res@cnLevels)
     res@cnLevels=(/0.5/)
     plot_blue3=gsn_csm_pres_hgt(wks, msf_djf_1(:,{270:}), res)
     delete(res@cnLevels)
     res@cnLevels=(/-1/)
     plot_blue4=gsn_csm_pres_hgt(wks, msf_djf_1(:,{180:}), res)

     overlay(plot, plot_blue1)
     overlay(plot, plot_blue2)
     overlay(plot, plot_blue3)
     overlay(plot, plot_blue4)

    res@cnLineColor="red"
    res@cnLinesOn=False
    ; plot1=gsn_csm_pres_hgt(wks, msf_djf_2(:,:), res)
     res@cnLinesOn=True
     delete(res@cnLevels)
     res@cnLevels=(/-1/)
     plot_red1=gsn_csm_pres_hgt(wks, msf_djf_2(:,{0:180}), res)
     delete(res@cnLevels)
     res@cnLevels=(/0/)
   
     plot_red2=gsn_csm_pres_hgt(wks, msf_djf_2(:,{120:270}), res)
     ; plot_red2=gsn_csm_pres_hgt_vector(wks, msf_djf_2(:,{120:270}),ud_djf_2(:,{140:170}),w_djf_2(:,{140:170}), res)
     delete(res@cnLevels)
     res@cnLevels=(/0.5/)
     plot_red3=gsn_csm_pres_hgt(wks, msf_djf_2(:,{270:}), res)
     delete(res@cnLevels)
     res@cnLevels=(/-1/)
     plot_red4=gsn_csm_pres_hgt(wks, msf_djf_2(:,{180:}), res)
     
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
  lgres@lgLineThicknessF   = res@cnLineThicknessF  ; legend line thickness = contour line thickness
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
  annoid1 = gsn_add_annotation(plot,lbid,amres)   ; add legend to plot

    ; ; overlay(plot, plot1)
    ;     xarrow1 = (/162,162/)
    ;   yarrow1 = (/  700,  200/)
    ;   xarrow2 = (/ 50,150/)
    ;   yarrow2 = (/  0,  0/)
      
    ;       ;--Customize each arrow and draw again.
    ;       resp=True
    ;   resp@gsLineColor       = "blue"
    ;   resp@gsLineDashPattern = 0    ; solid line
    ;   resp@gsLineThicknessF  = 2.0  ; 3x thicker
    ;   ; resp@gs
    ;   plot=arrow(wks,plot,xarrow1,yarrow1,resp)
    draw(plot)
    frame(wks)
    
 end
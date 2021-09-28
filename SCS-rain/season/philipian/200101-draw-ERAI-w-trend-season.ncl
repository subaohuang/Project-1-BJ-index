;;********************************************************
;;回归ERA的u,v,slp,z到水汽收支方程的不同的变化上
;;
;;lnshuheng 2019.12.24
;**********************************************************
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

;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"wwnd"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
    ; ;;;region 
    ;   lat1 = 0
    ;   lat2 = 15
    ;   lon1 = 110
    ;   lon2 = 155
    ; ;;;
    
    ;;;region 
      lat1 = -5
      lat2 = 45
      lon1 = 90
      lon2 = 170
    ;;;


    ;;;time span
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    w = short2flt(f1->w(timeind,:,{lat1:lat2},{lon1:lon2}))
    w = (/w*-1/)


    

    level = w&level


   ;;;处理成年异常数据
     MONTH = (/"MAM","JJA","SON","DJF"/)
   ;;;
  ;;;
  


  ;;;;
    wks2   = gsn_open_wks("pdf", "ERA-In中的W趋势平面图-季节-mk")
    plot=new(4,graphic)  
    str = (/"1","2","3","4"/)
    iplot = new(4,graphic)
   
    se = (/"(a)","(b)","(c)","(d)"/)

  do i = 0,3
    w_annual = month_to_season(w, MONTH(i))

    if MONTH(i).eq."DJF" then 
       w_annual(0:dim-2,:,:,:) = (/w_annual(1:dim-1,:,:,:)/)
       w_annual(dim-1,:,:,:) = (/w(dim*12-1,:,:,:)/)
    end if
    
    ; ;;;垂直积分
    ; dp=dpres_plevel(level,1013,min(level),0)              ; caculate dp
    ; dp_cfm=conform_dims(dimsizes(w_annual),dp,1)          ; confirm dp to caculate sum(var*dp)
    ; wgt_w = w_annual(:,1,:,:)
    ; wgt_w=(/dim_sum_n(dp_cfm*w_annual,1)/)               ; sum(var*dp)

    wgt_w = w_annual(:,{500},:,:)

    ; wgt_t = t_annual(:,1,:,:)
    ; wgt_t=(/dim_sum_n(dp_cfm*t_annual,1)/)               ; sum(var*dp)

    ; ;;;



   ;;;计算趋势

      
    wgt_w_rc = regCoef_n(ispan(1,yearnum,1), wgt_w, 0, 0)
    copy_VarCoords(wgt_w(1,:,:), wgt_w_rc)

    prob_w = rc_test(wgt_w_rc)
    

    ; mk = trend_manken(wgt_w, False, 0)
    ; wgt_w_rc := mk(1,:,:)
    ; copy_VarCoords(wgt_w(1,:,:), wgt_w_rc)

    ; prob_w := mk(0,:,:)
    ; prob_w = (/1 - prob_w/)
    ; copy_VarCoords(wgt_w(1,:,:), prob_w)

 


    ; prob_t  = prob_t*10^6
   ;;;
   wgt_w_rc = (/wgt_w_rc*yearnum/) *100  
    


;;;=================================================
;;;平面分布图
;;===================================================
   
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic = False
   
    ; res@cnFillPalette="BlueWhiteOrangeRed"
    ; res@cnFillPalette="BlueRed"
    ; res@cnFillPalette = "MPL_PuOr"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
    res@lbLabelStride=2
    res@lbOrientation="Vertical"
    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -5
    ; res@mpMaxLatF             = 45
    ; res@mpMinLonF             =100
    ; res@mpMaxLonF             = 170 

    res@mpMinLatF             = -5
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 100
    res@mpMaxLonF             = 170 
    ; res@mpFillOn = True
    ; res@mpFillColor = "red"

    ; res@mpFillDrawOrder = "PreDraw"


    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"


    res@tmXBLabelFontHeightF=0.025
    res@tmYLLabelFontHeightF=0.025
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0
    res@tmXBValues = (/100,120,140,160/) 
    res@tmXBLabels = (/100,120,140,160/) +""+"E" 


    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -0.02
    ; res@cnMaxLevelValF       =   0.02
    ; res@cnLevelSpacingF      =  0.005

    res@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
    gsn_define_colormap(wks2, "MPL_PuOr")
    
    res@cnLevels = (/-0.02,-0.015,-0.01,-0.005,0.005,0.01,0.015,0.02/)*100
    res@cnFillColors = (/22,32,45,60,-1,80,97,100,110/) 


    res@gsnLeftString = se(i) + " "+MONTH(i) + " W at 500 hPa "
    res@gsnLeftStringFontHeightF = 0.03
    plot(i) = gsn_csm_contour_map(wks2, wgt_w_rc, res)


  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=1

 
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
 
   plot@$str(i)$ = gsn_add_polyline(wks2, plot(i), pgx_pool, pgy_pool, pgres)
 

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0.0,0.1/)
    res_pro@lbLabelBarOn = False

    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.001
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True


    iplot(i)=gsn_csm_contour(wks2,prob_w,res_pro)
    overlay(plot(i), iplot(i))
    delete(wgt_w_rc)
end do 

 pres = True
 pres@gsnPanelLabelBar = True
 ; pres@lbOrientation = "Vertical"
 pres@lbLabelFontHeightF = 0.015

;  pres@lbOrientation = "Vertical"
; pres@pmLabelBarHeightF = 0.3
 gsn_panel(wks2, plot, (/2,2/), pres)

    
end 















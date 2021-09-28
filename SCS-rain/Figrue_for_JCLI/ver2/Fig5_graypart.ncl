;;********************************************************
;;画散点图,得到关系
;;
;;lsh 2020-11-04
;**********************************************************



begin
;;;read data
lat1 = -10
lat2 =  30
lon1 =  100
lon2 = 180

; startmon = 198201
; endmon = 201912

; startmon = 199501
; endmon = 201512

 startmon = 197901
 endmon = 201512

;;;读取数据

    f1 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
    date:=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    olr = short2flt(f1->olr(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

    printVarSummary(olr)

    olr = (/dtrend_msg_n(ispan(1,dimsizes(olr&time),1),olr,False,True,0)/)

    ; f2 = addfile("/home/yangsong3/data-observation/SST/NOAA_OISST/sst.mnmean.nc", "r")
    ; date:=cd_calendar(f2->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; sst = short2flt(f2->sst(timeind,{lat1:lat2+1.0},{lon1-1.0:lon2+1.0})) ;;;2.5 ×2.5 grid
    ; printVarSummary(sst)
    ;;;;


    f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
    date:=cd_calendar(f2->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_r =  lonFlip(f2->sst(timeind,:,:))

    sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

    sst = (/dtrend_msg_n(ispan(1,dimsizes(sst&time),1),sst,False,True,0)/)


    printVarSummary(sst)
   
   sst_tropical = avg(dim_avg_n_Wrap(sst_r(:,{-30:30},{0:360}),(/1,2/)))
    print("热带平均的SST= "+sst_tropical)


    fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
    f3 = addfiles(fils3, "r")

    date:=cd_calendar(f3[:]->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude
    lon = f3[1]->longitude

    indlat = ind(lat.ge.-10.0.and.lat.le.30.0)
    indlon = ind(lon.ge.100.0.and.lon.le.180.0)


    div_r = new((/dimsizes(time),7, dimsizes(indlat),dimsizes(indlon)/),float)
      

    do i = 0,dimsizes(timeind)/12-1

      div_r(i*12:(i+1)*12-1,:,:,:) = short2flt(f3[i]->d(:,{1000:850},indlat(::-1),indlon))
      
    end do
    ;

    div = dim_avg_n_Wrap(div_r,1)


    div = (/dtrend_msg_n(ispan(1,dimsizes(div&time),1),div,False,True,0)/)



    div = (/div*10^5/)
    div = (/-1*div/)           ;;;转化成con
;;;;



;;;;chazhi 
print("**********start interp************")
    ; xi = sst&lon
    ; yi = sst&lat

    xi = sst&longitude
    yi = sst&latitude

    xo = olr&lon
    yo = olr&lat

  
    ;;;;chazhi SST
    sst_chazhi = linint2_Wrap(xi, yi, sst, True, xo, yo, 0)

    sst_chazhi!1= "lat"
    sst_chazhi!2= "lon"

    sst_chazhi&lat = yo
    sst_chazhi&lon = xo
    ;;;
    


    xi := div&longitude
    yi := div&latitude

    ;;;;chazhi div
    div_chazhi = linint2_Wrap(xi, yi, div, True, xo, yo, 0)

    div_chazhi!1= "lat"
    div_chazhi!2= "lon"

    div_chazhi&lat = yo
    div_chazhi&lon = xo
    ;;;

    olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)
;;;





;;;;Monthly bin

    sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
    olr = (/runave_n_Wrap(olr, 3, 0, 0)/)
    div_chazhi = (/runave_n_Wrap(div_chazhi, 3, 0, 0)/)



    sst_one = ndtooned(sst_chazhi(:,{0:20},{125:170}))
    olr_one = ndtooned(olr(:,{0:20},{125:170}))

    div_one = ndtooned(div_chazhi(:,{0:20},{125:170}))
;;;;;




;;;bin month for OLR and SST
    bin_num =  (30.75-25.0)/0.25   ;;;;bin的个数，也是中间数

    bin_sst  = fspan(25.0-0.125, 30.75-0.125, toint(bin_num) +1)  ;;;bin最左边的数字
    olr_bin_sst = new(toint(bin_num),float)
    printMinMax(sst_one, 1)
    
    std_olr_sst = olr_bin_sst


    do i = 0,bin_num-1
        in := ind(sst_one.ge.bin_sst(i).and.sst_one.lt.(bin_sst(i)+0.25))
        olr_bin_sst(i) = dim_avg(olr_one(in))
        std_olr_sst(i) = dim_stddev(olr_one(in))
    end do 

   dolr_dsst = center_finite_diff(olr_bin_sst, bin_sst(1:) - 0.125, False, 0)


    mnmx = new((/2,dimsizes(olr_bin_sst)/),float)

    mnmx(0,:) = olr_bin_sst - std_olr_sst
    mnmx(1,:) = olr_bin_sst + std_olr_sst


    x = bin_sst(1:) - 0.125

    ind1 = ind(x.eq.26.5)
    ind2 = ind(x.eq.28.75)
    ind3 = ind(x.eq.29.75)

    rc1 = regline(fspan(26.5,28.75,ind2-ind1+1), dolr_dsst(ind1:ind2))
    print(rc1@yintercept)

    rc2 = regline(fspan(28.75,29.75,ind3-ind2+1), dolr_dsst(ind2:ind3))

    print(rc2@yintercept)
    print (rc1+"*x"+"+" +rc1@yintercept)
    print (rc2+"*x"+"+" +rc2@yintercept)
    

    ; y1 = rc1*fspan(26.5,28.75,ind2-ind1+1) + rc1@yintercept
    ; y2 = rc2*fspan(28.75,30.0,toint((30.0-28.75)/0.25)+1) + rc2@yintercept
;;;;;;;



;;;bin month for OLR and DIV  
   ; bin_div = ispan(-30-5, 50-5, 5)*0.01  ;;; 0.1 最左边的
   ; bin_num := dimsizes(bin_div) -1

   bin_num :=  (0.65+0.3)/0.05   ;;;;bin的个数，也是中间数
    
   bin_div  = fspan(-0.3-0.025, 0.65-0.025, toint(bin_num) +1)  ;;;bin最左边的数字

   olr_bin_div = new(toint(bin_num),float)

   std_olr_div = olr_bin_div


 do i = 0,bin_num-1
        in := ind(div_one.ge.bin_div(i).and.div_one.lt.(bin_div(i)+0.05))
        olr_bin_div(i) = dim_avg(olr_one(in))
        std_olr_div(i) = dim_stddev(olr_one(in))
 end do 
olr_bin_div(0) = 265.0
  dolr_ddiv = center_finite_diff(olr_bin_div, bin_div(1:) - 0.025, False, 0)
;;;;


    dolr_ddiv = (/dolr_ddiv/10.0/)


    mnmx1 = new((/2,dimsizes(olr_bin_div)/),float)

    mnmx1(0,:) = olr_bin_div - std_olr_div
    mnmx1(1,:) = olr_bin_div + std_olr_div
   

;;;bin SST and DIV
    bin_num :=  (30.75-25.0)/0.25   ;;;;bin的个数，也是中间数
    
    bin_sst  = fspan(25.0-0.125, 30.75-0.125, toint(bin_num) +1)  ;;;bin最左边的数字
    div_bin_sst = new(toint(bin_num),float)
    printMinMax(sst_one, 1)

    std_div_sst = div_bin_sst

    do i = 0,bin_num-1
        in := ind(sst_one.ge.bin_sst(i).and.sst_one.lt.(bin_sst(i)+0.25))
        div_bin_sst(i) = dim_avg(div_one(in))
        std_div_sst(i) = dim_stddev(div_one(in))
    end do 

   ddiv_dsst = center_finite_diff(div_bin_sst, bin_sst(1:) - 0.125, False, 0)
   ddiv_dsst  = (/ddiv_dsst*10.0/)
;;;;

    l = div_bin_sst(ind((bin_sst(1:)-0.125).eq.29.0))

  
    mnmx2 = new((/2,dimsizes(olr_bin_sst)/),float)

    mnmx2(0,:) = div_bin_sst - std_div_sst
    mnmx2(1,:) = div_bin_sst + std_div_sst


;;;;draw plot 
    wks=gsn_open_wks("pdf", "SST-OLR-DIV-relation-tropical-wnp-monlydata-panel")

    plot = new(6, graphic)


    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False

    res@xyMonoLineColor = True
    res@xyMarkerColor ="gray30"

    res@xyMarkLineModes   = (/"Markers"/)
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.001    ; Marker size (default 0.01)
   
    res@trYReverse =  True
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=290.0
    res@trYMinF=180.0

    
    res@tmYLValues    = fspan(290,180,12)
    res@tmYLMinorValues = fspan(290,180,23)
    res@tmYLLabels    = fspan(290,180,12)
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=31.0
    res@trXMinF=24.0
    res@tmXBValues    = fspan(24.0,31.0,8)
    res@tmXBLabels    = fspan(24.0,31.0,8)


    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

    res@tmXBLabelFontHeightF=0.02
    res@tmXBMajorLengthF=0.008        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0.008  ;调整刻度的长度，注意2个命令要同时使用
    res@tmXUseBottom=False
    res@tmXBLabelDeltaF=-0.8

    res@tmYLLabelFontHeightF=0.02
    res@vpXF=0.15
    res@vpYF=0.9
    res@vpHeightF=0.5
    res@vpWidthF=0.6
    
    res@tiMainString="(a) ~F25~OLR and SST"

    res@tiYAxisString="~F25~OLR [W/m~S~2~N~]"
    res@tiYAxisFontHeightF=0.02
    
    res@tiXAxisString="~F25~SST [~S~o~N~C]"
    res@tiXAxisFontHeightF=0.02
    
    res@xyCurveDrawOrder = "PreDraw"
  

    plot(0) =gsn_csm_xy(wks, sst_one, (/olr_one/), res)
    
    ;;;添加标准差
    res@xyMarkLineMode   = "Lines" 
    res@gsnXYFillColors = "blue"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    res@gsnXYFillOpacities=0.1

    res@gsnXRefLine = (/28.0,29.0/)
    res@gsnYRefLine = 240.0

    plot_std1= gsn_csm_xy(wks, bin_sst(1:) - 0.125,mnmx(0:1,:),res)       ; create line plot

    overlay(plot(0), plot_std1)
  

    res@trXMaxF= 0.8
    res@trXMinF= -0.6
    res@tmXBValues    := ispan(-6,8,2)*0.1
    res@tmXBLabels    := ispan(-6,8,2)*0.1


    res@tiMainString="(b) ~F25~OLR and Conv "
    res@tiYAxisString="~F25~OLR [W/m~S~2~N~]"
    res@tiXAxisString="~F25~Conv [10~S~-5~N~s~S~-1~N~]"

    delete(res@gsnXRefLine)
    delete(res@gsnYRefLine)


    plot(2) =gsn_csm_xy(wks, div_one, (/olr_one/), res)
    res@gsnXRefLine = (/0.0,l/)
    res@gsnYRefLine = 240.0
    
    plot_std2= gsn_csm_xy(wks, bin_div(1:) - 0.025,mnmx1(0:1,:),res)       ; create line plot

    overlay(plot(2), plot_std2)


    res@trXMaxF=31.0
    res@trXMinF=24.0
    res@tmXBValues    := fspan(24.0,31.0,8)
    res@tmXBLabels    := fspan(24.0,31.0,8)
    
    res@trYMaxF=.5
    res@trYMinF=-0.5
    res@tmYLValues    := ispan(-5,5,1)*0.1
    res@tmYLLabels    := ispan(-5,5,1)*0.1
    
    res@tiMainString="(c) ~F25~Conv and SST "
    res@tiXAxisString="~F25~SST [~S~o~N~C]"
    res@tiYAxisString="~F25~Conv [10~S~-5~N~s~S~-1~N~]"
    res@trYReverse =  False
   
   
    
    delete(res@gsnXRefLine)
    delete(res@gsnYRefLine)

    
    plot(4) =gsn_csm_xy(wks, sst_one, (/div_one/), res)
     
    res@gsnYRefLine = 0.0
    res@gsnXRefLine = (/28.0,29.0/)

    plot_std3= gsn_csm_xy(wks, bin_sst(1:) - 0.125,mnmx2(0:1,:),res)       ; create line plot

    overlay(plot(4), plot_std3)
    
    delete( res@gsnXYFillColors)
    delete(res@xyLineColor)    
    delete(res@gsnXRefLine)

      ; delete(res@gsnXRefLine)
    delete(res@gsnYRefLine)



;;;


;;;添加拟合的线

    plres=True
    plres@gsLineThicknessF=3.0
    plres@gsLineColor = "blue"
    plres@gsLineDashPattern=0

    plot_line=gsn_add_polyline(wks, plot(0), bin_sst(1:) - 0.125, olr_bin_sst, plres)
    plot_line1=gsn_add_polyline(wks, plot(2), bin_div(1:) - 0.025, olr_bin_div, plres)
    plot_line2=gsn_add_polyline(wks, plot(4), bin_sst(1:) - 0.125, div_bin_sst, plres)

    mkres = True
    mkres@gsMarkerColor = "blue"
    mkres@gsMarkerIndex = 1
    mkres@gsMarkerSizeF = 0.02

    ; plot_mk=gsn_add_polymarker(wks, plot(0), bin_sst(1:) - 0.125, olr_bin_sst, mkres)
    ; plot_mk1=gsn_add_polymarker(wks, plot(2), bin_div(1:) - 0.025, olr_bin_div, mkres)
    ; plot_mk2=gsn_add_polymarker(wks, plot(4), bin_sst(1:) - 0.125, div_bin_sst, mkres)
;;;


;;;;画偏导数的线

    res@trXMaxF = 31.0
    res@trXMinF = 24.0
    res@tmXBValues    := fspan(24.0,31.0,8)
    res@tmXBLabels    := fspan(24.0,31.0,8)


    res@tmYLValues    := fspan(20,-30,11)
    res@tmYLMinorValues := fspan(20,-30,11)
    res@tmYLLabels    := fspan(20,-30,11)
    res@trYMaxF=20.0
    res@trYMinF=-30.0


    res@xyMarkerColor ="black"
    res@xyMarkLineModes   := (/"Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=0

    res@tiMainString="(d) ~F34~6~F25~(OLR)/~F34~6~F25~(SST) from (a)"
    res@tiYAxisString="~F25~Chage rate [W/m~S~2~N~/~S~o~N~C]"
 
    res@gsnYRefLineThicknessF = 1.0
    
    ; res@gsnYRefLine = 0.0

    ; res@gsnXRefLine = (/28.0/) ;;;这里会报错，只能画一条，很奇怪
    res@trYReverse =  True
    
    m1 =  res@tmXBValues
    m2 = res@tmYLValues

    plot(1) =gsn_csm_xy(wks, bin_sst(1:) -0.125, (/dolr_dsst/), res)
   

    res@tmYLValues    := fspan(8,-20,8)
    res@tmYLMinorValues := fspan(8,-20,15)
    res@tmYLLabels    := fspan(8,-20,8)
    res@trYMaxF=8.0
    res@trYMinF=-20.0

    res@trXMaxF=0.8
    res@trXMinF=-0.6
    res@tmXBValues    := ispan(-6,8,2)*0.1
    res@tmXBLabels    := ispan(-6,8,2)*0.1

    res@tiMainString="(e) ~F34~6~F25~(OLR)/~F34~6~F25~(Conv)  from (b)"
    res@tiYAxisString="~F25~Chage rate [10~S~-6~N~s~S~-1~N~/~S~o~N~C]"
    res@tiXAxisString="~F25~Conv [10~S~-5~N~s~S~-1~N~]"

    ; res@gsnYRefLine = 0.0
    ; res@gsnXRefLine = (/0.0/)

    m3 =  res@tmXBValues
    m4 = res@tmYLValues

    plot(3) =gsn_csm_xy(wks, bin_div(1:) -0.025, (/dolr_ddiv/), res)
    ; plot_addref2 =gsn_add_polyline(wks,plot(3), conform(bin_sst(1:),l,0), res@tmYLValues, Lineres)
    ; plot_addref2 =gsn_add_polyline(wks,plot(3), res@tmXBValues, conform(res@tmYLValues,0.0,0), Lineres)

    
   
    res@tmYLValues    := fspan(-1,2,7)
    ; res@tmYLMinorValues := fspan(-2,2,9)
    res@tmYLLabels    := fspan(-1,2,7)
    res@trYMaxF=2.0
    res@trYMinF=-1.0
    res@trXMaxF=31.0
    res@trXMinF=24.0
    res@tmXBValues    := fspan(24.0,31.0,8)
    res@tmXBLabels    := fspan(24.0,31.0,8)


    res@tiMainString="(f) ~F34~6~F25~(Conv)/~F34~6~F25~(SST)  from (c)"
    res@tiYAxisString="~F25~Chage rate [10~S~-6~N~s~S~-1~N~/~S~o~N~C]"
    res@tiXAxisString="~F25~SST [~S~o~N~C]"

    res@trYReverse =  False
    res@gsnYRefLine = 0.0
    res@gsnXRefLine = 28.0

    res@gsnXRefLine = (/28.0/)

    m5 =  res@tmXBValues
    m6 = res@tmYLValues
    plot(5) =gsn_csm_xy(wks, bin_sst(1:) - 0.125, (/ddiv_dsst/), res)

    ; plot_addref3 =gsn_add_polyline(wks,plot(5), conform(bin_sst(1:),29.0,0), res@tmYLValues, Lineres)
    ; plot_addref3 =gsn_add_polyline(wks,plot(5), res@tmXBValues, conform(res@tmYLValues,0.0,0), Lineres)

;;;添加阴影部分

  ;;;;;; Set up arrays to hold bar coordinates
  xbars1 = (/28,29,29,28,28/)
  ybars1 = (/20,20,-30,-30,20/)

  xbars2 = (/0,0.142,0.142,0,0/)
  ybars2 = (/8,8,-20,-20,8/)
  
  xbars3 = (/28,29,29,28,28/)
  ybars3 = (/-1,-1,2,2,-1/)
  

;;;;;; Create the bars and attach to appropriate plot
  gres             = True
  ; gres@gsFillColor = "chartreuse1"
  gres@gsFillColor = "lightgray"
  gres@tfPolyDrawOrder =  "PreDraw"

  ; gres@gsFillOpacityF = 0.5

 
  poly1 = gsn_add_polygon(wks,plot(1),xbars1,ybars1,gres)
  poly2 = gsn_add_polygon(wks,plot(3),xbars2,ybars2,gres)
  poly3 = gsn_add_polygon(wks,plot(5),xbars3,ybars3,gres)
;;;

;;;画参考线
    Lineres =  True
    Lineres@tfPolyDrawOrder = "PreDraw"

    ; plot_addref1 =gsn_add_polyline(wks,plot(1), conform(bin_sst(1:),29.0,0), res@tmYLValues, Lineres)
    plot_addref1 = gsn_add_polyline(wks,plot(1), m1, conform(m2,0.0,0), Lineres)
    plot_addref2 = gsn_add_polyline(wks,plot(3), m3, conform(m4,0.0,0), Lineres)
    plot_addref3 = gsn_add_polyline(wks,plot(5), m5, conform(m6,0.0,0), Lineres)


 pres  = True
 gsn_panel(wks, (/plot(0),plot(2),plot(4),plot(1),plot(3),plot(5)/), (/2,3/), pres)

 end 




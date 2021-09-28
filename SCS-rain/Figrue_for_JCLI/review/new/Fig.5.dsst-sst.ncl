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


    ; startmon = 197812
    ; endmon = 201511

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_r =  lonFlip(f2->sst(timeind,:,:))

    sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

    sst = (/dtrend_msg_n(ispan(1,dimsizes(sst&time),1),sst,False,True,0)/)


    printVarSummary(sst)
   
   sst_tropical = avg(dim_avg_n_Wrap(sst_r(:,{-30:30},{0:360}),(/1,2/)))
    print("热带平均的SST= "+sst_tropical)


    fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
    f3 = addfiles(fils3, "r")

    startmon = 197901
    endmon = 201512

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


    olr = (/center_finite_diff_n(sst_chazhi, 1, False, 0, 0)/)
;;;





;;;;Monthly bin

    ; sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
    ; olr = (/runave_n_Wrap(olr, 3, 0, 0)/)
    ; div_chazhi = (/runave_n_Wrap(div_chazhi, 3, 0, 0)/)



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


;;;;draw plot 
    wks=gsn_open_wks("pdf", "norm3-dSST-sst-tropical-wnp-monlydata-panel")

    plot = new(1, graphic)


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
    res@trYMaxF=2.4
    res@trYMinF=-2.0

    
    res@tmYLValues    = ispan(-20,24,4)*0.1
    ; res@tmYLMinorValues = ispan(-20,24,4)*0.1
    res@tmYLLabels    = ispan(-20,24,4)*0.1
    
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

    res@gsnXRefLine = (/28.0/)
    res@gsnYRefLine = 0.0

    plot_std1= gsn_csm_xy(wks, bin_sst(1:) - 0.125,mnmx(0:1,:),res)       ; create line plot

    overlay(plot(0), plot_std1)
  



;;;添加拟合的线

    plres=True
    plres@gsLineThicknessF=3.0
    plres@gsLineColor = "blue"
    plres@gsLineDashPattern=0

    plot_line=gsn_add_polyline(wks, plot(0), bin_sst(1:) - 0.125, olr_bin_sst, plres)
   
  draw(plot)
  frame(wks)

 end 




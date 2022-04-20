;;********************************************************
;;画散点图,得到关系
;;
;;lsh 2020-11-04
;**********************************************************
; function chazhi(data)
;  begin
;    xi = data&longitude
;    yi = data&latitude
;    xo = fspan(0.5,359.5,360)
;    yo = fspan(-89.5,89.5,180)
 
;    data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
;    data_chazhi!1= "lat"
;    data_chazhi!2= "lon"
   
;    data_chazhi&lat = yo
;    data_chazhi&lon = xo
;    copy_VarAtts(data, data_chazhi)
 
;    return data_chazhi
; end 
begin
  ;;;read data
  lat1 = -10
  lat2 =  30
  lon1 =  100
  lon2 = 180

  ; startmon = 198201
  ; endmon = 201912

  startmon = 197901
  endmon = 201512

  f1 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
  date:=cd_calendar(f1->time, 1)
  timeind :=ind(date.ge.startmon.and.date.le.endmon)
  olr = short2flt(f1->olr(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

  printVarSummary(olr)



  f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
  date:=cd_calendar(f2->time, 1)

  timeind:=ind(date.le.201512.and.date.ge.197901)

  sst_r =  lonFlip(f2->sst(timeind,:,:))

  sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

  printVarSummary(sst)


  fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
  f3 = addfiles(fils3, "r")

  date:=cd_calendar(f3[:]->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)
  time = f3[:]->time(timeind)
  lat = f3[1]->latitude
  lon = f3[1]->longitude


  div_r = new((/dimsizes(time), dimsizes(lat),dimsizes(lon)/),float)


  do i = 0,dimsizes(timeind)/12-1

  div_r(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->d(:,{200},:,:))

  end do
 
  ; div = dim_avg_n_Wrap(div_r(:,:,{lat1-1:lat2+1},{lon1-1:lon2+1}),1)
  div = div_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})


  div = (/div*10^5/)
  div = (/-1*div/)

  ; div = (/rmMonAnnCycTLL(div)/)


;;;;chazhi 
   print("**********start interp************")
    ; xi = sst&lon
    ; yi = sst&lat

    xi = sst&longitude
    yi = sst&latitude

    xo = olr&lon
    yo = olr&lat



    sst_chazhi = linint2_Wrap(xi, yi, sst, True, xo, yo, 0)

    sst_chazhi!1= "lat"
    sst_chazhi!2= "lon"

    sst_chazhi&lat = yo
    sst_chazhi&lon = xo


   div_chazhi = linint2_Wrap(div&longitude, div&latitude, div, True, xo, yo, 0)

   div_chazhi!1= "lat"
   div_chazhi!2= "lon"

   div_chazhi&lat = yo
   div_chazhi&lon = xo
;;;;
   delete(sst)
   delete(div)
   
;;;;
    printVarSummary(sst_chazhi)
    olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)
;;;;




;;;;Climatology bin
    ; sst_clm = dim_avg_n_Wrap(sst_chazhi(:,{0:20},{125:160}), 0)   
    sst_clm = dim_avg_n_Wrap(month_to_season(sst_chazhi(:,{0:20},{125:160}),"MAM"), 0)   
    olr_clm = dim_avg_n_Wrap(olr(:,{0:20},{125:160}), 0)


    sst_clm_one = ndtooned(sst_clm)
    olr_clm_one = ndtooned(olr_clm)

;;;


;;;;anomaly 

    olr = (/rmMonAnnCycTLL(olr)/)
    sst_chazhi = (/rmMonAnnCycTLL(sst_chazhi)/)
    div_chazhi = (/rmMonAnnCycTLL(div_chazhi)/)
;;;

;;;;Monthly bin

   ; sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
   ; olr = (/runave_n_Wrap(olr, 3, 0, 0)/)

   dim = dimsizes(sst_chazhi)
   
sst_chazhi_june = sst_chazhi(2::12,:,:)
sst_chazhi_july = sst_chazhi(3::12,:,:)
sst_chazhi_aug = sst_chazhi(4::12,:,:)

sst_chazhi_jj =  array_append_record(sst_chazhi_june,sst_chazhi_july,0)
sst_chazhi_jja = array_append_record(sst_chazhi_jj , sst_chazhi_aug,0)


olr_june = olr(2::12,:,:)
olr_july = olr(3::12,:,:)
olr_aug = olr(4::12,:,:)

olr_jj =  array_append_record(olr_june,olr_july,0)
olr_jja = array_append_record(olr_jj , olr_aug,0)


div_june = div_chazhi(2::12,:,:)
div_july = div_chazhi(3::12,:,:)
div_aug = div_chazhi(4::12,:,:)

div_jj =  array_append_record(div_june,div_july,0)
div_jja = array_append_record(div_jj , div_aug,0)

 

  
   ;  sst_one = ndtooned(sst_chazhi_jja(:,{0:15},{125:160}))
   ;  olr_one = ndtooned(olr_jja(:,{0:15},{125:160}))

    ; sst_one = ndtooned(sst_chazhi(:,{0:20},{125:160}))
;;;;


;;;;计算单点回归系数

    corr = escorc_n(sst_chazhi_jja, olr_jja, 0, 0)                    ;;计算相关性
    copy_VarCoords(sst_chazhi(1,:,:), corr)
    printVarSummary(corr)
    std1 = dim_stddev_n(olr_jja, 0)
    std2 = dim_stddev_n(sst_chazhi_jja, 0)
    
    rc = corr
    rc = (std1/std2)*corr

    corr_div_sst = escorc_n(sst_chazhi_jja, div_jja, 0, 0) 
    corr_div_olr = escorc_n(olr_jja, div_jja, 0, 0) 

    rc_partcial   = corr 
    rc_partcial = ((/corr -  corr_div_olr*corr_div_sst/)/(/1 - corr_div_sst^2/))*(std1/std2)
    



;;;;;

    rc_one = ndtooned(rc({0:20},{125:160}))
    ; rc_one = ndtooned(rc_partcial({0:20},{125:160}))
;;;;



;;;bin 
   
   sstmax = 30.0
   sstmin = 27.5
   ; sstmax = 30.0
   ; sstmin = 27.5

    bin_num =  (sstmax- sstmin)/0.25   ;;;;bin的个数，也是中间数

    bin  = fspan(sstmin-0.125, sstmax-0.125, toint(bin_num) +1)  ;;;bin最左边的数字
    olr_bin = new(toint(bin_num),float)
printMinMax(sst_clm_one, 1)

do i = 0,bin_num-1
   in := ind(sst_clm_one.ge.bin(i).and.sst_clm_one.lt.(bin(i)+0.25))

   olr_bin(i) = dim_avg(rc_one(in))
end do 
   

    x = bin(1:) - 0.125
    m1 = 27.5
    m2 = 29.0
    m3 = 29.5

    ind1 = ind(x.eq.m1)
    ind2 = ind(x.eq.m2)
    ind3 = ind(x.eq.m3)

    rc1 = regline(fspan(m1,m2,ind2-ind1+1), olr_bin(ind1:ind2))
    print(rc1@yintercept)

    rc2 = regline(fspan(m2,m3,ind3-ind2+1), olr_bin(ind2:ind3))
   

    print(rc2@yintercept)
    print (rc1+"*x"+"+" +rc1@yintercept)
    print (rc2+"*x"+"+" +rc2@yintercept)

;;;;draw plot 
    wks=gsn_open_wks("pdf", "SST-OLR-regression-scatter-plot-mam")


    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False


    res@xyLineColor     = (/"red"/)          ; change line color
    res@xyMonoLineColor = True

    res@xyMarkerColor ="gray30"
    res@xyMarkLineModes   = (/"Markers","Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=1
    ; res@xyExplicitLegendLabels = (/"", "mon_change", "Linear trend"/) 
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.001     ; Marker size (default 0.01)
   
    res@trYReverse =  True
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=30.0
    res@trYMinF=-30.0

    
    res@tmYLValues    = ispan(-30,30, 5)
    res@tmYLMinorValues = ispan(-30,30, 1)
    res@tmYLLabels    = ispan(-30,30, 5)
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=31.0
    res@trXMinF=26.0
    res@tmXBValues    = fspan(26.0,31.0,11)
    res@tmXBLabels    = fspan(26.0,31.0,11)


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
    
    res@tiMainString="~F22~SST and OLR"

    res@tiYAxisString="~F22~OLR [W/m^2]"
    res@tiYAxisFontHeightF=0.02
    
    res@tiXAxisString="~F22~SST [~S~o~N~C]"
    res@tiXAxisFontHeightF=0.02
    plot = new(2, graphic)
    

    plot(0) =gsn_csm_xy(wks, sst_clm_one, rc_one, res)


    plres=True
    plres@gsLineThicknessF=3.0
    plres@gsLineColor = "red"
    plres@gsLineDashPattern=0

    printVarSummary(bin)
    printVarSummary(olr_bin)
    plot_line=gsn_add_polyline(wks, plot(0), bin(1:) - 0.125, olr_bin, plres)

    mkres = True
    mkres@gsMarkerColor = "blue"
    mkres@gsMarkerIndex = 1
    mkres@gsMarkerSizeF = 0.02

    plot_mk=gsn_add_polymarker(wks, plot(0), bin(1:) - 0.125, olr_bin, mkres)
 

 draw(plot(0))
 frame(wks)


 end 




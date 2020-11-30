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

; f2 = addfile("/home/yangsong3/data-observation/SST/NOAA_OISST/sst.mnmean.nc", "r")
; date:=cd_calendar(f2->time, 1)
; timeind :=ind(date.ge.startmon.and.date.le.endmon)
; sst = short2flt(f2->sst(timeind,{lat1:lat2+1.0},{lon1-1.0:lon2+1.0})) ;;;2.5 ×2.5 grid
; printVarSummary(sst)
; ;;;;


  f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
  date:=cd_calendar(f2->time, 1)

  timeind:=ind(date.le.201512.and.date.ge.197901)

  sst_r =  lonFlip(f2->sst(timeind,:,:))
  
  sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

  printVarSummary(sst)



;;;;chazhi 
   print("**********start interp************")
    ; xi = sst&lon
    ; yi = sst&lat

    xi = sst&longitude
    yi = sst&latitude

    ;   xi = sst&lon
    ; yi = sst&lat


    xo = olr&lon
    yo = olr&lat



    sst_chazhi = linint2_Wrap(xi, yi, sst, True, xo, yo, 0)

    sst_chazhi!1= "lat"
    sst_chazhi!2= "lon"

    sst_chazhi&lat = yo
    sst_chazhi&lon = xo
;;;;
   
   
;;;;
    printVarSummary(sst_chazhi)
     olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)



;;;;Climatology bin
    ; sst_clm = dim_avg_n_Wrap(sst_chazhi(:,{0:20},{125:170}), 0)   
    ; olr_clm = dim_avg_n_Wrap(olr(:,{0:20},{125:170}), 0)
    


    ; sst_clm_MAM = dim_avg_n_Wrap( month_to_season(sst_chazhi(:,{0:20},{125:170}),"MAM"), 0)   
    ; olr_clm_MAM = dim_avg_n_Wrap(month_to_season(olr(:,{0:20},{125:170}),"MAM"), 0)

    ; sst_clm_JJA = dim_avg_n_Wrap( month_to_season(sst_chazhi(:,{0:20},{125:170}),"JJA"), 0)   
    ; olr_clm_JJA = dim_avg_n_Wrap(month_to_season(olr(:,{0:20},{125:170}),"JJA"), 0)
    
    sst_clm = dim_avg_n_Wrap( month_to_season(sst_chazhi(:,{0:20},{125:170}),"SON"), 0)   
    olr_clm = dim_avg_n_Wrap(month_to_season(olr(:,{0:20},{125:170}),"SON"), 0)

    ; sst_clm_DJF = dim_avg_n_Wrap( month_to_season(sst_chazhi(:,{0:20},{125:170}),"DJF"), 0)   
    ; olr_clm_DJF = dim_avg_n_Wrap(month_to_season(olr(:,{0:20},{125:170}),"DJF"), 0)

    ; sst_mam_jja = array_append_record(sst_clm_MAM,sst_clm_JJA,0)
    ; sst_son_djf  = array_append_record(sst_clm_SON,sst_clm_DJF,0)
    
    ; sst_clm = array_append_record(sst_mam_jja,sst_son_djf,0)

    ; olr_mam_jja = array_append_record(olr_clm_MAM,olr_clm_JJA,0)
    ; olr_son_djf  = array_append_record(olr_clm_SON,olr_clm_DJF,0)
    
    ; olr_clm = array_append_record(olr_mam_jja,olr_son_djf,0)
    
  

    sst_clm_one = ndtooned(sst_clm)
    olr_clm_one = ndtooned(olr_clm)

    mono = 1                            ; ascending=1 , descending=-1
    ii   =  dim_pqsort_n(sst_clm_one,mono,0)
    xx   = sst_clm_one(ii)                        ; ascending order 
    yy  = olr_clm_one(ii)
 
    n = 4
    c1 = lspoly(xx, yy, 1,n)    ; all weights are set to one 
 
    Y1 = c1(0) + c1(1)*xx + c1(2)*xx^2 + c1(3)*xx^3 
;;;

;;;;Monthly bin

   sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
   olr = (/runave_n_Wrap(olr, 3, 0, 0)/)

   ; dim = dimsizes(sst_chazhi)
   
   ; sst_chazhi_june = sst_chazhi(11:dim(0)-1:12,:,:)
   ; sst_chazhi_july = sst_chazhi(0:dim(0)-1:12,:,:)
   ; sst_chazhi_aug = sst_chazhi(1:dim(0)-1:12,:,:)

   ; sst_chazhi_jj =  array_append_record(sst_chazhi_june,sst_chazhi_july,0)
   ; sst_chazhi_jja = array_append_record(sst_chazhi_jj , sst_chazhi_aug,0)


   ;  olr_june = olr(11:dim(0)-1:12,:,:)
   ;  olr_july = olr(0:dim(0)-1:12,:,:)
   ;  olr_aug = olr(1:dim(0)-1:12,:,:)

   ;  olr_jj =  array_append_record(olr_june,olr_july,0)
   ;  olr_jja = array_append_record(olr_jj , olr_aug,0)
    
   ;  sst_one = ndtooned(sst_chazhi_jja(:,{0:15},{125:160}))
   ;  olr_one = ndtooned(olr_jja(:,{0:15},{125:160}))

    sst_one = ndtooned(sst_chazhi(:,{0:20},{100:170}))
    olr_one = ndtooned(olr(:,{0:20},{100:170}))

    ; sst_one = ndtooned(sst_chazhi(:,{-5:30},{100:170}))
    ; olr_one = ndtooned(olr(:,{-5:30},{100:170}))   

    jj   =  dim_pqsort_n(sst_one,mono,0)

    xx2 = sst_one(jj)
    yy2 = olr_one(jj)

    n = 3
    c2 = lspoly(xx2, yy2, 1,n)    ; all weights are set to one 

    Y2 = c2(0) + c2(1)*xx2 + c2(2)*xx2^2 
    ;;;


;;;bin 
;    bin_num =  (30.0-27.5)/0.25
;    bin  = fspan(27.5-0.125, 30.0 - 0.125, toint(bin_num) +1)   ;;; 27.5左右 为一个bin
;    olr_bin = new(toint(bin_num),float)
;    printMinMax(sst_clm_one, 1)

; do i = 0,bin_num-1
;    in := ind(sst_clm_one.ge.bin(i).and.sst_clm_one.lt.(bin(i)+0.25))

;    olr_bin(i) = dim_avg(olr_clm_one(in))
; end do 

; dolr_dsst = center_finite_diff(olr_bin, bin(1:) -0.125, False, 0)




;;;bin SON
   bin_num =  (30.0-28.5)/0.25
   bin  = fspan(28.5-0.125, 30.0 - 0.125, toint(bin_num) +1)   ;;; 27.5左右 为一个bin
   olr_bin = new(toint(bin_num),float)
   printMinMax(sst_clm_one, 1)

do i = 0,bin_num-1
   in := ind(sst_clm_one.ge.bin(i).and.sst_clm_one.lt.(bin(i)+0.25))

   olr_bin(i) = dim_avg(olr_clm_one(in))
end do 

dolr_dsst = center_finite_diff(olr_bin, bin(1:) -0.125, False, 0)



    rc1 = regline(fspan(28.5,29.0,3), dolr_dsst(0:2))
    print(rc1@yintercept)

    rc2 = regline(fspan(29.0,29.75,toint((29.75-29.0)/0.25)+1), dolr_dsst(2:5))

    print(rc2@yintercept)
    print (rc1+"*x"+"+" +rc1@yintercept)
    print (rc2+"*x"+"+" +rc2@yintercept)
exit()
    
;;;;draw plot 
    wks=gsn_open_wks("pdf", "SST-OLR-relation-clmdata-TWNP-HADISSTDATA-SON")


    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False


    res@xyLineColor     = (/"red"/)          ; change line color
    res@xyMonoLineColor = True

    ; res@xyMarkerColor ="gray30"
    res@xyMarkerColor ="black"

    res@xyMarkLineModes   = (/"Markers","Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=1
    ; res@xyExplicitLegendLabels = (/"", "mon_change", "Linear trend"/) 
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.005     ; Marker size (default 0.01)
   
    res@trYReverse =  True
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=280.0
    res@trYMinF=210.0

    
    res@tmYLValues    = fspan(280,210,8)
    res@tmYLMinorValues = fspan(280,210,8)
    res@tmYLLabels    = fspan(280,210,8)
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=31.0
    res@trXMinF=25.0
    res@tmXBValues    = fspan(25.0,31.0,13)
    res@tmXBLabels    = fspan(25.0,31.0,13)


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
    

    plot(0) =gsn_csm_xy(wks, sst_clm_one, (/olr_clm_one/), res)


    plres=True
    plres@gsLineThicknessF=3.0
    plres@gsLineColor = "red"
    plres@gsLineDashPattern=0

    printVarSummary(bin)
    printVarSummary(olr_bin)
    plot_line=gsn_add_polyline(wks, plot(0), bin(1:) - 0.25/2, olr_bin, plres)
    ; plot_line=gsn_add_polyline(wks, plot(0), bin(1:) , olr_bin, plres)

    mkres = True
    mkres@gsMarkerColor = "blue"
    mkres@gsMarkerIndex = 1
    mkres@gsMarkerSizeF = 0.02

    plot_mk=gsn_add_polymarker(wks, plot(0), bin(1:) - 0.25/2, olr_bin, mkres)
    ; plot_mk=gsn_add_polymarker(wks, plot(0), bin(1:) , olr_bin, mkres)



    res@tmYLValues    := fspan(30,-40,15)
    res@tmYLMinorValues := fspan(30,-40,15)
    res@tmYLLabels    := fspan(30,-40,15)
    res@trYMaxF=30.0
    res@trYMinF=-40.0
    res@xyMarkerColor ="black"
    res@xyMarkLineModes   := (/"Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=0

    res@tiMainString="~F34~6~F22~(OLR)/~F34~6~F22~(SST)"
    res@tiYAxisString="~F22~Chage rate [W/m^2/~S~o~N~C]"
    res@gsnYRefLine = 0.0
    res@gsnYRefLineThicknessF = 1.0
    plot(1) =gsn_csm_xy(wks, bin(1:) - 0.25/2, (/dolr_dsst/), res)


 pres  = True
 gsn_panel(wks, plot, (/1,2/), pres)

 end 




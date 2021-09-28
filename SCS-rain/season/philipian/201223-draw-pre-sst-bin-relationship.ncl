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
lat1 = -5
lat2 =  30
lon1 =  120
lon2 = 170

; startmon = 198201
; endmon = 201912

startmon = 197901
endmon = 201512



f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/gpcp_pre_chazhi_1979-2015.nc", "r")
date:=cd_calendar(f1->time, 1)
timeind :=ind(date.ge.startmon.and.date.le.endmon)
olr = short2flt(f1->p_chazhi(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

printVarSummary(olr)

; f2 = addfile("/home/yangsong3/data-observation/SST/NOAA_OISST/sst.mnmean.nc", "r")
; date:=cd_calendar(f2->time, 1)
; timeind :=ind(date.ge.startmon.and.date.le.endmon)
; sst = short2flt(f2->sst(timeind,{lat1:lat2+1.0},{lon1-1.0:lon2+1.0})) ;;;2.5 ×2.5 grid
; printVarSummary(sst)
;;;;


  f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
  date:=cd_calendar(f2->time, 1)

  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  sst_r =  lonFlip(f2->sst(timeind,:,:))
  
  sst = sst_r(:,{lat1:lat2},{lon1:lon2})

  printVarSummary(sst)


  sst_chazhi = sst

   
   
;;;;
    printVarSummary(sst_chazhi)
     olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)
;;;;

;;;;Monthly bin

   sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
   olr = (/runave_n_Wrap(olr, 3, 0, 0)/)

   dim = dimsizes(sst_chazhi)
   
   sst_chazhi_june = sst_chazhi(8:dim(0)-1:12,:,:)
   sst_chazhi_july = sst_chazhi(9:dim(0)-1:12,:,:)
   sst_chazhi_aug = sst_chazhi(10:dim(0)-1:12,:,:)

   sst_chazhi_jj =  array_append_record(sst_chazhi_june,sst_chazhi_july,0)
   sst_chazhi_jja = array_append_record(sst_chazhi_jj , sst_chazhi_aug,0)


    olr_june = olr(2:dim(0)-1:12,:,:)
    olr_july = olr(3:dim(0)-1:12,:,:)
    olr_aug = olr(4:dim(0)-1:12,:,:)

    olr_jj =  array_append_record(olr_june,olr_july,0)
    olr_jja = array_append_record(olr_jj , olr_aug,0)
    
    sst_one = ndtooned(sst_chazhi_jja(:,{5},{145}))
    olr_one = ndtooned(olr_jja(:,{5},{145}))

    ; sst_one = ndtooned(sst_chazhi(:,{-20:20},{125:160}))
    ; olr_one = ndtooned(olr(:,{-20:20},{125:160}))


;;;bin 
   
   sstmax = 30.0
   sstmin = 28.5
   ; sstmax = 30.5
   ; sstmin = 26.5

   bin_num =  (sstmax- sstmin)/0.25   ;;;;bin的个数，也是中间数

    bin  = fspan(sstmin-0.125, sstmax-0.125, toint(bin_num) +1)  ;;;bin最左边的数字
    olr_bin = new(toint(bin_num),float)
printMinMax(sst_one, 1)

do i = 0,bin_num-1
   in := ind(sst_one.ge.bin(i).and.sst_one.lt.(bin(i)+0.25))

   olr_bin(i) = dim_avg(olr_one(in))
end do 

dolr_dsst = center_finite_diff(olr_bin, bin(1:) - 0.125, False, 0)

    x = bin(1:) - 0.125

    ; ind1 = ind(x.eq.26.5)
    ; ind2 = ind(x.eq.28.25)
    ; ind3 = ind(x.eq.29.25)

    ; rc1 = regline(fspan(26.5,28.25,ind2-ind1+1), dolr_dsst(ind1:ind2))
    ; print(rc1@yintercept)

    ; rc2 = regline(fspan(28.25,29.25,ind3-ind2+1), dolr_dsst(ind2:ind3))

    ; print(rc2@yintercept)
    ; print (rc1+"*x"+"+" +rc1@yintercept)
    ; print (rc2+"*x"+"+" +rc2@yintercept)




;;;;draw plot 
    wks=gsn_open_wks("pdf", "SST-pre-relathionship-test")


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
    res@xyMarkerSizeF     = 0.00001     ; Marker size (default 0.01)
   
    res@trYReverse =  False
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=20.0
    res@trYMinF=0.0

    
    res@tmYLValues    = fspan(0,20,11)
    res@tmYLMinorValues = fspan(0,20,21)
    res@tmYLLabels    = fspan(0,20,11)
    
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
    

    plot(0) =gsn_csm_xy(wks,sst_one, (/olr_one/), res)


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


    res@tmYLValues    := ispan(-5,5,1)
    res@tmYLMinorValues := ispan(-5,5,1)
    res@tmYLLabels    := ispan(-5,5,1)
    res@trYMaxF=5.0
    res@trYMinF=-5.0
    res@xyMarkerColor ="black"
    res@xyMarkLineModes   := (/"Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=0

    res@tiMainString="~F34~6~F22~(OLR)/~F34~6~F22~(SST)"
    res@tiYAxisString="~F22~Chage rate [W/m^2/~S~o~N~C]"
    res@gsnYRefLine = 0.0
    res@gsnYRefLineThicknessF = 1.0
    plot(1) =gsn_csm_xy(wks, bin(1:) -0.125, (/dolr_dsst/), res)


 pres  = True
 gsn_panel(wks, plot, (/1,2/), pres)

 end 




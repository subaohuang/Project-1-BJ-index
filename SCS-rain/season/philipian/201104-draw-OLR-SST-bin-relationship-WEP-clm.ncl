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

startmon = 198201
endmon = 201912


f1 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
date:=cd_calendar(f1->time, 1)
timeind :=ind(date.ge.startmon.and.date.le.endmon)
olr = short2flt(f1->olr(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

printVarSummary(olr)

f2 = addfile("/home/yangsong3/data-observation/SST/NOAA_OISST/sst.mnmean.nc", "r")
date:=cd_calendar(f2->time, 1)
timeind :=ind(date.ge.startmon.and.date.le.endmon)
sst = short2flt(f2->sst(timeind,{lat1:lat2+1.0},{lon1-1.0:lon2+1.0})) ;;;2.5 ×2.5 grid
printVarSummary(sst)
;;;;


;;;;chazhi 
   print("**********start interp************")
    xi = sst&lon
    yi = sst&lat

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


;;;;Climatology bin
    sst_clm = dim_avg_n_Wrap(sst_chazhi(:,{0:10},{125:160}), 0)   
    olr_clm = dim_avg_n_Wrap(olr(:,{0:10},{125:160}), 0)

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
    
   sst_one = ndtooned(sst_chazhi(:,{0:15},{125:160}))
   olr_one = ndtooned(olr(:,{0:15},{125:160}))

    jj   =  dim_pqsort_n(sst_one,mono,0)

    xx2 = sst_one(jj)
    yy2 = olr_one(jj)

    n = 4
    c2 = lspoly(xx2, yy2, 1,n)    ; all weights are set to one 
 
    Y2 = c2(0) + c2(1)*xx2 + c2(2)*xx2^2 + c2(3)*xx2^3 
;;;


;;;;bin 
   bin_num =  (29.75-28.0)/0.25 

   bin  = fspan(28.0, 29.75, 8)
   olr_bin = new(toint(bin_num),float)
   
do i = 0,bin_num-1
   in := ind(sst_clm_one.ge.bin(i).and.sst_clm_one.lt.(bin(i)+0.25))
   olr_bin(i) = dim_avg(olr_clm_one(in))
end do 




 
;;;;draw plot 
    wks=gsn_open_wks("pdf", "SST-OLR-relation-Climatology")


    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False


    res@xyLineColor     = (/"red"/)          ; change line color
    res@xyMonoLineColor = True

    res@xyMarkerColor ="black"
    res@xyMarkLineModes   = (/"Markers","Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=1
    ; res@xyExplicitLegendLabels = (/"", "mon_change", "Linear trend"/) 
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.005      ; Marker size (default 0.01)
   
    res@trYReverse =  True
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=270.0
    res@trYMinF=210.0

    
    res@tmYLValues    = fspan(270,210,7)
    res@tmYLMinorValues = fspan(270,210,16)
    res@tmYLLabels    = fspan(270,210,7)
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=30.0
    res@trXMinF=26.0
    res@tmXBValues    = fspan(26.0,30.0,9)
    res@tmXBLabels    = fspan(26.0,30.0,9)


    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
    res@tmXBLabelFontHeightF=0.015

    res@tmXBMajorLengthF=0.008        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0.008  ;调整刻度的长度，注意2个命令要同时使用
    res@tmXUseBottom=False
    res@tmXBLabelDeltaF=-0.8

    res@tmYLLabelFontHeightF=0.015
    res@vpXF=0.15
    res@vpYF=0.9
    ; res@vpHeightF=0.6
    ; res@vpWidthF=0.6
    
    res@tiMainString="Mean OLR in each SST bin (from climatological data)"
    res@tiYAxisString="~F22~OLR [W/m^2]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~SST [~S~o~N~C]"
    res@tiXAxisFontHeightF=0.012
    plot = new(2, graphic)
    
   plot(0) =gsn_csm_xy(wks, sst_clm_one, (/olr_clm_one/), res)


  ;  plot(1) =gsn_csm_xy(wks, xx2, (/yy2,Y2/), res)

    
    plres=True
    plres@gsLineThicknessF=2.0
    plres@gsLineColor = "blue"
    plres@gsLineDashPattern=0
    plres@xyMarkLineMode = "line"


    printVarSummary(bin)
    printVarSummary(olr_bin)
    plot_line=gsn_add_polyline(wks, plot(0), bin(1:) , olr_bin, plres)

    mkres = True
    mkres@gsMarkerColor = "red"
    mkres@gsMarkerIndex = 1
    mkres@gsMarkerSizeF = 0.02

    plot_mk=gsn_add_polymarker(wks, plot(0), bin(1:), olr_bin, mkres)


    ; 90~S~o~N~E

    draw(plot)
    frame(wks)

 end 




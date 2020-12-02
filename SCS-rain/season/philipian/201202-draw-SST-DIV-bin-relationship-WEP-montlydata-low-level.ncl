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

startmon = 197901
endmon = 201612


f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
date:=cd_calendar(f2->time, 1)

timeind:=ind(date.le.endmon.and.date.ge.startmon)

sst_r =  lonFlip(f2->sst(timeind,:,:))

sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

; sst = (/dtrend_msg_n(ispan(1,dimsizes(sst&time),1),sst,False,True,0)/)

printVarSummary(sst)



fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
f3 = addfiles(fils3, "r")

date:=cd_calendar(f3[:]->time, 1)
timeind:=ind(date.ge.startmon.and.date.le.endmon)
time = f3[:]->time(timeind)
lat = f3[1]->latitude
lon = f3[1]->longitude


div_r = new((/dimsizes(time),7, dimsizes(lat),dimsizes(lon)/),float)
  

do i = 0,dimsizes(timeind)/12-1

  div_r(i*12:(i+1)*12-1,:,:,:) = short2flt(f3[i]->d(:,{1000:850},:,:))
  
end do
;


div = dim_avg_n_Wrap(div_r(:,:,{-10:30},{100:180}),1)

div = (/div*10^5/)
div = (/-1*div/)


;;;;


;;;;chazhi 
   print("**********start interp************")
    xi = div&longitude
    yi = div&latitude

    xo = sst&longitude
    yo = sst&latitude


    div_chazhi = linint2_Wrap(xi, yi, div, True, xo, yo, 0)

    div_chazhi!1= "lat"
    div_chazhi!2= "lon"

    div_chazhi&lat = yo
    div_chazhi&lon = xo
;;;;
   
   
;;;;
    printVarSummary(div_chazhi)
    sst = where(ismissing(div_chazhi), sst@_FillValue, sst)



;;;;Climatology bin
    div_clm = dim_avg_n_Wrap(div_chazhi(:,{0:15},{125:160}), 0)   
    sst_clm = dim_avg_n_Wrap(sst(:,{0:15},{125:160}), 0)


    div_clm_one = ndtooned(div_clm)
    sst_clm_one = ndtooned(sst_clm)

;;;

;;;;Monthly bin

   div_chazhi = (/runave_n_Wrap(div_chazhi, 3, 0, 0)/)
   sst = (/runave_n_Wrap(sst, 3, 0, 0)/)

   ; dim = dimsizes(div_chazhi)
   
   ; div_chazhi_june = div_chazhi(11:dim(0)-1:12,:,:)
   ; div_chazhi_july = div_chazhi(0:dim(0)-1:12,:,:)
   ; div_chazhi_aug = div_chazhi(1:dim(0)-1:12,:,:)

   ; div_chazhi_jj =  array_append_record(div_chazhi_june,div_chazhi_july,0)
   ; div_chazhi_jja = array_append_record(div_chazhi_jj , div_chazhi_aug,0)


   ;  sst_june = sst(11:dim(0)-1:12,:,:)
   ;  sst_july = sst(0:dim(0)-1:12,:,:)
   ;  sst_aug = sst(1:dim(0)-1:12,:,:)

   ;  sst_jj =  array_append_record(sst_june,sst_july,0)
   ;  sst_jja = array_append_record(sst_jj , sst_aug,0)
    
   ;  div_one = ndtooned(div_chazhi_jja(:,{0:15},{125:160}))
   ;  sst_one = ndtooned(sst_jja(:,{0:15},{125:160}))

    div_one = ndtooned(div_chazhi(:,{0:20},{125:160}))
    sst_one = ndtooned(sst(:,{0:20},{125:160}))



    ; div_one = ndtooned(div_chazhi(:,{-5:30},{100:170}))
    ; sst_one = ndtooned(sst(:,{-5:30},{100:170}))   

;;;bin month
    bin_num =  (30.5-25.0)/0.25   ;;;;bin的个数，也是中间数
    
    bin_month  = fspan(25.0-0.125, 30.5-0.125, toint(bin_num) +1)  ;;;bin最左边的数字
    div_bin_month = new(toint(bin_num),float)
    printMinMax(sst_one, 1)

do i = 0,bin_num-1
   in := ind(sst_one.ge.bin_month(i).and.sst_one.lt.(bin_month(i)+0.25))

   div_bin_month(i) = dim_avg(div_one(in))
end do 

   ddiv_dsst = center_finite_diff(div_bin_month, bin_month(1:) - 0.125, False, 0)


   ddiv_dsst  = (/ddiv_dsst*10.0/)



;    bin_num =  (30.75-26.75)/0.25

;    bin  = fspan(26.75, 30.75, toint(bin_num) +1)
;    sst_bin = new(toint(bin_num),float)

; do i = 0,bin_num-1
   
;    in := ind(div_one.ge.bin(i).and.div_one.lt.(bin(i)+0.25))
  
;    sst_bin(i) = dim_avg(sst_one(in))
; end do 
; dsst_ddiv = center_finite_diff(sst_bin, bin(1:), False, 0)

; rc1 = regline(fspan(28.0,28.75,4), dsst_ddiv(4:7))
; print(rc1@yintercept)

; rc2 = regline(fspan(28.75,30.25,toint((30.25-28.75)/0.25)+1), dsst_ddiv(7:13))

; print(rc2@yintercept)
; print (rc1+"*x"+"+" +rc1@yintercept)
; print (rc2+"*x"+"+" +rc2@yintercept)
; rc1 = regline(fspan(28.0,30.25,10), dsst_ddiv(4:13))
; print(rc1@yintercept)
; print (rc1+"*x"+"+" +rc1@yintercept)



 ; exit()
;;;;draw plot 
    ; wks=gsn_open_wks("pdf", "div-sst-relation-monlydata-850hpa")
    wks=gsn_open_wks("pdf", "div-sst-relation-clmdata-low-level(850-1000)")


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
   
    res@trYReverse =  False
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=.5
    res@trYMinF=-.5

    
    res@tmYLValues    = ispan(-5,5,1)*0.1
    ; res@tmYLMinorValues = ispan(-5,5,1)*0.1
    res@tmYLLabels    = ispan(-5,5,1)*0.1
    
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
    
    res@tiMainString="(a) ~F25~mean DIV in each SST bin"

    res@tiYAxisString="~F25~DIV [W/m~S~2~N~]"
    res@tiYAxisFontHeightF=0.02
    
    res@tiXAxisString="~F25~SST [~S~o~N~C]"
    res@tiXAxisFontHeightF=0.02
    plot = new(2, graphic)
    

    plot(0) =gsn_csm_xy(wks, sst_one, (/div_one/), res)


    plres=True
    plres@gsLineThicknessF=3.0
    plres@gsLineColor = "red"
    plres@gsLineDashPattern=0

 
    plot_line=gsn_add_polyline(wks, plot(0), bin_month(1:) - 0.125,div_bin_month , plres)

    mkres = True
    mkres@gsMarkerColor = "blue"
    mkres@gsMarkerIndex = 1
    mkres@gsMarkerSizeF = 0.02

    plot_mk=gsn_add_polymarker(wks, plot(0), bin_month(1:) -0.125, div_bin_month, mkres)


    res@tmYLValues    := fspan(-2,2,9)
    ; res@tmYLMinorValues := fspan(-2,2,9)
    res@tmYLLabels    := fspan(-2,2,9)
    res@trYMaxF=2.0
    res@trYMinF=-2.0



    res@xyMarkerColor ="black"
    res@xyMarkLineModes   := (/"Lines"/)
    res@xyMonoLineThickness = True
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=0

    res@tiMainString="~F34~6~F22~(DIV)/~F34~6~F22~(SST)"
    res@tiYAxisString="~F22~Chage rate [W/m^2/~S~o~N~C]"
    res@gsnYRefLine = 0.0
    res@gsnYRefLineThicknessF = 1.0
    plot(1) =gsn_csm_xy(wks, bin_month(1:) - 0.125, (/ddiv_dsst/), res)


 pres  = True
 gsn_panel(wks, plot, (/1,2/), pres)

 end 




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
function trend(data)
begin
  rc1 =regCoef_n(ispan(1, dimsizes(data&time),1), data,0,0)
  rc1 = rc1*dimsizes(data&time)
  copy_VarCoords(data(1,:,:), rc1)
  return rc1
end 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

begin

yearstart = 198001
yearend = 202012

lat1 = -50
lat2 = 50
lon1 = 0
lon2 = 360

;;;;read data
var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc"/)

diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"

f1 = addfile(diri+var(0), "r")
f2 = addfile(diri+var(1), "r")
f3 = addfile(diri+var(2), "r")

time = f1->time 
date = cd_calendar(time, 1)
delete(time)

timeind :=ind(date.le.yearend.and.date.ge.yearstart)

t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

t=t+273.15

q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
q = q/1000.0  ;;;convert to kg/kg

hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;

Cp = 1004
Lv = 2.5*10^6
g = 9.8


;;;MSE 
 plev = t&level
 ptop = min(plev)
 psfc = 1013.0
 dz = plev
 dz = dpres_plevel(plev, psfc, ptop, 0)
 dz&level = plev
 printVarSummary(dz)
 

 A = t
 A =Cp*t

 B = q
 B = Lv*q
;;;;


;;;;;;;;;;;;;;;;;;;
 A_col =  dim_sum_wgt_n_Wrap(A(:,{1000:500},:,:), dz({1000:500}), 1, 1)
 A_col = (/rmMonAnnCycTLL(A_col)/)

 A_col_12 = split_12(A_col)
 A_col_12 = (/A_col_12/10^6/)
;;;;;


;;;;;;;;;;;;;;;;;;;
 B_col =  dim_sum_wgt_n_Wrap(B(:,{1000:500},:,:), dz({1000:500}), 1, 1)
 B_col = (/rmMonAnnCycTLL(B_col)/)

 B_col_12 = split_12(B_col)
 B_col_12 = (/B_col_12/10^6/)
 ;;;;;
 A_col_annual = month_to_season(A_col, "MAM")
 B_col_annual = month_to_season(B_col, "MAM")

 printVarSummary(A_col_annual)
 A_col_annual = (/dim_avg_n_Wrap(A_col_12(0:4,:,:,:), 0)/)
 B_col_annual = (/dim_avg_n_Wrap(B_col_12(0:4,:,:,:), 0)/)



 rc_A_mam = trend(A_col_annual)
 
 rc_B_mam = trend(B_col_annual)

 prob_A = rc_test(rc_A_mam)
 prob_B = rc_test(rc_B_mam)
 printVarSummary(rc_A_mam)


  wks = gsn_open_wks("pdf","CpT和Lq的趋势1980-2020")          ; send graphics to PNG file

;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=0.01
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn= False
  ; res@lbLabelBarOn=True

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -15
  res@mpMaxLatF             = 20
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  

 
  res@tmYLValues = (/-10,0,10,20/) 
  res@tmYLLabels = (/"10S","0","10N","20N"/) 

  res@tmYLMinorValues = ispan(-15,20,5)


  res@mpFillOn = True
  ; res@mpLandFillColor = "gray"
  ; res@mpOceanFillColor = "Transparent"
  ; res@mpMonoFillColor = True
  res@mpFillDrawOrder = "PostDraw"


 
  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.01
  res@tmYLLabelFontHeightF=0.01
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.6
  res@cnMaxLevelValF       =   1.6
  res@cnLevelSpacingF      =  0.4

  ; res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour level
  ; res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)

  res@lbLabelFontHeightF = 0.01
  res@lbLabelStride = 1
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueDarkRed18"
  res@cnFillPalette = "BlueWhiteOrangeRed"

;   res@mpShapeMode = "FreeAspect"
;   res@vpHeightF = 0.2
;   res@vpWidthF = 0.7
; ; 
  ; res@gsnLeftString="Annual mean MSE trend " 
  plot = new(2,graphic)

  res@gsnLeftString="(a) JFMAM 1000-300hPa integrated CpT trend " 
  plot(0) = gsn_csm_contour_map(wks, rc_A_mam, res)
  res@gsnLeftString="(b) JFMAM 1000-500hPa integrated Lq trend " 
  plot(1) = gsn_csm_contour_map(wks, rc_B_mam, res)

  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.0
  pgres@gsLineDashPattern=0

  ;;;画box
  lat1 = 0
  lat2 = 20
  lon1 = 250
  lon2 = 280

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;
  add  = ispan(1,2,1)+""
  do i=0,1
  plot@$add(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do
  ;;;打点

  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnFillOn = True

  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLevels = (/0,0.1/)
  res_pro@lbLabelBarOn = False
  res_pro@cnMonoFillColor = False
  res_pro@cnFillColors = (/"Transparent","White","Transparent"/)


  res_pro@cnLinesOn=False
  ; res_pro@cnMonoLineColor =False
  ; res_pro@cnLineColors = ((/"Transparent","black","Transparent"/))

  res_pro@cnMonoFillPattern = True
  res_pro@cnFillPattern = 17
  res_pro@cnFillDotSizeF = 0.002
  
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True
   
  iplot = plot 
  iplot(0) = gsn_csm_contour(wks,prob_A,res_pro)
  iplot(1) = gsn_csm_contour(wks,prob_B,res_pro)
  
  do i=0,1
  overlay(plot(i), iplot(i))
  end do

    ;;  

  pres = True
  pres@gsnPanelLabelBar  = True
  pres@pmLabelBarOrthogonalPosF = -0.02
  ; pres@pmLabelBarParallelPosF = 0.02
  pres@pmLabelBarHeightF = 0.05
  pres@pmLabelBarWidthF = 0.6

  gsn_panel(wks, plot, (/2,1/), pres)

end




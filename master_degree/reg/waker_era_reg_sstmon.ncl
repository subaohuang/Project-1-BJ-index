begin 
diri="/media/sf_share/era/"
fils=systemfunc("ls "+diri+"era_u_sst_slp.nc")
f1=addfile(fils, "r")

sst=f1->sst(0:455,:,:)
sst_a=sst@scale_factor*sst+sst@add_offset-274.15
sst_a=doubletofloat(sst_a)
copy_VarAtts(sst, sst_a)
copy_VarCoords(sst, sst_a)
sst_a@units="°C"

sst_a!0="month"
sst_a!1="latitude"
sst_a!2="longitude"
sst_a&month=ispan(1,456,1)
sst_a&latitude=sst&latitude
sst_a&longitude=sst&longitude

sst_a&latitude@units = "degrees_north"
sst_a&longitude@units = "degrees_east"

sst_a_new = sst_a(latitude|:,longitude|:,month|:)

printVarSummary(sst_a_new)

;计算rc




 rc = regCoef(sst_a_new&month,sst_a_new)
rc!0="lat"
rc!1="lon"
rc&lat=sst_a_new&latitude
rc&lon=sst_a_new&longitude
rc&lat@units = "degrees_north"
rc&lon@units = "degrees_east"

printVarSummary(rc)
;---------------------------------------------------
;画图 
;-------------------------------------------------------
wks=gsn_open_wks("png", "/media/sf_share/era_sstrc_month")
res=True
res = True
res@gsnAddCyclic=False
res@gsnDraw=False
res@gsnFrame=False

; res@cnFillPalette  = "temp_diff_18lev"
 res@cnFillPalette="BlueDarkRed18"
res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
res@cnFillOn = True 
res@cnLinesOn = False
res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
res@cnMinLevelValF       =   -0.003
res@cnMaxLevelValF       =   0.003
res@cnLevelSpacingF      =  0.00025
res@lbOrientation = "Vertical"

res@mpLimitMode           = "LatLon"
res@mpMinLatF             = min(rc&lat)
res@mpMaxLatF             = max(rc&lat)
res@mpMinLonF             = min(rc&lon)
res@mpMaxLonF             = max(rc&lon)
res@mpCenterLonF          =avg(rc&lon)



plot = gsn_csm_contour_map(wks, rc, res)

draw(plot)
frame(wks)


end
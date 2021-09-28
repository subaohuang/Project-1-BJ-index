load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;***************************************

begin

;;;;
;;; read data from grib file
;;;;


outfile="/home/ys17-19/lsh/Project/ENSO/CESM/B/heating_profile"
CTRL_file="/home/ys17-19/cjx/all_output_data/B2000_CTRL_output/CTRL_climatology.nc"
fin_ctrl = addfile(CTRL_file,"r")
g_font=0.02
l_thick = 2

;;; get variables: temperature, longitude, latitude

heat_IN = fin_ctrl->Z3(0,:,:)
lat = fin_ctrl->lat
lon = fin_ctrl->lon
heat = heat_IN*0
copy_VarMeta(heat_IN,heat)
;printVarSummary(lat)
printVarSummary(heat)


latS = 0.0
latN = 20.0
lonW = 125.0
lonE = 160.0

latC = (latS+latN)*0.5
lonC = (lonW+lonE)*0.5
;modified heat

do i=0,dimsizes(lon)-1,1
do j=0,dimsizes(lat)-1,1
if (lon(i).ge.lonW .and. lon(i).le.lonE .and. lat(j).ge.latS .and. lat(j).le.latN) then
stat_lat=lat(j)
stat_lon=lon(i)
heat(j,i) = tofloat((1 - abs(stat_lat-latC)/abs(latN-latC))*(1 - abs(stat_lon-lonC)/abs(lonE-lonC))) 
end if
end do
end do

print(heat({latC},{lonC}))
;print(heat)
;;;;
;;; create plot (png)
;;;;


wks = gsn_open_wks("pdf",outfile)                 ; open wk station
gsn_define_colormap(wks,"precip_diff_12lev")

res                        = True         ; mods desired  
res@cnFillOn             = False               ; turn on color for contours
res@cnLinesOn            = True             ; turn off contour lines
res@cnLineLabelsOn       = True              ; turn off contour line labels
res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
res@cnMinLevelValF        = 0.1         ; set min contour level
res@cnMaxLevelValF        = 1.0          ; set max contour level
res@cnLevelSpacingF       = 0.1          ; set contour interval
res@mpFillDrawOrder = "PreDraw"
res@cnInfoLabelOn = False
res@cnLineLabelFormat="0@;*.1f"
;;  res@cnLevels      =(/-1.5,-1,-0.5,0,0.5,1,1.5/)
; res@cnFillColors =(/2,4,5,6,7,9,10,11,13/)


res@gsnFrame         = False
res@gsnDraw  = False

res@gsnSpreadColors      = True               ; use full color map
res@gsnSpreadColorStart  = 2               ; start at color 17
res@gsnSpreadColorEnd    = 14                ; end at color 200

;  res@gsnPanelLabelBar    = True   ;Turn on panel labelbar
res@gsnStringFontHeightF = g_font

res@gsnMaximize      = True       ; fill up the page
res@gsnAddCyclic   = True;False
res@gsnPaperOrientation = "portrait"
res@cnLineThicknessF= l_thick
res@lbLabelBarOn = True   ; turn off the label bar
res@lbOrientation          = "vertical"
res@lbLabelFontHeightF  = g_font              ; make labels smaller
res@lbLabelStride = 2 
res@cnLineDrawOrder = "PostDraw"
res@mpFillOn   = True     ; turn off gray continents
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor      = "black";"white" 
;res@mpFillDrawOrder="PostDraw"
;res@mpCenterLonF = 180.
res@mpMaxLatF  = 30
res@mpMinLatF  = -30
res@mpMaxLonF  = 180
res@mpMinLonF  = 60

;  res@mpShapeMode = "FreeAspect" 

;  res@vpWidthF = 0.8
;  res@vpHeightF = 0.6

;  res@tmXBMajorOutwardLengthF = -0.
;  res@tmYLMajorOutwardLengthF = -0.
;  res@tmYRMajorOutwardLengthF = -0.

;  res@tmXBMinorOutwardLengthF = -0.
;  res@tmYLMinorOutwardLengthF = -0.
;  res@tmYRMinorOutwardLengthF = -0.

; now change the size of the tickmark labels
res@tmXBTickSpacingF = 60 
res@tmYLTickSpacingF = 15
res@tmXBLabelFontHeightF = g_font
; res@tmXBLabelFontAspectF = 1.6
res@tmXBLabelFontThicknessF = 2.0
res@tmYLLabelFontHeightF = g_font
; res@tmYLLabelFontAspectF = 1.3
res@tmYLLabelFontThicknessF = 2.0
; res@tmYLLabelFontQuality =
;-----------------------------------------------------------------------------

res@gsnLeftString = ""
res@gsnRightString = "Additional Heating"
plot = gsn_csm_contour_map_ce(wks,heat,res) ; create a default plot
;x_loc=(/90,160,160,90,90/)
;y_loc=(/10,10,-10,-10,10/)
;lnres   =    True
;lnres@gsLineColor="black"
;lnres@gsLineThicknessF= l_thick
;lnres@gsLineDashPattern=0
;dum=gsn_add_polyline(wks,plot,x_loc,y_loc,lnres)

;polyres                   = True
;polyres@gsMarkerIndex     = 2          ; polymarker style
;polyres@gsMarkerSizeF     = 16.          ; polymarker size
;polyres@gsMarkerThicknessF     = 6.          ; polymarker size
;polyres@gsMarkerColor     = "black"          ; polymarker size
;plot@$unique_string("fxxk ncl")$    = gsn_add_polymarker(wks,plot,110,10,polyres)
draw(plot)
frame(wks)
delete([/latS,latN,lonW,lonE,latC,lonC,heat,res/])
end

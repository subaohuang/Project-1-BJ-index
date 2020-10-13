

begin
; load "$NCARG_ROOTb/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOTb/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOTb/ncarg/nclscripts/csm/contributed.ncl"

fname = "/mnt/f/share/MRE/2014.txt"

nclos=10
nrows=1500
data_r = asciiread(fname,(/nrows,nclos/),"double")               
; print(data_r(:,3))
; exit()
lon=data_r(:,1)
lat=data_r(:,2)
print(min(lon))
print(max(lon))
print(min(lat))
print(max(lat))
; exit()
; print(lon)
; exit()
data=data_r(:,7)
; print(data)
data@_FillValue=-999
data@_MissingValue=-999
data=where(data.le.1000, data, data@_FillValue)

; print("data="+min(data))
; print("data"+max(data))
; exit()

 xo=fspan(73,136,100)
 yo=fspan(4,54,50)
; t_new=natgrid(lon, lat, data, xo,yo ) ;natgrid插值,就近插值，注意xo与yo必须是double或者float型
  rscan = (/10, 5, 2/)
  opt  = False          
  t_new = obj_anal_ic(lon,lat,data,xo,yo,rscan, opt)
; ; print(t_new)
;

 t_new!0="lat"
 t_new!1="lon"
 t_new&lon = xo
 t_new&lat = yo
 t_new&lat@units = "degrees_north"
 t_new&lon@units = "degrees_east"
 ; t_new2=t_new(lat|:,lon|:)
 printVarSummary(t_new)
 ;**********************
 wks=gsn_open_wks("pdf","/mnt/f/share/era/data/ceshi")
  
  res =True

 res@gsnDraw          = False             ; so we can add poly stuff
 res@gsnFrame         = False
 res@gsnAddCyclic =False

 res@cnFillOn             = True               ; Fill contours
 res@cnLinesOn=False
 res@cnFillPalette        = "BkBlAqGrYeOrReViWh200"            ; choose colormap 
 res@gsnSpreadColorStart = 11
 res@gsnSpreadColorEnd = 161
 res@cnLineColor = "black"
 res@cnFillDrawOrder="PreDraw"
 res@cnInfoLabelOn=False
 res@mpDataSetName = "Earth..4"
 res@mpDataBaseVersion = "mediumRes"
 ;设置zoom范围
 ; res@mpOutlineDrawOrder="PreDraw"
   res@mpLimitMode           = "LatLon"
   res@mpMinLatF             =18
   res@mpMaxLatF             =max(t_new&lat)
   res@mpMinLonF             = min(t_new&lon)
   res@mpMaxLonF             = max(t_new&lon)
   res@mpCenterLonF   = (res@mpMinLonF+res@mpMaxLonF)/2.0

 ;画地图，覆盖contouf图
   res@mpFillAreaSpecifiers=(/"Land","china","taiwan","water"/)      ;确定需要填充的具体位置
   res@mpSpecifiedFillColors=(/"white","transparent","transparent","white"/) ;填充位置相应的具体颜色
   res@mpLandFillColor="white"   ;上面设置能完成，没什么用处
  
 ;上面画的没有轮廓，就是没有线，需要通过轮廓加入线
   res@mpOutlineOn = True; turn off outline 
   res@mpGeophysicalLineThicknessF= 4. ; double the thickness of geophysical boundaries 
   res@mpNationalLineThicknessF= 2. ; double the thickness of national boundaries 
   res@mpOutlineBoundarySets="NoBoundaries"
  res@mpGeophysicalLineThicknessF= 2.       ; double the thickness of geophysical boundaries
  res@mpNationalLineThicknessF= 2.         ; double the thickness of national boundaries
  ; delete(res@mpOutlineSpecifiers)          ; about to change size of attribute, must delete
  res@mpOutlineSpecifiers = (/"China:states"/)
  res@tmXTOn=False
  res@tmYROn=False

  plots = gsn_csm_contour_map(wks,t_new,res)

 draw(plots)
 frame(wks)

 end
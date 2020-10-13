

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
  gsn_define_colormap(wks,"WhViBlGrYeOrRe")      ; define a different color map
  res =True

 res@gsnDraw          = False             ; so we can add poly stuff
 res@gsnFrame         = False
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

  plots = gsn_csm_map(wks,res)


     arr =fspan(10,120,7)
     colors = (/10,30,38,48,56,66,74,94/)  ; marker colors, dimsizes must 
                                        ; be equal to dimsizes(arr)+1

  R   = data
                                        ; how to color the markers.
  labels = new(dimsizes(arr)+1,string)  ; Labels for legend.

  num_distinct_markers = dimsizes(arr)+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,dimsizes(R)/),typeof(lat),-999)
  lon_new = new((/num_distinct_markers,dimsizes(R)/),typeof(lon),-999)

 ;       
  do i = 0, num_distinct_markers-1
    if (i.eq.0) then
      indexes = ind(R.lt.arr(0))
      labels(i) = "x < " + arr(0)
    end if
    if (i.eq.num_distinct_markers-1) then
      indexes = ind(R.ge.max(arr))
      labels(i) = "x >= " + max(arr)
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then       
      indexes = ind(R.ge.arr(i-1).and.R.lt.arr(i))
      labels(i) = arr(i-1) + " <= x < " + arr(i)
    end if

   if (.not.any(ismissing(indexes))) then 
      npts_range = dimsizes(indexes)   ; # of points in this range.
      lat_new(i,0:npts_range-1) = lat(indexes)
      lon_new(i,0:npts_range-1) = lon(indexes)
    end if
    delete(indexes)            ; Necessary b/c "indexes" may be a different                          ; size next time.
  end do
   

  
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.

  txres               = True
  txres@txFontHeightF = 0.015

;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
; At the same time, draw a legend showing the meaning of the markers.
;

  xleg = (/0.07,0.07,0.32,0.32,0.57,0.57,0.82,0.82/)   ; Location of

  xtxt = (/0.16,0.16,0.41,0.41,0.66,0.66,0.91,0.91/)   ; legend markers

  yleg = (/0.22,0.17,0.22,0.17,0.22,0.17,0.22,0.17/)   ; and text

  ytxt = (/0.22,0.17,0.22,0.17,0.22,0.17,0.22,0.17/)   ; strings.



    do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerThicknessF = 0.7*(i+1)
      gsn_polymarker(wks,plots,lon_new(i,:),lat_new(i,:),gsres)
    ;
    ; Add marker and text for the legend.
    ;
      gsn_polymarker_ndc(wks, xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks,labels(i),xtxt(i),ytxt(i),txres)
    end if
    end do

 draw(plots)
 frame(wks)

 end
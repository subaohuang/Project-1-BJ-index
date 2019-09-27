;==========================================================
; User Specifications
;==========================================================
nlat = 128
mlon = 256
typeMask = "integer"     ; "byte" [default] or "integer"

grdType = 0              ; 0=gaussian , 1=regular , 2=regular offset
lonStrt = 0              ; 0 [GM] or -180 [Dateline]

diro = "./"
filo = "SampleLandSeaMask.nc"

PLOT = True              ; create a plot [ps]
;==========================================================
; End Required User Specifications
; Note: User may edit the mask below
;==========================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

function localByte2Int ( x:byte )
local dimx, X
begin
  delete(x@_FillValue)

  dimx = dimsizes(x)
  X    = new ( dimx, "integer", "No_FillValue")
  X    = x
  return(X)
end

begin

  if (grdType.eq.0) then      ; Gaussian [ S->N, GM eastward]
      lat  = latGau(nlat, "lat", "latitude", "degrees_north")
      lon  = lonGlobeF(mlon, "lon", "longitude", "degrees_east")
      gw   = latGauWgt(nlat, "lat", "gaussian weights", "")
  end if

  if (grdType.eq.1) then      ; Regular  [ 90S->90N, GM eastward]
      lat  = latGlobeF(nlat, "lat", "latitude", "degrees_north")
      lon  = lonGlobeF(mlon, "lon", "longitude", "degrees_east")
  end if

  if (grdType.eq.2) then      ; Regular offset from pole and GM
      lat  = latGlobeFo(nlat, "lat", "latitude", "degrees_north")
      lon  = lonGlobeFo(mlon, "lon", "longitude", "degrees_east")
  end if

  if (lonStrt.lt.0) then
      lon = (/ lon - 180. /)  ; subtract 180 from all values
      lon&lon = lon           ; update coordinates
  end if

                              ; read in land/sea mask basemap file
  b    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")

; lsm: five possible values: 0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf

  if (isvar("typeMask") .and. typeMask.eq."integer") then
      lsm  = localByte2Int( landsea_mask(b->LSMASK,lat,lon) ) ; type is "integer"
  else
      lsm  = landsea_mask(b->LSMASK,lat,lon)                  ; type is "byte"
  end if
  lsm@long_name = "Land-Sea Mask"
  lsm@units     = "0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf"
  printVarSummary(lsm)

;=================== OPTIONAL =============================
; Manually edit lsm [type byte] to create a custom land mask. [range 0-7 if byte]
;==========================================================
 ;lsm  = where (lsm.eq.2, 0, lsm) ; make lakes "ocean"
 ;lsm  = where (lsm.eq.2, 1, lsm) ; make lakes "land"
 ;lsm  = where (lsm.eq.3, 1, lsm) ; make island "land"
 ;lsm  = where (lsm.eq.4, 1, lsm) ; make ice "land"
 ;lsm@units     = "change as required"

;==========================================================
; Create netCDF
;==========================================================
   system("/bin/rm -f "+diro+filo)
   ncdf = addfile (diro+filo , "c")

   ncdf@title = "Land-Sea Mask"
   ncdf@creation_date = systemfunc("date")
   ncdf->LSMASK = lsm

;==========================================================
; Create netCDF
;==========================================================
if (PLOT) then

   wks = gsn_open_wks("ps","lsmask")
   gsn_define_colormap(wks,"wh-bl-gr-ye-re")

   res = True
   res@mpFillOn         = False                ; do not color-fill the map
   res@gsnTickMarksOn   = False                ; turn off all tick marks
   res@mpPerimOn        = True                 ; turn the map perimeter on
   res@mpPerimDrawOrder = "PostDraw"           ; draw the map perimeter last
   res@cnLinesOn        = False                ; turn off the contour lines
   res@cnLineLabelsOn   = False                ; turn off the contour line labels
   res@cnLevelSelectionMode = "ExplicitLevels" ; explicitly set the levels via cnLevels
   res@cnLevels         = (/1.,2.,3.,4./)      ; set the levels
   res@cnFillOn         = True                 ; turn on color fill
   res@cnFillMode       = "RasterFill"         ; use raster fill
   res@cnFillColors     = (/60,100,20,140,5/)  ; set the colors that will be used to color fill
   res@lbLabelStrings   = ispan(0,4,1)         ; labels for the labelbar boxes
   res@lbLabelAlignment = "BoxCenters"         ; put the labels in the center of the label bar boxes
   res@lbTitleString    = "0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf"  ; labelbar title
   res@lbTitleFontHeightF = 0.0125             ; labelbar title font height

   plot = gsn_csm_contour_map_ce(wks,lsm,res)
end if
end

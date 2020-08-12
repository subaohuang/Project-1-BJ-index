;;********************************************************
;;将模式的海温数据的大西洋地区的海温增加，Hadley SST
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;用来跑F2000实验
;;lnshuheng 2019.12.07 
;**********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
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

;;;;;;
undef("create_mask_array")
function create_mask_array(data,latlon_polygon:list)
local lat1d,lon1d
begin
  lat1d = ndtooned(conform(data,data&latitude,1))
  lon1d = ndtooned(conform(data,data&longitude,2))
  imask = reshape(gc_inout(lat1d,lon1d,latlon_polygon[0],latlon_polygon[1]),\
                  dimsizes(data))
  return(imask)
end
;----------------------------------------------------------------------
; Given a data array and a mask array of the same size, mask the
; data array and copy metadata.
;----------------------------------------------------------------------
undef("apply_mask")
function apply_mask(data:numeric,imask:logical)
begin
  if(.not.isatt(data,"_FillValue"))
    data@_FillValue = default_fillvalue(typeof(data))
  end if
  data_mask = where(imask,data,data@_FillValue)
  copy_VarMeta(data,data_mask)
  return(data_mask)
end


;;------------------------------------------------------------------
begin


;;;---------------读入Had变化的海温数据----------------------------
    fils1 = "/home/ys17-19/lsh/data/sst/HadISST_sst.nc"
    f1 = addfile(fils1, "r")
    time = f1->time
    date = cd_calendar(time, 1) 
    lon = f1->longitude
    lat = f1->latitude
    timeind = ind(date.le.201512.and.date.ge.197901)
    ; SST_Had  = f1->sst(timeind,:,:)
    SST_Had_r = lonFlip(f1->sst(timeind,::-1,:))   ;;;;0-360
    SST_Had_r(:,{35:},:) = SST_Had_r@_FillValue
    SST_Had_r(:,{:-35},:) = SST_Had_r@_FillValue

;;;;   
printVarSummary(SST_Had_r)


;;;;mask,只有留下太平洋三角形的溫度變化

    ; lat_polygon = (/0.0,-15,-30.0,30.0,15.0,0.0/)
    ; lon_polygon = (/165.0,227.5,290.0,255.0,210.0,165.0/)


    lat_polygon = (/0.0,-16.0,-30.0,30.0,15.0,0.0/)
    lon_polygon = (/165.0,226.5,290.0,255.0,210.0,165.0/)

    latlon_polygon = [/lat_polygon,lon_polygon/]
    printVarSummary(latlon_polygon)
    imask = create_mask_array(SST_Had_r,latlon_polygon)
    printVarSummary(imask)
    SST_Had_h = apply_mask(SST_Had_r,imask)

    print("do this")
;;;;






;;;差值到pop的坐标上
    TEMP_file = addfile ("/home/ys17-19/lsh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.TEMP.mon.0251-0280.nc", "r")
    TLONG = TEMP_file->TLONG
    TLAT = TEMP_file->TLAT
    TEMP = TEMP_file->TEMP
    TEMP@lat2d = TLAT
    TEMP@lon2d = TLONG
    printVarSummary(TEMP)
    

    SST_Had = PopLatLon(SST_Had_h,"1x1d","gx1v6","bilin","da","100716")

    SST_Had!1 = "nlat"
    SST_Had!2 = "nlon"
    SST_Had@lat2d   = TLAT
    SST_Had@lon2d   = TLONG
  
    printVarSummary(SST_Had)
;;;


;;;;计算每个月的趋势
  dim = dimsizes(SST_Had)
  SST_atlantic_rc1 = new((/12,dim(1),dim(2)/), typeof(SST_Had))



  SST_Had_season1 = month_to_season(SST_Had, "MAM")
  do i=0,11
    SST_atlantic_rc1(i,:,:) = regCoef_n(ispan(1,dimsizes(SST_Had_season1&time),1),SST_Had(i:dim(0)-1:12,:,:),0,0)
  end do 

  ; prob1 = rc_test(SST_atlantic_rc1)
  ; SST_atlantic_rc1  = where(prob1.le.0.05, SST_atlantic_rc1, SST_atlantic_rc1@_FillValue)
  printVarSummary(SST_atlantic_rc1)
  print(SST_atlantic_rc1@_FillValue)

  ssta_atlantic1 = SST_atlantic_rc1 * dimsizes(SST_Had_season1&time)
  ssta_atlantic1 = where(ismissing(ssta_atlantic1) ,0, ssta_atlantic1)                         
;;; 

;;;;;
  ssta = ssta_atlantic1
  ssta!0 = "month"
  ssta!1 = "nlat"
  ssta!2 = "nlon"
  ssta@lat2d   = TLAT
  ssta@lon2d   = TLONG

  copy_VarAtts(SST_Had, ssta)
;;;;

; plot
  plot = new(1, graphic)
  wks  = gsn_open_wks("pdf","SST_Had_r_mask")            ; send graphics to PNG file

  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False

  res@cnFillOn             = True                 ; turn on color
  res@cnLinesOn            = False                ; turn off contour lines

  res@gsnAddCyclic         = True                 ; force cyclic value
  res@mpFillOn             = False
  res@mpFillDrawOrder      = "PreDraw"           ; color of land
  ; res@mpLandFillColor      = "grey"               ; color of land

  res@lbLabelsOn           = True

  res@cnFillPalette        = "BlueWhiteOrangeRed"               ; set color map
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF       = -1.0             ; set min contour level
  res@cnMaxLevelValF       = 1.0              ; set max contour level      
  res@cnLevelSpacingF      = 0.1 
  ; res@cnMinLevelValF       = 22             ; set min contour level
  ; res@cnMaxLevelValF       =  36             ; set max contour level      
  ; res@cnLevelSpacingF      = 2 

  res@mpMinLatF            = -90.                ; range to zoom in on
  res@mpMaxLatF            =  90.
  res@mpMinLonF            = 0.
  res@mpMaxLonF            =  360.
  res@mpCenterLonF = 180.0

  ; WGT_365 = (/WGT_365*TEMPMC_DAY_FINAL/)

plot(0) = gsn_csm_contour_map(wks,ssta(1,:,:),res)

draw(plot)
frame(wks)
print("done draw")
exit()

;;;;;


  ;;;;
;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f /home/ys17-19/lsh/CESM-data/B/modify_sst/Had_sst_addatlanticwarming_"+"year.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/CESM-data/B/modify_sst/Had_sst_addatlanticwarming_"+"year.nc","c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
    ncdf->ssta   =  (/ssta/)                           ; 4D               
    ncdf->TLAT = TLAT
    ncdf->TLONG = TLONG
       

end 

     







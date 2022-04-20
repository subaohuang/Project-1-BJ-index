
;----------------------------------------------------------------------
undef("read_ocean_land_colormap")
function read_ocean_land_colormap(num_ocean)
local cmap_ocn, cmap_lnd
begin
  cmap_ocn = read_colormap_file("GMT_ocean")
  cmap_lnd = read_colormap_file("OceanLakeLandSnow")
  ; cmap_lnd = read_colormap_file("GMT_globe")
  ; newcmap = array_append_record(cmap_ocn(0:num_ocean-1,:),cmap_lnd(2::2,:),0)
  newcmap = array_append_record(cmap_ocn(0:num_ocean-1,:),cmap_lnd(120:240:4,:),0)

  ; newcmap = array_append_record(cmap_ocn(0:num_ocean-1,:),cmap_lnd(120:225:3,:),0)
  ; newcmap =cmap_lnd(136:241,:)
  
  return(newcmap)
end

;----------------------------------------------------------------------
; Given min/max elevation value and the number of colors you want
; for ocean, this function creates a custom color map AND calculates 
; levels for both land and ocean, based on an elevation pivot point.
;
; The assumption is that the end ocean index location is one
; less than the begining of the land index location.
;----------------------------------------------------------------------
undef("calc_levels_and_colors")
function calc_levels_and_colors(wks,emin,emax,split_elev,num_ocean_values)
local start_ocean, ocean_range, land_range, olevels, llevels, nol, nll, clen
begin
  cmap = read_ocean_land_colormap(num_ocean_values)
  clen = dimsizes(cmap(:,0))

  start_ocean = 0
  end_ocean   = num_ocean_values-1
  start_land  = end_ocean+1
  ocean_range = end_ocean-start_ocean+1
  land_range  = (clen-start_land+1)

  olevels     = fspan(emin,split_elev,ocean_range)
  ; llevels     = fspan(split_elev,emax,land_range)
  ; levels1=ispan(-65,535,60)
  ; levels2=ispan(560, 2025, 58)
  ; levels3=ispan(2525, 5025,5 )
  ; llevels =   array_append_record(array_append_record(levels1,levels2,0),levels3,0)
 ; llevels=(/-65,50,60,70,80,90,100,200,300,400,500,1000,1500,2000,3000,4000/)

 ; llevels=(/-65,0,10,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,1500,2000,3000,4000/)
   levels1=fspan(-65,0,5)
  levels2=fspan(10 ,100, 10)
  levels3=fspan(110, 1000,10)
  levels4=fspan(1100, 4000,10)
  ; levels5=fspan(3100, 4000,10)
  llevels = array_append_record(array_append_record(array_append_record(levels1,levels2,0),levels3,0),levels4,0)

  ; llevels=(/-65,-  0,10,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,1500,2000,3000,4000/)

  nol         = dimsizes(olevels)
  nll         = dimsizes(llevels)
  levels      = new((nol-1)+(nll-2),float)
  levels(0:nol-2) = olevels(1:)
  levels(nol-1:)  = llevels(1:nll-2)
  return([/levels,cmap/])
end

;----------------------------------------------------------------------
; This function reads a binary file containing elevation data and
; generates the necessary lat/lon coordinate arrays for plotting later.
; The information on the binary file is provided at the beginning of
; this script.
;
; The binary file was downloaded from:
;     http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO5/TOPO/ETOPO5/
;----------------------------------------------------------------------
undef("read_elev_data")
function read_elev_data(topo_file)
local nlat, nlon, topo_file, lat, lon
begin
;---Read data as a straight binary file
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev = cbinread(topo_file,(/nlat,nlon/),"short")

;---Create 1D coordinate arrays
  lat       = fspan(90,-90,nlat)
  lon       = fspan(0,360,nlon)
  lat!0     = "lat"
  lon!0     = "lon"
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  lat&lat   = lat
  lon&lon   = lon

;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat
  elev&lon  = lon
  elev@_FillValue=integertoshort( 9999999)

 printVarSummary(elev)

  elev=integertoshort(where(elev.le.0, elev@_FillValue, elev))
  ; printVarSummary(elev_m)
  return(elev)
end

;----------------------------------------------------------------------
; This procedure draws a global 5' topographic map by contouring the 
; given elevation data.
;----------------------------------------------------------------------
undef("draw_topo_map")
procedure draw_topo_map(wks,elev,title)
local res, labels, nlevels
begin

;;;读入海温数据
 
  ;;;读入海温数据
   diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201312.and.date.ge.197901)
  ; sst=lonFlip(f1->sst(timeind,:,:))
  sst=f1->sst(timeind,:,:)
  printVarSummary(sst)

  sst!0="time"
  sst!1="lat"
  sst!2="lon"



   ;******************************************************************
    ; use poisson_grid_fill to fill missing value
    ;******************************************************************
    nscan = 1500           ; usually *much* fewer
    eps=0.01
    gtype = False            ; "gendat" does not generate cyclic fields
    guess = 1                ; use zonal means
    relc  = 0.6              ; standard relaxation coef
    opt   = 0
    
    poisson_grid_fill(sst, gtype, guess, nscan, eps, relc, opt)
    
    fils2=systemfunc("ls /mnt/f/share/shixiong/landsea.nc")



  ssta=(/rmMonAnnCycTLL(sst(:,:,:))/)
  copy_VarCoords(sst(:,:,:), ssta)
 

 
  sst_djf=month_to_season(sst, "DJF")
  sst_s=(/sst_djf(1:,:,:)/)*1.0
  copy_VarCoords(sst_djf(0:dimsizes(sst_djf&time)-2,:,:), sst_s)
  delete(sst_s&time)
  sst_s&time=ispan(1979, 2012, 1)  
  ;;;;


  ssta_djf=month_to_season(ssta, "DJF")
  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
  copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
  delete(ssta_s&time)
  ssta_s&time=ispan(1979, 2012, 1)  
  ;;;;


  ;;;两个时间段的海温异常
  nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
  nino3_std=dim_standardize_n(nino3, 1, 0)
  
  ;;;挑选年份 canino  enino year
  caind=ind(nino3_std.ge.1)
  laind=ind(nino3_std.le.-1)
  midind=ind(nino3_std.le.1.and.nino3_std.ge.-1)
  
  ssta_e=dim_avg_n_Wrap(ssta_s(caind,:,:),0)
  ssta_l=dim_avg_n_Wrap(ssta_s(laind,:,:),0)
  ssta_n=dim_avg_n_Wrap(sst_s(midind,:,:),0)
;;;


  ;;;

;---Set some resources for contouring and mapping
  res                    = True

  res@gsnMaximize        = True     ; maximize plot in frame

  res@cnFillOn           = True     ; turn on contour fill
  res@cnFillMode         = "MeshFill"
  res@cnLinesOn          = False    ; turn off contour lines
  res@cnLineLabelsOn     = False    ; turn off line labels
  res@cnInfoLabelOn      = False    ; turn off info label
  res@lbBoxLinesOn       = False    ; turn off labelbar box lines

  res@lbLabelBarOn=False
  res@tmXBOn=False
  res@tmYLOn=False

  
;---Calculate "nice" contour levels, and create a color map to match
  ; split_elev               = -62; -68    ; meters
  split_elev               = -65; -68    ; meters
  num_ocean_colors         = 43
  ; levels_and_colors        = calc_levels_and_colors(wks,min(elev),max(elev),split_elev,num_ocean_colors)
  levels_and_colors        = calc_levels_and_colors(wks,min(elev),5000,split_elev,num_ocean_colors)
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels_and_colors[0]
  res@cnFillPalette        = levels_and_colors[1]
  res@gsnAddCyclic         = False           ; don't add longitude cyclic point

  res@mpFillOn                    = False
  res@pmTickMarkDisplayMode       = "Always"   ; Nicer map labels
  res@mpGeophysicalLineThicknessF = 2          ; Thicker map outlines
  ; res@mpGeophysicalLineColor="white"
  ; res@tiMainString                = title      ; Main title
 
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF                = min(elev&lat)
  res@mpMaxLatF                = max(elev&lat)
  res@mpMinLonF                = min(elev&lon)
  res@mpMaxLonF                = max(elev&lon)
  res@mpCenterLonF             = avg(elev&lon) 

  res@cnLineDrawOrder="PreDraw"
  res@cnFillDrawOrder="PreDraw"
  res@gsnDraw=False
  res@gsnFrame=False
  res@lbLabelBarOn=False


  res@tmXBLabelFontHeightF=0.05


    plot = gsn_csm_contour(wks,elev,res)
 ;;;******************************************************************************
 ;; 叠上去的海温
 ; colormap=read_colormap_file("matlab_jet")
 ; colormap=read_colormap_file("BlueYellowRed")




    ; ;;;厄尔尼诺
    ; cmap=new((/17,4/),typeof(colormap))
    ; cmap(0:4,:)=colormap(2:90:22,:)
    ; cmap(5:16,:)=colormap(125:246:11,:)
    ; ;;;

    ; ;;;lanina
    ; cmap=new((/13,4/),typeof(colormap))
    ; cmap(0:7,:)=colormap(2:93:13,:)
    ; cmap(8:12,:)=colormap(125:225:25,:)
    ; ;;;



  res2=True
  res2@gsnFrame=False
  res2@gsnDraw=False
  res2@cnFillOn=True
  res2@cnLinesOn=False
  res2@cnInfoLabelOn=False
  res2@cnLineLabelsOn=False

  res2@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res2@cnMinLevelValF       =  -1.5
  res2@cnMaxLevelValF       =  1.5
  res2@cnLevelSpacingF      = 0.1
  res2@mpLimitMode           = "LatLon"

  res2@mpMinLatF                = min(elev&lat)
  res2@mpMaxLatF                = max(elev&lat)
  res2@mpMinLonF                = min(elev&lon)
  res2@mpMaxLonF                = max(elev&lon)
  res2@mpCenterLonF              =avg(elev&lon) 

  res2@mpFillOn=True

  res2@mpLandFillColor="black"
  res2@mpGeophysicalLineColor ="black"
  res2@mpOutlineOn = True ; turn off outline 
  res2@mpDataSetName = "Earth..4"
  res2@mpDataBaseVersion = "MediumRes"
   res2@mpFillDrawOrder="PostDraw"
   res2@tmXBOn=False
   res2@tmYLOn=False
   res2@tmYROn=False
   res2@tmXTOn=False
  
  res2@cnMinLevelValF       =  -10
  res2@cnMaxLevelValF       =  32
  res2@cnLevelSpacingF      = 1

  res2@cnFillPalette ="matlab_jet"
  res2@lbLabelBarOn=False
  res2@vpHeightF=1
  res2@vpWidthF=1
  plot1=gsn_csm_contour_map(wks, ssta_n, res2)

 ; overlay(plot, plot1)
  ; overlay(plot1, plot)
  draw(plot1)
  frame(wks)



end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
  diri="/mnt/f/share/" 
     wks_type = "png"
   wks_type@wkWidth = 8000
   wks_type@wkHeight = 8000
  wks = gsn_open_wks(wks_type,diri+"地形加中性海温")          ; send graphics to PNG file
  topo_filename=diri+"ETOPO5.DAT"
  elev = read_elev_data(diri+"ETOPO5.DAT")
  elev = read_elev_data(topo_filename)
  draw_topo_map(wks,elev(:,:),topo_filename)
  printVarSummary(elev)


  




end

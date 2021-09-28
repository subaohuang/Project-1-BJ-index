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

  ; startmon = 198201
  ; endmon = 201912

  startmon = 197901
  endmon = 201512

  ; f1 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
  ; date:=cd_calendar(f1->time, 1)
  ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
  ; olr = short2flt(f1->olr(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

  ; printVarSummary(olr)


  f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")

  date:=cd_calendar(f1->time, 1)
  timeind :=ind(date.ge.startmon.and.date.le.endmon)
  olr = f1->precip(timeind,{lat1:lat2},{lon1:lon2}) ;;;2.5 ×2.5 grid


  f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
  date:=cd_calendar(f2->time, 1)

  timeind:=ind(date.le.201512.and.date.ge.197901)

  sst_r =  lonFlip(f2->sst(timeind,:,:))

  sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

  printVarSummary(sst)


  fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
  f3 = addfiles(fils3, "r")

  date:=cd_calendar(f3[:]->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)
  time = f3[:]->time(timeind)
  lat = f3[1]->latitude
  lon = f3[1]->longitude


  div_r = new((/dimsizes(time), dimsizes(lat),dimsizes(lon)/),float)


  do i = 0,dimsizes(timeind)/12-1

  div_r(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->d(:,{850},:,:))

  end do
 
  ; div = dim_avg_n_Wrap(div_r(:,:,{lat1-1:lat2+1},{lon1-1:lon2+1}),1)
  div = div_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})


  div = (/div*10^5/)
  div = (/-1*div/)

  ; div = (/rmMonAnnCycTLL(div)/)


;;;;chazhi 
   print("**********start interp************")
    ; xi = sst&lon
    ; yi = sst&lat

    xi = sst&longitude
    yi = sst&latitude

    xo = olr&lon
    yo = olr&lat



    sst_chazhi = linint2_Wrap(xi, yi, sst, True, xo, yo, 0)

    sst_chazhi!1= "lat"
    sst_chazhi!2= "lon"

    sst_chazhi&lat = yo
    sst_chazhi&lon = xo


    div_chazhi = linint2_Wrap(div&longitude, div&latitude, div, True, xo, yo, 0)

    div_chazhi!1= "lat"
    div_chazhi!2= "lon"

    div_chazhi&lat = yo
    div_chazhi&lon = xo
;;;;
   delete(sst)
   delete(div)
   
;;;;
    printVarSummary(sst_chazhi)
    olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)
;;;;

;;;;anomaly 

    olr = (/rmMonAnnCycTLL(olr)/)
    sst_chazhi = (/rmMonAnnCycTLL(sst_chazhi)/)
    div_chazhi = (/rmMonAnnCycTLL(div_chazhi)/)
;;;

;;;;Monthly bin

   sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
   olr = (/runave_n_Wrap(olr, 3, 0, 0)/)

   dim = dimsizes(sst_chazhi)
   
  sst_chazhi_june = sst_chazhi(5::12,:,:)
  sst_chazhi_july = sst_chazhi(6::12,:,:)
  sst_chazhi_aug = sst_chazhi(7::12,:,:)

  sst_chazhi_jj =  array_append_record(sst_chazhi_june,sst_chazhi_july,0)
  sst_chazhi_jja = array_append_record(sst_chazhi_jj , sst_chazhi_aug,0)


  olr_june = olr(5::12,:,:)
  olr_july = olr(6::12,:,:)
  olr_aug = olr(7::12,:,:)

  olr_jj =  array_append_record(olr_june,olr_july,0)
  olr_jja = array_append_record(olr_jj, olr_aug,0)


  ; div_june = div_chazhi(2::12,:,:)
  ; div_july = div_chazhi(3::12,:,:)
  ; div_aug = div_chazhi(4::12,:,:)

  ; div_jj =  array_append_record(div_june,div_july,0)
  ; div_jja = array_append_record(div_jj , div_aug,0)

   ;  sst_one = ndtooned(sst_chazhi_jja(:,{0:15},{125:160}))
   ;  olr_one = ndtooned(olr_jja(:,{0:15},{125:160}))

    ; sst_one = ndtooned(sst_chazhi(:,{0:20},{125:160}))
;;;;

;;;;计算区域平均

   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi_jja(:,{0:20},{125:160}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr_jja(:,{0:20},{125:160}), (/1,2/))


   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi_jja(:,{5:25},{125:150}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr_jja(:,{5:25},{125:10}), (/1,2/))


   sst_wnp = sst_chazhi(:,{0:20},{125:160})

   olr_wnp = olr(:,{0:20},{125:160})
   div_wnp = div_chazhi(:,{0:20},{125:160})

   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi(:,{5:25},{120:150}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr(:,{5:25},{120:150}), (/1,2/))
; ;;


;;;计算超前滞后回归
    
    x = sst_wnp  
    ; y = olr_wnp
    ; y = div_wnp
    y=olr_wnp

    mxlag    = 4


    x_Lead_y = esccr(sst_wnp(lat|:,lon|:,time|:),olr_wnp(lat|:,lon|:,time|:),mxlag)
    y_Lead_x = esccr(olr_wnp(lat|:,lon|:,time|:),sst_wnp(lat|:,lon|:,time|:),mxlag)    ; switch the order of the series
    
    
    dim = dimsizes(sst_wnp)
    ccr = new ( (/dim(1),dim(2),2*mxlag+1/), float)    
    ccr(:,:,0:mxlag-1) = y_Lead_x(:,:,1:mxlag:-1)  ; "negative lag", -1 reverses order
    ccr(:,:,mxlag:)    = x_Lead_y(:,:,0:mxlag)     ; "positive lag"
    
    ccr!0 = "lat"
    ccr!1 = "lon"
    ccr&lat  = sst_wnp&lat
    ccr&lon  = sst_wnp&lon
    ccr&lat@units = "degree_north"
    ccr&lon@units = "degree_east"

    plot = new(mxlag*2+1,graphic)
    wks=gsn_open_wks("pdf", "region-WNP-sst&div超前滞后相关")

    res = True
    res@gsnAddCyclic = False
    res@gsnFrame=False
    res@gsnDraw=False
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=1
    res@pmLabelBarOrthogonalPosF = 0.17
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = 0
    res@mpMaxLatF             = 20
    res@mpMinLonF             = 125
    res@mpMaxLonF             = 160
    res@lbOrientation="horizontal"
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpFillOn = False

    res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF      = -0.5        ; set min contour level
    res@cnMaxLevelValF      = 0.5        ; set max contour level
    res@cnLevelSpacingF     = 0.1       ; set contour spacing


    res@lbLabelFontHeightF = 0.02
    res@lbAutoManage = False

   do i=0, dimsizes(plot)-1
     plot(i) = gsn_csm_contour_map(wks, ccr(:,:,i), res)
   end do 

   pres = True
   gsn_panel(wks, plot, (/3,3/), pres)

 end 




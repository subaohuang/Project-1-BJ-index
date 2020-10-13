begin
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201312.and.date.ge.197901)
  sst=f1->sst(timeind,:,:)
  ; sst=f1->sst(timeind,:,:)
  
  sst!0="time"
  sst!1="lat"
  sst!2="lon"
  delete(date)
  delete(f1)
  
  ssta=(/rmMonAnnCycTLL(sst(:,:,:))/)
  copy_VarCoords(sst(:,:,:), ssta)

  ssta_djf=month_to_season(ssta, "DJF")
  ssta_son=month_to_season(ssta, "SON")

  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    ; ssta_s=(/ssta_son(0:33,:,:)/)
    copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
   delete(ssta_s&time)
   ssta_s&time=ispan(1979, 2012, 1)  
 ;;;两个时间段的海温异常
 ssta_s1=dim_avg_n_Wrap(ssta_s({1979:1994},:,:),0)
 ssta_s2=dim_avg_n_Wrap(ssta_s({1995:},:,:),0)
 ;;;

;;;处理成冬季平均场
  sst_djf=month_to_season(sst, "DJF")
  sst_son=month_to_season(sst, "SON")

  sst_s=(/sst_djf(1:,:,:)/)*1.0
    ; sst_s=(/sst_son(0:33,:,:)/)
  copy_VarCoords(sst_djf(0:dimsizes(sst_djf&time)-2,:,:), sst_s)
   delete(sst_s&time)
   sst_s&time=ispan(1979, 2012, 1)  
 ;;;两个时间段的海温异常
 sst_s1=dim_avg_n_Wrap(sst_s({1979:1994},:,:),0)
 sst_s2=dim_avg_n_Wrap(sst_s({1995:},:,:),0)
 sst_s3=sst_s2 - sst_s1
 copy_VarCoords(sst_s1, sst_s3)

 ;;;


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/ 2个时间段的海温_整个地球_mid120")


  res                   = True
  ; res@mpProjection      = "Mollweide"       ; choose projection
  res@mpProjection      = "Orthographic"       ; choose projection
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
  res@mpCenterLonF      = 120.             ; choose center lon
  res@mpCenterLatF      = 10.               ; choose center lat





  res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpGridLatSpacingF = 30.               ; spacing for lat lines
  res@mpGridLonSpacingF = 30.               ; spacing for lon lines
  res@mpFillOn          = False


  res@cnFillOn          = True              ; color plot desired
  res@cnFillPalette     = "NCV_blu_red"
  res@cnLineLabelsOn    = False             ; turn off contour lines
  res@txFontHeightF     = 0.015 
  res@cnLinesOn         =False

  res@vpXF            = 0.1                 ; make plot bigger
  res@vpYF            = 0.9         
  res@vpWidthF        = 0.8
  res@vpHeightF       = 0.8
  
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =  -1
  res@cnMaxLevelValF       =  1
  res@cnLevelSpacingF      =  0.1

  res@lbLabelFontHeightF  = 0.015           ; label bar font height
  res@lbLabelBarOn =True


  contour = gsn_csm_contour_map(wks,sst_s3,res)  ; create the plot

  
end




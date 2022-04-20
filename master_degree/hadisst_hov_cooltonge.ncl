;hadisst_eof分析
begin
  diri="/mnt/f/share/ersst/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  yearstart=195801
  yearend=201212
  timeind=ind(date.le.yearend.and.date.ge.yearstart)

  sst=lonFlip(f1->sst(timeind,{10:-10},:))

  
  sst_annual=month_to_annual(sst(:,:,{240:270}), 1)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  printVarSummary(sst_annual)
  
  sst_hov=dim_avg_n_Wrap(sst_annual, 1)
  
  sst_sm=runave_n_Wrap(sst_hov, 5, 1, 0)

  ;;;画图

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/HAD_sst冷舌的变化")
  gsn_define_colormap(wks, "MPL_PuBu")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  

  res@cnLineColors="Transparent"  ;
  res@lbLabelBarOn=True
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  res@cnFillMode="RasterFill"
  res@cnRasterSmoothingOn=True

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

  ; res@cnMinLevelValF       =   24.0
  ; res@cnMaxLevelValF       =   28.0
  ; res@cnLevelSpacingF      =  0.25

  res@tmXTOn=False
  res@tmYROn=False
  
  res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  res@tmYLLabelsOn=True
  res@tmYUseLeft=True

  res@tmYLValues    = ispan(3,dimsizes(timeind)/12-1,5)    
  res@tmYLLabels    = ispan(1960,2010,5)

    ; res@tmYLValues    = ispan(2,dimsizes(timeind)/12-1,5)
    ; res@tmYLLabels=ispan(1980, 2010, 5)

  res@tmXBValues=ispan(240,270,5)
  res@tmXBLabels=ispan(240,270,5)

  res@tmXBMinorOn=True
  res@tmXBMinorValues=ispan(240,270,5)
  res@tmXBMinorOn=True
  
  res@cnInfoLabelOn=False
  res@vpWidthF = 0.8
  res@vpHeightF = 0.2
 
  plot=new(1,graphic)
  plot(0)= gsn_csm_contour(wks,sst_sm,res)
  
  ;;;28.75作为warm pool 
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = 26
    res@cnLineThicknessF = 4.
    res@cnFillOn        = False

    plot1 = gsn_csm_contour(wks, sst_sm, res)
    
    overlay(plot, plot1)

    draw(plot)
    frame(wks)
    end
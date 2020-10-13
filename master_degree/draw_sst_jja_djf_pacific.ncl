;************************************************************
;annual cycle SST
;************************************************************
begin
  ;;;read sst
    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    time = f2->time
    date = cd_calendar(time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    sst=lonFlip(f2->sst(timeind,{60:-60},:))
  ;;;
  
  
    sst_jja = month_to_season(sst, "JJA")
    sst_djf = month_to_season(sst, "DJF")
    sst_jja_clm = dim_avg_n_Wrap(sst_jja,0)
    sst_djf_clm = dim_avg_n_Wrap(sst_djf,0)

    sst_clm = dim_avg_n_Wrap(sst, 0)

    sst_jja_a = sst_jja_clm - sst_clm

    sst_djf_a = sst_djf_clm - sst_clm

    copy_VarCoords(sst_clm, sst_jja_a)
    copy_VarCoords(sst_clm, sst_djf_a)
;************************************************************
;draw plot
;***********************************************************
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/SST的冬夏季异常场")

    gsn_define_colormap(wks,"hotcolr_19lev")
  ; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
    cmap=read_colormap_file("MPL_PuOr")
    cmap=(/cmap(::-1,0:3)/)

    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn =False

    res@cnFillPalette="GMT_polar"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnFillDrawOrder="PreDraw"
    res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.1
    
    res@tmXTOn = False
    res@tmYROn = False
    res@tmXBLabelFontHeightF=0.012
    res@tmYLLabelFontHeightF=0.012

    res@tmXBMajorOutwardLengthF=0
    res@tmXBMajorLengthF=0.007


    res@tmYLMajorOutwardLengthF=0
    res@tmYLMajorLengthF=0.01

    res@tmYLMajorThicknessF=2.0
    res@tmXBMajorThicknessF=2.0
    res@tmBorderThicknessF=2.5

 

    res@mpFillDrawOrder="PostDraw"
    res@mpFillOn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -10
    res@mpMaxLatF             = 10
    res@mpMinLonF             = 0
    res@mpMaxLonF             = 360
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@mpFillDrawOrder  = "PostDraw"

    ;---draw plot 
    plot=new(2,graphic)
    res@gsnLeftString="(a)JJA"
    plot(0) = gsn_csm_contour_map(wks,sst_jja_a,res)
    res@gsnLeftString="(b)DJF"
    plot(1) = gsn_csm_contour_map(wks,sst_djf_a,res)
   
    
    ;;;panel
    pres=True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar=True
    pres@lbOrientation ="Vertical"
    ; pres@pmLabelBarWidthF=0.5
    gsn_panel(wks, plot, (/2,1/), pres)
    frame(wks)
    end
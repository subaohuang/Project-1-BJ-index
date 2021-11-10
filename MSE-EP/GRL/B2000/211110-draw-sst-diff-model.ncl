;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"


function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 


function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

begin


  ; fils1 = "/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.TEMP.mon.0250-0299.nc"  
  fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h1.0251-0280.nc"
  fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_db_year_80_20/lsh_B2000_alt_north_db_year_80_20.TEMP.h1.0251-0280.nc"
 
  modelname = "lsh_B2000_alt_north_db_year_80_20"


    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    
    time1 = f1->time
    time2 = f2->time
    time1 = (/time1 -1/)
    time2 = (/time2 -1/)
    date1=cd_calendar(time1, 1)
    date2=cd_calendar(time2, 1)


    timeind1=ind(date1.ge.25601.and.date1.le.28012)
    timeind2=ind(date2.ge.25601.and.date2.le.28012)


    PSL_CT  =  f1->TEMP(timeind1,0,:,:)

  

    PSL_indian  =  f2->TEMP(timeind2,0,:,:)
    TLAT = f2->TLAT
    TLONG = f2->TLONG
   
    PSL_indian!1 = "nlat"
    PSL_indian!2 = "nlon"
    PSL_indian@lat2d   = TLAT
    PSL_indian@lon2d   = TLONG


    printVarSummary(PSL_indian)

    
    season = (/"MAM"/)
    wks=gsn_open_wks("pdf", modelname+"diff-Potensial-temperature-fourseason-0280.pdf")

    plot = new(1,graphic)
    plot1 = plot
    iplot = plot
    ; str3 = ispan(1, 4, 1) +""
    ; plot1 =new(4,graphic)
    ; iplot = new(4,graphic)

     ; order  = (/"(a)","(b)","(c)","(d)"/)
    
     ; plot_add = new(4,graphic)
    
    PSL_CT_split12 = split_12(PSL_CT)
    PSL_indian_split12 = split_12(PSL_indian)
    printVarSummary(PSL_CT_split12)


    PSL_CT_clm = dim_avg_n_Wrap(month_to_season(PSL_CT, season), 0)
    PSL_indian_clm = dim_avg_n_Wrap(month_to_season(PSL_indian, season), 0)
    

        PSL_CT_clm = dim_avg_n_Wrap(dim_avg_n_Wrap(PSL_CT_split12(0:4,:,:,:), 0),0)
    PSL_indian_clm = dim_avg_n_Wrap(dim_avg_n_Wrap(PSL_indian_split12(0:4,:,:,:), 0),0)
   
    
    PSL_diff_indian = PSL_indian_clm
    PSL_diff_indian = (/PSL_indian_clm - PSL_CT_clm/)

    

    ;;;t_test 
    prob_slp_indian = t_test(PSL_indian,PSL_CT)
   
    prob_slp_indian!0 = "nlat"
    prob_slp_indian!1 = "nlon"
    prob_slp_indian@lat2d   = TLAT
    prob_slp_indian@lon2d   = TLONG


    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    ; res@cnFillMode  = 
    ; res@cnRasterSmoothingOn=True
    res@trGridType = "TriangularMesh"

    ; res@cnFillPalette = "BlueDarkRed18"


    gsn_define_colormap(wks, "BlueWhiteOrangeRed")

    ; gsn_define_colormap(wks, "BlRe")
    ; gsn_define_colormap(wks, "BlueDarkRed18")

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -1.2
    res@cnMaxLevelValF       =   1.2
    res@cnLevelSpacingF      =  0.2


    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    res@mpMinLatF             = -50.0
    res@mpMaxLatF             = 50.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF = 180.0

    ; res@tmXBMode="Explicit"
    ; ; res@tmYLMode="Explicit"

    ; res@tmXBValues = (/-60,0,60,120,180,240/) 
    ; res@tmXBLabels = (/"60W","0","60E","120E","180","120W"/) 


    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3

    ; res@tmXBLabelFontHeightF=0.02
    ; res@tmYLLabelFontHeightF=0.02
    ; res@tmXBMajorThicknessF=1.0
    ; res@tmYLMajorThicknessF=1.0
    ; res@tmBorderThicknessF=1.0
    res@gsnLeftStringFontHeightF = 0.02

    res@pmLabelBarOrthogonalPosF = 0.2   ;;;上下位置

    
    res@gsnRightString = ""
    res@gsnLeftString = "JFMAM SST diff (ATL warming - CTRL)"

 
    plot = gsn_csm_contour_map(wks, PSL_diff_indian, res)

    ; res2 = res
    ; res2@cnLevelSelectionMode = "ExplicitLevels"
    ; res2@cnFillOn = False
    ; res2@cnLinesOn = True
    ; res2@cnLevels = (/0.5/)
    ; res2@cnLineColors = "black"
    ; res2@cnMonoLineColor  = True
    ; res2@cnLineThicknessF = 2
    ; res2@cnMonoLineThickness  = True
    ; res2@gsnLeftString = ""
    ; res2@cnInfoLabelOn = False

  

    ; plot1 = gsn_csm_contour(wks, PSL_diff_indian, res2)
   
   
    ; overlay(plot, plot1)
 
     

    ;;;画box
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   
    lat1 = 0
    lat2 = 20
    lon1 = 250
    lon2 = 280

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    
    plot1 = gsn_add_polyline(wks, plot, pgx_pool, pgy_pool, pgres) 



    ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1

    plot_add = gsn_add_polyline(wks, plot, pgx, pgy, pgres)


     ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
   
    res_pro@gsnCenterString = ""


    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    

    iplot=gsn_csm_contour(wks,prob_slp_indian,res_pro)
    overlay(plot, iplot)
draw(plot)
frame(wks)
    ; pres = True
    ; pres@gsnPanelLabelBar  = True
    ; pres@gsnMaximize = True
    ; pres@lbLabelStride = 1
    ; pres@lbOrientation = "Horizontal"
    ; ; pres@lbOrientation = "Vertical"
    ; ; pres@pmLabelBarHeightF = 0.15
    ; pres@pmLabelBarWidthF = 0.45
    ; pres@pmLabelBarHeightF = 0.04
    ; pres@pmLabelBarOrthogonalPosF = -0.008   ;;;上下位置
    ; pres@pmLabelBarParallelPosF = 0.02    ;;;左右位置 
    ; pres@lbLabelFontHeightF = 0.01
    ; pres@lbAutoManage =False

    ; gsn_panel(wks, plot, (/2,2/), pres)
    ; print("done script")
  end
  


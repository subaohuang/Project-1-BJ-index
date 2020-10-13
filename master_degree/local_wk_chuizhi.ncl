begin
    diri="/mnt/f/share/ncep/ncep1/"
    diri_pi="/mnt/f/share/era/data/"
    fils1=systemfunc("ls "+diri+"ncep1_m_zonal_mer.nc")
    f1=addfile(fils1, "r")
    m_zonal=f1->m_zonal(:,:,{-5:5},:)
    m_meri=f1->m_meri(:,:,{-5:5},:)
    printVarSummary(m_zonal)
    
    fi_zonal=f1->fi_zonal(:,:,{-5:5},:)
    fi_meri=f1->fi_meri(:,:,{-5:5},:)
    
 
    m_zonal_mean=dim_avg_n_Wrap(m_zonal(time|:,lat|:,level|:,lon|:), (/0,1/))
    m_meri_mean=dim_avg_n_Wrap(m_meri(time|:,lat|:,level|:,lon|:), (/0,1/))
    pi=get_pi("float")

    fi=cos((m_zonal&lat/180)*pi)
   
    w_zonal=m_zonal
    do i=0,dimsizes(m_zonal&lat)-1
    w_zonal(:,:,i,:)=m_zonal(:,:,i,:)*9.8*50/fi(i)
    end do
    w_zonal_mean=dim_avg_n_Wrap(w_zonal(time|:,lat|:,level|:,lon|:), (/0,1/))
   
    level=fi_zonal&level*100

    ud=-1*center_finite_diff_n (fi_zonal, level, False, 0, 1)
    vd=-1*center_finite_diff_n (fi_meri, level, False, 0, 1)

   copy_VarCoords(fi_zonal, ud)
   copy_VarCoords(fi_zonal, vd)

    ud_zonal_mean=dim_avg_n_Wrap(ud(time|:,lat|:,level|:,lon|:), (/0,1/))
    vd_zonal_mean=dim_avg_n_Wrap(vd(time|:,lat|:,level|:,lon|:), (/0,1/))
    
   m_zonal_lon=dim_avg_n_Wrap(m_zonal, 2)

    rc=regCoef_n(ispan(1,dimsizes(m_zonal&time),1), m_zonal_lon, 0, 0)
    rc=(/rc*120/)
    copy_VarCoords(m_zonal(1,:,1,:), rc)
     
   wks2=gsn_open_wks("pdf", ""+diri_pi+"local沃克经向(-5-5)平均和趋势")

;;;画图
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False

    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -0.004
    res2@cnMaxLevelValF = 0.004
    res2@cnLevelSpacingF = 0.0005
    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    res2@cnFillPalette="temp_19lev"
    ; res2@cnMonoFillScale=True
    res2@cnInfoLabelOn=False

   
    res2@lbOrientation="Horizontal"
    
    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    res2@tmXBValues=(/0,90,180,270,357.5/)
    res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(0,360,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7


    res2@vcRefMagnitudeF = 2.0                ; define vector ref mag
    res2@vcRefLengthF    = 0.04        ; define length of vec ref
    res2@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res2@vcMinDistanceF  = 0.015       ; thin out vectors
    res2@vcMapDirection  = False
   

        ; res@vcFillArrowHeadXF=0.1
    res2@vcFillArrowMinFracWidthF=1.0
    res2@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    ; res2@vcFillArrowWidthF=0.05
    ; res2@vcFillArrowHeadYF=0.15




    plot=new(2,graphic)



    plot(0)=gsn_csm_pres_hgt_vector(wks2,m_zonal_mean,ud_zonal_mean,w_zonal_mean, res2)


    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -0.0004
    res2@cnMaxLevelValF = 0.0004
    res2@cnLevelSpacingF = 0.00005
    plot(1)=gsn_csm_pres_hgt(wks2,rc,res2)
  pres=True
  gsn_panel(wks2, plot, (/2,1/), pres)


;;;;
 wks=gsn_open_wks("pdf", ""+diri_pi+"local沃克w")
  res = True
  res@gsnAddCyclic=False
  res@gsnDraw=False
  res@gsnFrame=False
  res@cnFillPalette="BlueRed"
  ; res@cnFillPalette="precip_diff_12lev"
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False
 
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  ; res@cnMinLevelValF       =   -0.005
  ; res@cnMaxLevelValF       =   0.005
  ; res@cnLevelSpacingF      =  0.001

  res@lbLabelBarOn=True
  res@lbOrientation   ="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -60
  res@mpMaxLatF             = 60
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 357.5
  res@mpCenterLonF          =180.0
  
  res@tmXTOn=False
  res@tmYROn=False
  res@tmYLMode="Explicit"
  res@tmYLValues=ispan(-60,60,30)
  res@tmYLMinorValues=ispan(-60,60,10)
  res@tmYLMinorOn=True
  res@tmYLLabels=(/"60~S~o~N~S","30~S~o~N~S","0","30~S~o~N~N","60~S~o~N~N"/)
  res@tmYLLabelFontHeightF=0.01

  res@tmXBMode="Explicit"
  res@tmXBValues=fspan(0.5,358,13)
  res@tmXBMinorOn=False  
  res@tmXBLabelFontHeightF=0.01
  res@tmXBLabelFontThicknessF=0.02
  res@tmXBLabels=(/"0","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
            "180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0"/)

  ; plot_w=gsn_csm_contour_map(wks, w_zonal_mean({-60:60},:), res)
 
  
  ; draw(plot_w)
  ; frame(wks)


   end

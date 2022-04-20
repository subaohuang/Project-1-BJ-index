function one_test(pv_ano,east_time)
begin


    pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)

    t_pv_ano_east = pv_east/var_pv_ano_east

    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)

    confi_pv_east_se=confi_pv_east(::4,::4)

 return confi_pv_east
 end


begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_global_u_v_1.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)

fils1=systemfunc("ls "+diri+"velocity_potensial.nc")
f2=addfile(fils1, "r")
pv=f2->pv

pv_anu=month_to_annual(pv(:,{20:-20},{30:290}), 1)

printVarSummary(pv_anu)

u_r=f1->u(tind,{850},{20:-20},{30:290})
v_r=f1->v(tind,{850},{20:-20},{30:290})

u=u_r*u_r@scale_factor+u_r@add_offset
v=v_r*v_r@scale_factor+v_r@add_offset
copy_VarCoords(u_r, u)
copy_VarCoords(v_r, v)
delete(u_r)
delete(v_r)

u_anu = month_to_annual(u,1)
v_anu = month_to_annual(v,1)

pv_ano  =  dim_rmvmean_n_Wrap(pv_anu, 0)
u_ano  =  dim_rmvmean_n_Wrap(u_anu, 0)
v_ano  =  dim_rmvmean_n_Wrap(v_anu, 0)

; u_ano=rmMonAnnCycTLL(f1->u(tind,{850},{20:-20},{30:290}))
; v_ano=rmMonAnnCycTLL(f1->v(tind,{850},{20:-20},{30:290}))

printVarSummary(u_ano)

diri3=diri+"locationxi.ascii"
diri4=diri+"locationdong.ascii"
locationxi=asciiread(diri3, 432, "float") 
locationdong=asciiread(diri4, 432, "float")

locationxi_anu=month_to_annual(locationxi, 1)
locationxi_anu_std=dim_standardize(locationxi_anu, 1)
locationxi_anu_std!0="time"
printVarSummary(locationxi_anu_std)
locationxi_anu_std&time=ispan(1,36,1)
;;;提取数据
   
  east_time = ind(locationxi_anu_std.gt.0.5)
  west_time=ind(locationxi_anu_std.lt.0.5)

 
  pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
  u_east=dim_avg_n_Wrap(  u_ano(east_time,:,:),0)
  v_east=dim_avg_n_Wrap( v_ano(east_time,:,:),0)

  pv_west=dim_avg_n_Wrap( pv_ano(west_time,:,:),0)
  u_west=dim_avg_n_Wrap( u_ano(west_time,:,:),0)
  v_west=dim_avg_n_Wrap( v_ano(west_time,:,:),0)


  pv_clm=dim_avg_n_Wrap( pv_ano(:,:,:),0)
  u_clm=dim_avg_n_Wrap( u_anu(:,:,:),0)
  v_clm=dim_avg_n_Wrap( v_anu(:,:,:),0)

  pv_clm2=dim_avg_n_Wrap(pv_anu, 0)

;;;计算t值
    siglvl=0.05

    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)

    t_pv_ano_east = pv_east/var_pv_ano_east

    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)

    confi_pv_east_se=confi_pv_east(::4,::4)


    prob1d=ndtooned(confi_pv_east_se)
    i1d=ind(prob1d.lt.0.05)
    i=ind_resolve(i1d, dimsizes(confi_pv_east_se))

    ilat=confi_pv_east_se&lat(i(:,0))
    ilon=confi_pv_east_se&lon(i(:,1))

;;------
  
    var_pv_ano_west    =  dim_variance_n(pv_ano(west_time,:,:), 0)    ; 13.477
    var_pv_ano_west  =  sqrt(var_pv_ano_west/n) 
    var_pv_ano_west  =where(var_pv_ano_west.eq.0, var_pv_ano_west@_FillValue, var_pv_ano_west)

    t_pv_ano_west = pv_west/var_pv_ano_west

    confi_pv_west=student_t(t_pv_ano_west , n-1)

    copy_VarCoords(pv_west, confi_pv_west)

    confi_pv_west_se=confi_pv_west(::4,::4)


    prob2d=ndtooned(confi_pv_west_se)
    i21d=ind(prob2d.lt.0.05)
    i2=ind_resolve(i21d, dimsizes(confi_pv_west_se))

    ilat2=confi_pv_west_se&lat(i2(:,0))
    ilon2=confi_pv_west_se&lon(i2(:,1))
;;;------
;风场的T检验
   confi_u_east= one_test(u_ano,east_time)
   ; confi_v_east= one_test(v_ano,east_time)
   

   u_east=where(confi_u_east.le.0.2, u_east, u_east@_FillValue)
   ; v_east=where(confi_v_east.le.0.2, v_east, v_east@_FillValue)

   confi_u_west= one_test(u_ano,west_time)
   ; confi_v_west= one_test(v_ano,west_time)

   u_west=where(confi_u_west.le.0.2, u_west, u_west@_FillValue)
   ; v_west=where(confi_v_west.le.0.2, v_west, v_west@_FillValue)

  
;---Create plot
  wks = gsn_open_wks("pdf",diri+"速度势场")      ; open a workstation

  
    res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnScalarContour        = True
    res@gsnLeftString=""
    res@gsnRightString="East Phase"
    ; gsn_define_colormap(wks,"BlWhRe")
    gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  
    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn = False

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnMinLevelValF       =  cnMinLevelValF
    ; res@cnMaxLevelValF       =   cnMaxLevelValF
    ; res@cnLevelSpacingF      =  cnLevelSpacingF
    res@cnFillDrawOrder="PreDraw"

    res@lbOrientation = "Vertical"
    res@lbBoxEndCapStyle="TriangleBothEnds"

    res@pmLabelBarWidthF=0.1     ;colorbar的宽度

    res@tmYLMode="Explicit"
    res@tmYLValues=ispan(-20,20,10)
    res@tmYLLabels=(/"20~S~o~N~S","10~S~o~N~S","0~S~o~N","10~S~o~N~N","20~S~o~N~N"/)
    res@tmYLLabelFontHeightF=0.01

    res@tmYROn=False
    res@tmXTOn=False
    

    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(30,290,30)
    res@tmXBLabelFontHeightF=0.01
    res@tmXBLabelFontThicknessF=0.02
    res@tmXBMinorOn=True
    res@tmXBLabels=(/"30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
                    "180~S~o~N~W","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W"/)
    


    res@mpDataSetName = "Earth..4"
    res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"white","transparent"/) ;填充位置相应的具体颜色

    res@mpOutlineOn = True; turn off outline 
    ; res@mpGeophysicalLineThicknessF= 4.          ; double the thickness of geophysical boundaries 
    res@mpNationalLineThicknessF= 2.               ; double the thickness of national boundaries 
    res@mpOutlineSpecifiers = (/"Land"/)
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(pv_anu&lat)
    res@mpMaxLatF             = max(pv_anu&lat)
    res@mpMinLonF             = min(pv_anu&lon)
    res@mpMaxLonF             = max(pv_anu&lon)
    res@mpCenterLonF          =avg(pv_anu&lon)
;
     
; Usually, when data is placed onto a map, it is TRANSFORMED to the
; specified projection. Since this model is already on a native lambert
; conformal grid, we want to turn OFF the tranformation.
;
  ; res@tfDoNDCOverlay           = True

   res@vcRefMagnitudeF          = 1.5           ; define vector ref mag
   res@vcRefLengthF             = 0.04          ; define length of vec ref
   res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
   res@vcMinDistanceF           = 0.015           ; thin vectors
   res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角
   
    ; res@vcRefMagnitudeF = 1.5                ; define vector ref mag
    ; res@vcRefLengthF    = 0.04        ; define length of vec ref
    ; res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    ; res@vcMinDistanceF  = 0.015       ; thin out vectors
    ; res@vcMapDirection  = False

  
    ; ; res@vcFillArrowHeadXF=0.1
    ; res@vcFillArrowMinFracWidthF=0.25
    ; res@vcFillArrowHeadMinFracYF=1.0
 
    ; res@vcFillArrowWidthF=0.05
    ; res@vcFillArrowHeadYF=0.1

 



   res@lbLabelStrings=ispan(-10, 10, 2)


 
 plot=new(3,graphic)
  plot(0) = gsn_csm_vector_scalar_map(wks,u_east,v_east,pv_east,res)
  delete(res@lbLabelStrings)
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =  -4000000
    res@cnMaxLevelValF       =   6000000
    res@cnLevelSpacingF      =  1000000
      res@lbLabelStrings=ispan(-4, 6, 1)
     res@vcRefMagnitudeF=8.0
        res@gsnRightString="Climate mean state"
  plot(1) = gsn_csm_vector_scalar_map(wks,u_clm,v_clm,pv_clm2,res)
  delete(res@lbLabelStrings)
  delete(res@cnMinLevelValF)
  delete(res@cnMaxLevelValF)
  delete(res@cnLevelSelectionMode)
  delete(res@cnLevelSpacingF)
   res@vcRefMagnitudeF=1.5
   res@lbLabelStrings=ispan(-10, 10, 2)
          res@gsnRightString="West Phase"
  plot(2) = gsn_csm_vector_scalar_map(wks,u_west,v_west,pv_west,res)

 

;;;;做t检验

   res_pro=True
   res_pro@gsnDraw=False
   res_pro@gsnFrame=False
   res_pro@cnLevelSelectionMode = "ExplicitLevels"
   res_pro@cnLinesOn=True
   
   res_pro@cnInfoLabelOn=False
   res_pro@cnLevels =0.05
   res_pro@cnLineThicknessF=4.0
   res_pro@cnLineColor="white"
   ; res_pro@cnLineLabelsOn=False

 ; res_pro@cnLineLabelBackgroundColor="Transparent"
 ; res_pro@cnLineThicknessF = 3.
   ; res_pro@cnFillOn        = True
   iplot1=gsn_csm_contour(wks, confi_pv_east , res_pro)
   iplot2=gsn_csm_contour(wks, confi_pv_west , res_pro)


  
   polyres                   = True
   polyres@gsMarkerIndex     = 16          ; polymarker style
   polyres@gsMarkerSizeF     = 3.         ; polymarker size
   polyres@gsMarkerColor   ="white"
    dum1=gsn_add_polymarker(wks, plot(0), ilon(:), ilat(:), polyres)
    dum2=gsn_add_polymarker(wks, plot(2), ilon2(:), ilat2(:), polyres)


   overlay(plot(0), iplot1)
   overlay(plot(2), iplot2)


  pres=True
  gsn_panel(wks, plot, (/3,1/), pres)







end




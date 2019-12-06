begin
; fils1  =  "/home/ys17-19/lsh/data/era/era_pre_global.nc"
; fils1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/uwnd.mon.mean.nc"
; fils1  =  "/home/ys17-19/lsh/data/era/era_sur_global.nc"
    fils1  =  "/home/ys17-19/data-observation/20thCRV2/uwnd.mon.mean.nc"
    fils2  =  "/home/ys17-19/data-observation/20thCRV2/vwnd.mon.mean.nc"
; 

    f1 = addfile(fils1, "r")
    f2 = addfile(fils2, "r")
    

    date  =  cd_calendar(f1->time, 1)

    timeind  =  ind(date.le.201212.and.date.ge.194001)
   

    u_heat = f1->uwnd(timeind,{850},:,:)
    v_heat = f2->vwnd(timeind,{850},:,:)

    u_heat = (/dtrend_msg_n(ispan(1,dimsizes(u_heat&time),1),u_heat,False,True,0)/)
    v_heat = (/dtrend_msg_n(ispan(1,dimsizes(v_heat&time),1),v_heat,False,True,0)/)
    ; v_heat = short2flt(f1->vwnd(timeind,{1000},:,:))
    u_heat&time = date(timeind)
    v_heat&time = date(timeind)
    delete(timeind)
    delete(date)
   
    ;;;
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    u_djf = month_to_season(u_heat, "DJF")
    v_djf = month_to_season(v_heat, "DJF")
    printVarSummary(u_djf)
    

    u_djf&time=ispan(1940,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    v_djf&time=ispan(1940,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    u_djf(0:dimsizes(u_djf&time)-2,:,:)=(/u_djf(1:,:,:)/)
    u_djf(dimsizes(u_djf&time)-1,:,:)=(/u_heat({201212},:,:)/)
    v_djf(0:dimsizes(u_djf&time)-2,:,:)=(/v_djf(1:,:,:)/)
    v_djf(dimsizes(u_djf&time)-1,:,:)=(/v_heat({201212},:,:)/)
    ;;;
     
    ;;;读入nino3.4指数
      f3=addfile("/home/ys17-19/lsh/data/sst/HadISST_sst.nc", "r")
      date=cd_calendar(f3->time, 1)
      timeind=ind(date.ge.194001.and.date.le.201212)
      sst=lonFlip(f3->sst(timeind,{20:-20},:))
      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)
      ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
      copy_VarCoords(ssta, ssta_de)
      ssta_de&time := date(timeind)
      nino3_4=dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))   
  
      nino3_4_djf=(/month_to_season(nino3_4, "DJF")/)

      nino3_4_djf!0 = "time"
      nino3_4_djf&time=ispan(1940,2012,1)  
      nino3_4&time = date(timeind)

      nino3_4_djf(0:dimsizes(nino3_4_djf&time)-2)=(/nino3_4_djf(1:)/)
      nino3_4_djf(dimsizes(nino3_4_djf&time)-1)=(/nino3_4({201212})/)

      nino3_4_anu_std=dim_standardize(nino3_4_djf, 1)
      nino3_4_anu_std!0="time" 
      nino3_4_anu_std&time=ispan(1940, 2012, 1)
      
      sst_djf=month_to_season(ssta_de, "DJF")
       sst_djf!0 = "time"
      sst_djf&time=ispan(1940,2012,1)  
      sst_djf(0:dimsizes(sst_djf&time)-2,:,:)=(/sst_djf(1:,:,:)/)
      sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/ssta_de({201212},:,:)/)
   

    ;;;

    
    ;;;回归风场和SST场
       rc_sst = regCoef_n(nino3_4_anu_std,sst_djf , 0, 0)
       rc_u = regCoef_n(nino3_4_anu_std,u_djf , 0, 0)
       rc_v = regCoef_n(nino3_4_anu_std,v_djf , 0, 0)

       copy_VarCoords(sst_djf(1,:,:), rc_sst)
       copy_VarCoords(u_djf(1,:,:), rc_u)
       copy_VarCoords(v_djf(1,:,:), rc_v)

       printVarSummary(rc_sst)
     
    ;;;

;;;;draw plot
;;---------------画趋势图--------------------------------------
    wks=gsn_open_wks("pdf", "U,V,SSTregressiononnino3.4")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "BlRe"

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    ; res@cnMinLevelValF       =   -0.3
    ; res@cnMaxLevelValF       =   0.3
    ; res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.6
    res@vpHeightF = 0.3

    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.03           ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefMagnitudeF          = 0.6
    res@vcRefLengthF             = 0.04     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"

    res@gsnLeftString = "Observation"
    plot = gsn_csm_vector_scalar_map(wks,rc_u,rc_v,rc_sst,res)
   
   draw(plot)
   frame(wks)
end    





;;;

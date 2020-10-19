;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
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

begin
season_r = (/"spring","summer","autumn","winter"/)
season = (/"MAM","JJA","SON","DJF"/)
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_atl_" + season_r +"_notest_wholeatl/"+"F_2000_atl_"+season_r+"_notest_wholeatl.cam.h1.0101-4012.nc" 
  
  ; wks=gsn_open_wks("pdf", "CESM-F2000-u,v,OMEGA-alt-warming-foul-season-global")
  
  ; plot = new(4,graphic)
  ; iplot = new(4,graphic)
  ; str3 = ispan(1,4,1) + ""
  m = (/2,5,8,0/)
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind:=ind(date.ge.601.and.date.le.4012)

  ;;;;U,V
  u_CT  =  f1->UBOT(timeind,:,:)
  v_CT  =  f1->VBOT(timeind,:,:)

  u_diff = u_CT(0:11,:,:)
  v_diff = u_diff

  printVarSummary(u_diff)

do i = 0,3

  f2=addfile(fils2(i), "r")
  print("读入file"+fils2(i))

  u_id  =  f2->UBOT(timeind,:,:)
  v_id  =  f2->VBOT(timeind,:,:)

  
  dim  = dimsizes(u_id)
  
  if season(i).ne."DJF" then
    print(m(i))
    do j = 0,2
      print(j+m(i))
      u_diff(j+m(i),:,:) = dim_avg_n_Wrap(u_id(j+m(i):dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(j+m(i):dim(0)-1:12,:,:),0)
      v_diff(j+m(i),:,:) = dim_avg_n_Wrap(v_id(j+m(i):dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(j+m(i):dim(0)-1:12,:,:),0)
    end do 
  else 
    ;   u_diff(2,:,:) = dim_avg_n_Wrap(u_id(2:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(2:dim(0)-1:12,:,:),0)
    ;   u_diff(3,:,:) = dim_avg_n_Wrap(u_id(3:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(3:dim(0)-1:12,:,:),0)
    ;   u_diff(4,:,:) = dim_avg_n_Wrap(u_id(4:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(4:dim(0)-1:12,:,:),0)
      
    ;   v_diff(2,:,:) = dim_avg_n_Wrap(v_id(2:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(2:dim(0)-1:12,:,:),0)
    ;   v_diff(3,:,:) = dim_avg_n_Wrap(v_id(3:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(3:dim(0)-1:12,:,:),0)
    ;   v_diff(4,:,:) = dim_avg_n_Wrap(v_id(4:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(4:dim(0)-1:12,:,:),0)
      
      u_diff(11,:,:) = dim_avg_n_Wrap(u_id(11:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(11:dim(0)-1:12,:,:),0)
      u_diff(0,:,:) = dim_avg_n_Wrap(u_id(0:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(0:dim(0)-1:12,:,:),0)
      u_diff(1,:,:) = dim_avg_n_Wrap(u_id(1:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(1:dim(0)-1:12,:,:),0)

      v_diff(11,:,:) = dim_avg_n_Wrap(v_id(11:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(11:dim(0)-1:12,:,:),0)
      v_diff(0,:,:) = dim_avg_n_Wrap(v_id(0:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(0:dim(0)-1:12,:,:),0)
      v_diff(1,:,:) = dim_avg_n_Wrap(v_id(1:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(1:dim(0)-1:12,:,:),0)
  end if 


    ; u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, season(i)), 0)
    ; u_id_clm = dim_avg_n_Wrap(month_to_season(u_id, season(i)), 0)

    ; u_diff = u_CT_clm
    ; u_diff = (/u_id_clm - u_CT_clm/)

    ; v_CT_clm = dim_avg_n_Wrap(month_to_season(v_CT, season(i)), 0)
    ; v_id_clm = dim_avg_n_Wrap(month_to_season(v_id, season(i)), 0)

    ; v_diff = v_CT_clm
    ; v_diff = (/v_id_clm - v_CT_clm/)
    u_diff@_FillValue = 10^9
    v_diff@_FillValue = 10^9
;;;;;

end do 

;;;;画图
 wks = gsn_open_wks("pdf", "diff-modify-wind")

;;;draw_plot 
  plot  =new(4,graphic)

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
 
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueRed"

  res@cnFillOn=True
  res@cnLinesOn= False
  
  res@cnLineThicknessF = 2.0
  res@gsnContourNegLineDashPattern = 14
  res@gsnContourPosLineDashPattern = 0

  ; res@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
  ; gsn_define_colormap(wks, "MPL_PuOr")

  ; res@cnLevels = (/-0.02,-0.015,-0.01,-0.005,0.005,0.01,0.015,0.02/)*100
  ; res@cnFillColors = (/22,32,45,60,-1,80,97,100,110/) 


  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@lbLabelStride = 2 
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -60
  res@mpMaxLatF             = 60
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.025
  res@tmYLLabelFontHeightF=0.025
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 

  
  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.022           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  ; res@vcRefAnnoJust = "TopLeft"
  res@vcRefAnnoBackgroundColor = "Transparent"
  res@vcRefAnnoPerimOn = False
  res@vcRefAnnoOrthogonalPosF = -1.2
  res@vcRefMagnitudeF          = 2.0
  res@vcRefLengthF             = 0.1     ; define length of v
  res@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =True
  res@mpFillColor = "White"
  res@mpFillDrawOrder = "PostDraw"


  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8

  ; res@gsnLeftString  = se(i)+" " + season(i)+" U,V at 992hPa trend"
  res@gsnLeftStringFontHeightF = 0.02
  res@gsnLeftStringOrthogonalPosF = -.15
  ; plot(i) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,slp_rc ,res)



  plot(0) = gsn_csm_vector_map(wks,dim_avg_n_Wrap(u_diff(2:4,:,:),0),dim_avg_n_Wrap(v_diff(2:4,:,:),0) ,res)
  plot(1) = gsn_csm_vector_map(wks,dim_avg_n_Wrap(u_diff(5:7,:,:),0),dim_avg_n_Wrap(v_diff(5:7,:,:),0) ,res)
  plot(2) = gsn_csm_vector_map(wks,dim_avg_n_Wrap(u_diff(8:10,:,:),0),dim_avg_n_Wrap(v_diff(8:10,:,:),0) ,res)
  plot(3) = gsn_csm_vector_map(wks,dim_avg_n_Wrap(u_diff(0:1,:,:),0),dim_avg_n_Wrap(v_diff(0:1,:,:),0) ,res)


pres =True
gsn_panel(wks, plot, (/2,2/), pres)





;;-----------------------写入nc文件---------------------------
    
    ; system("/bin/rm -f /home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind_"+season(i)+".nc")   ; remove any pre-existing file
    ; ncdf = addfile("/home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind_"+season(i)+".nc","c")  ; open output netCDF file


    system("/bin/rm -f /home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind_10m.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind_10m.nc","c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->u_diff  =  u_diff                           ; 4D               
       ncdf->v_diff  = v_diff                       ; 3D
       
end 

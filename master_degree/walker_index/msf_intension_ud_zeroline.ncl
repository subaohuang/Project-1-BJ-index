 
function gain_data(data1)
 begin
 
 rc1=regCoef(data1&time,data1)


 data_fuc=new((/3,dimsizes(data1&time)/),double)


 data_fuc(0,:) = data1(:)
 data_fuc(1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(2,:) =0.0

return(data_fuc)
end 

 begin
 diri="/media/sf_share/era/"
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 fils2=systemfunc("ls "+diri+"msf_era_level_location1979-2012.nc")
 f2=addfile(fils2, "r")
 locationxi=f2->locationxi
 locationdong=f2->locationdong
 printVarSummary(locationxi)
 level=f1->level(:)
 lon=f1->lon(:)

 ; lonind1=ind(lon.eq.150)
 ; lonind2=ind(lon.eq.210)
  ud_time=f1->ud_time

  ud_int=int2p_n_Wrap(level, ud_time, locationxi&level  , 1, 0)
  printVarSummary(ud_int)
  ud_new_up=ud_int(ind(ud_int&level.eq.300):ind(ud_int&level.eq.150),:,:)
  ud_new_down=ud_int(ind(ud_int&level.eq.1000):ind(ud_int&level.eq.750),:,:)
 ;**************************************************************
 ; m=ind_nearest_coord(level, locationxi&level, 0)
 ;  print("i="+m+"  rlat="+level+"  lat(i)="+locationxi&level(m))
 ; locationxi_s=locationxi(m,:)
 ; print(locationxi_s&level)  
 ;***************************************************************
 ; intension_ud=dim_avg_n_Wrap(ud_time(ind(level.eq.300):ind(level.eq.150),lonind1:lonind2,:),(/0,1/))-dim_avg_n_Wrap(ud_time(ind(level.eq.1000):ind(level.eq.750),lonind1:lonind2,:),(/0,1/))
 ; print(intension_ud)
 ;基于平均的结果，150E到210度E，150到300和1000到750的平均结果
 ;***************************************************************************
 locationxi_avg=dim_avg_n_Wrap(locationxi, 1)
 locationdong_avg=dim_avg_n_Wrap(locationdong, 1)

  do i=0,dimsizes(locationxi&level)-1
 locationxi(i,:)=where(ismissing(locationxi(i,:)), locationxi_avg(i),locationxi(i,:))                  ;没有位置的年份用气候态来代替
 locationdong(i,:)=where(ismissing(locationdong(i,:)), locationdong_avg(i),locationdong(i,:))          ;没有位置的年份用气候态来代替
 end do 

 locationxi_up=locationxi(ind(locationxi&level.eq.300):ind(locationxi&level.eq.150),:)
 locationdong_up= locationdong(ind(locationdong&level.eq.300):ind(locationdong&level.eq.150),:)

 locationxi_down=locationxi(ind(locationxi&level.eq.1000):ind(locationxi&level.eq.750),:)
 locationdong_down= locationdong(ind(locationdong&level.eq.1000):ind(locationdong&level.eq.750),:)

 i_xi_up=new((/dimsizes(locationxi_up&level),dimsizes(locationxi_up&time)/),integer)
 i_dong_up=new((/dimsizes(locationxi_up&level),dimsizes(locationxi_up&time)/),integer)

 i_xi_down=new((/dimsizes(locationxi_down&level),dimsizes(locationxi_down&time)/),integer)
 i_dong_down=new((/dimsizes(locationxi_down&level),dimsizes(locationxi_down&time)/),integer)

  do i=0,dimsizes(locationxi_up&level)-1
   i_xi_up(i,:)=ind_nearest_coord(locationxi_up(i,:), lon, 0)
   i_dong_up(i,:)=ind_nearest_coord(locationdong_up(i,:), lon, 0)
  end do
    print("i="+i_xi_up(1,:)+"  rlat="+locationxi_up(1,:)+"  lat(i)="+lon(i_xi_up(1,:)))

 ud_int_up=new((/dimsizes(ud_new_up&level),dimsizes(ud_new_up&time)/),float)
  do i=0,dimsizes(ud_new_up&time)-1
    do j=0,dimsizes(ud_new_up&level)-1  
    ud_int_up(j,i)=dim_avg_n(ud_new_up(j,i_xi_up(j,i):i_dong_up(j,i),i),0)
   end do
 end do


  do i=0,dimsizes(locationxi_down&level)-1
   i_xi_down(i,:)=ind_nearest_coord(locationxi_down(i,:), lon, 0)
   i_dong_down(i,:)=ind_nearest_coord(locationdong_down(i,:), lon, 0)
  end do
   
 ud_int_down=new((/dimsizes(ud_new_down&level),dimsizes(ud_new_down&time)/),float)
    do i=0,dimsizes(ud_new_down&time)-1
   do j=0,dimsizes(ud_new_down&level)-1  
    ud_int_down(j,i)=dim_avg_n(ud_new_down(j,i_xi_down(j,i):i_dong_down(j,i),i),0)
end do
end do

printVarSummary(ud_int_up)
printVarSummary(ud_int_down)
 
 intension_ud=dim_avg_n(ud_int_up, 0)-dim_avg_n(ud_int_down, 0)

 ; copy_VarCoords(ud_new(), var_to)
 intension_ud!0="time"
 intension_ud&time=ispan(1,457,1)
 intesion_an=runave_n_Wrap(dim_rmvmean(month_to_annual(intension_ud(0:455), 1)),3,0,0)
 intesion_an!0="time"
 intesion_an&time=ispan(1,38,1)
 data_an=gain_data(intesion_an)
 print(data_an)
 wks1=gsn_open_wks("png", ""+diri+"msf_intension_era_ud_chazhi")
 ;************************************************************
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnLeftString="era"
    resxy@gsnRightString="intesion"
    resxy@tmXBMode="Explicit"
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(intension_ud&time)
    resxy@tmXBValues=ispan(13,dimsizes(intension_ud&time),60)
    resxy@tmXBLabels=ispan(1980,2015,5)

    resxy@vpHeightF=0.7
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9

    resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
    resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    resxy@xyDashPatterns    =(/0,0,2/)
    

    resxy@tmXBLabelFontAspectF=2.5
    resxy@tmYLLabelFontAspectF=2.5
    resxy@tmXTOn=False
    resxy@tmYROn=False
    ; resxy@tiYAxisString = title

    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy@pmLegendSide           = "Top"               ; Change location of 
    resxy@pmLegendParallelPosF   = 0.85                  ; move units right
    resxy@pmLegendOrthogonalPosF = -0.18                ; move units down
    resxy@pmLegendWidthF         = 0.10                ; Change width and
    resxy@pmLegendHeightF        = 0.1               ; height of legend.
    resxy@lgPerimOn              = False               ; turn off box around
    resxy@lgLabelFontHeightF     = .008                 ; label font height
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "linear trend"+"("+sprintf("%5.3f",rc1(1)*120)+"~S~o~N~/decade"+")","average"/)

    delete(resxy@tmXBValues)
    delete(resxy@tmXBLabels)
    delete(resxy@trXMaxF)
    resxy@tmXBValues=ispan(1,38,5)
    resxy@tmXBLabels=ispan(1979,2016,5)
    resxy@tmXBMinorOn=True
    resxy@tmXBMinorValues=ispan(1,38,1)

;*****************************************************************************************
   plot=gsn_csm_xy(wks1,intesion_an&time,data_an,resxy)
   draw(plot)
   frame(wks1)



 end 
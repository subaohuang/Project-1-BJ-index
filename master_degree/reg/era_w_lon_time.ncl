function gain_data(data1,data2)
 begin
 
 rc1=regCoef(data1&time,data1)

 rc2=regCoef(data2&time,data2)
  print(rc1)
  print(rc2)
 data_fuc=new((/2,4,dimsizes(data1&time)/),float)


 data_fuc(0,0,:) = data1(:)
 data_fuc(0,1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(0,2,:) =avg(data1(:))
 data_fuc(0,3,:)= rc1
 data_fuc(1,0,:) = data2(:)
 data_fuc(1,1,:) = rc2*data2&time+rc2@yintercept
 data_fuc(1,2,:) =avg(data2(:))
 data_fuc(1,3,:) =rc2
return(data_fuc)
end 

function draw_panel(wks1,wks2,time,data1,data2,title,rc1,rc2)
 begin
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnLeftString="era"
    resxy@gsnRightString="upward location"
    resxy@tmXBMode="Explicit"
    ; resxy@trXMaxF=456
    resxy@tmXBValues=ispan(100,240,20)
    resxy@tmXBLabels=ispan(100,240,20)+"~S~o~N~E"
    resxy@tmXBMinorOn=False
    
    
    resxy@tmYLMode="Explicit"
    
    resxy@trYMaxF=dimsizes(data1&time)
    resxy@tmYLValues=ispan(13,dimsizes(data1&time),60)
    resxy@tmYLLabels=ispan(1980,2015,5)

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
    resxy@tiYAxisString = title

    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy@pmLegendSide           = "Top"               ; Change location of 
    resxy@pmLegendParallelPosF   = 0.85                  ; move units right
    resxy@pmLegendOrthogonalPosF = -0.18                ; move units down
    resxy@pmLegendWidthF         = 0.10                ; Change width and
    resxy@pmLegendHeightF        = 0.1               ; height of legend.
    resxy@lgPerimOn              = False               ; turn off box around
    resxy@lgLabelFontHeightF     = .008                 ; label font height
    resxy@xyExplicitLegendLabels = (/"Anual change", "linear trend"+"("+sprintf("%5.3f",rc1(1)*120)+"~S~o~N~/decade"+")","average"/)         ; create explicit labels
    




    ; resxy@tmLabelAutoStride = True                    ; nice stride on labels
    contour=new(2,graphic)
    resxy2=resxy
    delete(resxy2@tmXBValues)
    delete(resxy2@tmXBLabels)
    delete(resxy2@xyExplicitLegendLabels)

    delete(resxy2@pmLegendDisplayMode)
    delete(resxy2@pmLegendSide)
    delete(resxy2@pmLegendParallelPosF   )
    delete(resxy2@pmLegendOrthogonalPosF )
    delete(resxy2@pmLegendWidthF        )
    delete(resxy2@pmLegendHeightF       )
    delete(resxy2@lgPerimOn              )
    delete(resxy2@lgLabelFontHeightF    )
    delete(resxy2@gsnLeftString)
    delete(resxy2@gsnRightString)


    resxy2@gsnLeftString="era"
    resxy2@gsnRightString="downward location"
    resxy2@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy2@pmLegendSide           = "Top"               ; Change location of 
    resxy2@pmLegendParallelPosF   = 0.15                  ; move units right
    resxy2@pmLegendOrthogonalPosF = -0.18                ; move units down
    resxy2@pmLegendWidthF         = 0.10                ; Change width and
    resxy2@pmLegendHeightF        = 0.1               ; height of legend.
    resxy2@lgPerimOn              = False               ; turn off box around
    resxy2@lgLabelFontHeightF     = .008                 ; label font height
    resxy2@xyExplicitLegendLabels = (/"Anual change", "linear trend"+"("+sprintf("%5.3f",rc2(1)*120)+"~S~o~N~/decade"+")","average"/)         ; create explicit labels
    

    resxy2@tmXBValues=ispan(180,280,20)
    resxy2@tmXBLabels=""+ispan(180,280,20)+"~S~o~N~E"
  

    contour(0)=gsn_csm_xy(wks1,data1(:,:),time,resxy)
    contour(1)=gsn_csm_xy(wks2,data2(:,:),time,resxy2)
    
 
    ; contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    ; contour(1)=gsn_csm_xy(wks, time,data2,resxy)

return contour
end 
;*****************************************************************************************
begin 

    diri="/media/sf_share/era/"     ;读取数据  
    fils1=systemfunc("ls "+diri+"era_w_rq.nc")
    f1=addfile(fils1, "r")
   
    time=f1->time
    date=cd_calendar(time, 1)
    lat=f1->latitude
    latind1=ind(lat.eq.5)
    latind2=ind(lat.eq.-5.0)         ;没去气流的都是-5到5°,去气流的是从20到-20度。
    tind1=ind(date.eq.197901)
    tind2=ind(date.eq.201701)

    w_r=f1->w(tind1:tind2,1,latind1:latind2,:)     ;到2014年12月
    w=w_r@scale_factor*w_r+w_r@add_offset
    w_new=doubletofloat(w)
    copy_VarMeta(w_r, w_new)

    ;去年循环-------------------------------------------------

      ; w_ano= rmMonAnnCycLLT(w_new(latitude|:,longitude|:,time|:))
      ; w_ano_new=w_ano(time|:,latitude|:,longitude|:)
      ;  w_avg=dim_avg_n_Wrap(w_ano_new, 1)
     ;-----------------------------------------------------------

     w_avg=dim_avg_n_Wrap(w_new, 1)
    ; w_avg=where(w_avg.gt.-0.025, w_avg@missing_value, w_avg)
    ; write_matrix(w_avg,"381f8.4",False)
    w_xi=w_avg(:,0:ind(w_avg&longitude.eq.240.0))
    w_dong=w_avg(:,ind(w_avg&longitude.eq.180.0):ind(w_avg&longitude.eq.280.0))


    ;去气流的是从20到-20度。------------------------------
    ; w_xi=where(w_xi.gt.-0.025,w_xi@missing_value , w_xi)
    ; w_dong=where(w_dong.lt.-0.013,w_dong@missing_value , w_dong)
    ;-------------------------------------------------------

    

     print(w_dong(365,:))
    
     location1=new((/dimsizes(w_r&time)/),float)
     location2=new((/dimsizes(w_r&time)/),float)

   
     do i=0,dimsizes(w_avg&time)-1,1
          location1(i)=w_xi&longitude(minind(w_xi(i,:)))
          location2(i)=w_dong&longitude(maxind(w_dong(i,:)))   
     end do
    
     

      location1!0="time"
      location2!0="time"
     
      location1&time=ispan(1,dimsizes(w_r&time),1)
      location2&time=ispan(1,dimsizes(w_r&time),1)
    
      data= gain_data(location1,location2)



    wks1=gsn_open_wks("png", "/media/sf_share/era/era_w_series_xi.2017_没去气流")
    wks2=gsn_open_wks("png", "/media/sf_share/era/era_w_series_dong.2017_没去气流")

    pres1=True
    pres2=True
    pres1@gsnMaximize=True
    contour=draw_panel(wks1,wks2,location1&time,data(0,0:2,:),data(1,0:2,:),"",data(0,3,:),data(1,3,:) )

    draw(contour(0))
    draw(contour(1))
    frame(wks1)
    frame(wks2)

  end
 
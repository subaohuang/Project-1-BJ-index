function gain_data(data1,data2)
 begin
 
 rc1=regCoef(data1&time,data1)  ;; 第一个位置是x,第二个位置是y。

 rc2=regCoef(data2&time,data2)

 data_fuc=new((/2,3,dimsizes(data1&time)/),float)


 data_fuc(0,0,:) = data1(:)
 data_fuc(0,1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(0,2,:) =avg(data1(:))

 data_fuc(1,0,:) = data2(:)
 data_fuc(1,1,:) = rc2*data2&time+rc2@yintercept
 data_fuc(1,2,:) =avg(data2(:))
return(data_fuc)
end 

function draw_panel(wks1,wks2,time,data1,data2,title)
 begin
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnLeftString="(a)"

    resxy@tmXBMode="Explicit"
    ; resxy@trYMaxF=456
    resxy@tmXBValues=ispan(100,240,20)
    resxy@tmXBLabels=ispan(100,240,20)+"~S~o~N~E"
    resxy@tmXBMinorOn=False 
    resxy@tmYLMode="Explicit"
    
    resxy@trYMaxF=dimsizes(data1&time)
    resxy@tmYLValues=ispan(13,dimsizes(data1&time),60)
    resxy@tmYLLabels=ispan(1980,2015,5)
    

    resxy@vpHeightF=0.7
    resxy@vpWidthF=0.5
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
    
    ; resxy@tmLabelAutoStride = True                    ; nice stride on labels
    contour=new(2,graphic)
    resxy2=resxy
    delete(resxy2@tmXBValues)
    delete(resxy2@tmXBLabels)

    resxy2@tmXBValues=ispan(180,280,20)
    resxy2@tmXBLabels=""+ispan(180,280,20)+"~S~o~N~E"
    resxy2@gsnLeftString="(b)"
    contour(0)=gsn_csm_xy(wks1,data1,time,resxy)

    contour(1)=gsn_csm_xy(wks2,data2,time,resxy2)


    ; contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    ; contour(1)=gsn_csm_xy(wks, time,data2,resxy)

return contour
end 
;*****************************************************************************************
begin 

    diri="/mnt/f/share/era/data/shixiong/"     ;读取数据  
  
    fils1=systemfunc("ls "+diri+"era_iod_pacific_w.nc")
    f1=addfile(fils1, "r")
   
    time=f1->time
    date=cd_calendar(time, 1)
    lat=f1->latitude
    latind1=ind(lat.eq.10)
    latind2=ind(lat.eq.-10)
     ;w_r=f1->w(0:455,1,20:100,:)    ;到20171月
    w_r=f1->w(0:431,latind1:latind2,{100:290})     ;到2014年12月
    w=w_r@scale_factor*w_r+w_r@add_offset
    w_new=doubletofloat(w)
    copy_VarMeta(w_r, w_new)
    delete(w_r)
    ; w_avg=where(w_avg.gt.-0.025, w_avg@missing_value, w_avg)
     ; write_matrix(w_avg,"381f8.4",False)
     printVarSummary(w_new)


    ;;;去除掉陆地上的纬向风的影响
     fils2=systemfunc("ls "+diri+"landsea.nc")
     a=addfile(fils2, "r")
     lsdata = a->LSMASK
     lsm  = landsea_mask(lsdata,w_new&latitude,w_new&longitude)   ;;得到的结果，1表示land   
     w_new = mask(w_new,lsm.eq.1,False)      
   ;;;
      ; w_xi=w_new(:,:,{100:240})
      ; w_dong=w_new(:,:,{180:280})
       w_avg=dim_avg_n_Wrap(w_new, 1)
      ; printVarSummary(w_avg)
     ;------------------------------------------
     
      ; w_xi=where(w_xi.gt.-0.025, w_xi@missing_value, w_xi)
      ; w_dong=where(w_dong.lt.-0.013, w_dong@missing_value, w_dong)
      
      ; w_xi_avg=dim_avg_n_Wrap(w_xi, 1)
      ; w_dong_avg=dim_avg_n_Wrap(w_dong, 1)

      w_avg_xi=w_avg(:,{100:240})
      w_avg_dong=w_avg(:,{180:280})
       ; w_avg_xi=where(w_avg_xi.gt.-0.025,w_avg@missing_value,w_avg_xi)
       ; w_avg_dong=where(w_avg_dong.lt.0.013,w_avg@missing_value,w_avg_dong)

      location1=new((/dimsizes(w_new&time)/),float)
      location2=new((/dimsizes(w_new&time)/),float)

    do i=0,dimsizes(w_new&time)-1,1
          location1(i)=w_avg_xi&longitude(minind(w_avg_xi(i,:)))
          location2(i)=w_avg_dong&longitude(maxind(w_avg_dong(i,:)))   
    end do

    
      location1!0="time"
      location2!0="time"
     
      location1&time=ispan(1,dimsizes(w_new&time),1)
      location2&time=ispan(1,dimsizes(w_new&time),1)
      
      data= gain_data(location1,location2)



    wks1=gsn_open_wks("png", "/mnt/f/share/era/era_w_series_xi.ceshi")
    wks2=gsn_open_wks("png", "/mnt/f/share/era/era_w_series_dong.ceshi")

    pres1=True
    pres2=True
    pres1@gsnMaximize=True
    contour=draw_panel(wks1,wks2,location1&time,data(0,:,:),data(1,:,:),"year" )

    draw(contour(0))
    draw(contour(1))
    frame(wks1)
    frame(wks2)



   opt=diri+"locationxi.ascii"

   asciiwrite(opt, sprintf ("%5.2f,", location1))

  opt=diri+"locationdong.ascii"

   asciiwrite(opt, sprintf ("%5.2f,", location2))


  end
 
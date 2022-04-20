function  ploy(locationxi_rea)
begin
dim=dimsizes(locationxi_rea)
dim1=dim(1)
mnmx      = new ((/3,dim(1)/), typeof(locationxi_rea))
mnmx(0,:) = dim_avg_n_Wrap(locationxi_rea, 0)
mnmx(1,0:dim1-1) = dim_min_n_Wrap(locationxi_rea, 0)
mnmx(2,0:dim1-1) = dim_max_n_Wrap(locationxi_rea, 0)
return mnmx
end 

function level_mean(data)
local data,msf_mean
begin
plev=data&level
psfc=conform_dims((/dimsizes(data&time),dimsizes(data&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)

msfdp=data(level|:,time|:,lon|:)*dp
copy_VarCoords(data(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(data(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 

undef("search_location")
function  search_location(msf_mean_annual,mode)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:290})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+8).gt.0))  then 
         a=(msf_mean_annual_xi(t,i+1)-msf_mean_annual_xi(t,i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi(t)=tofloat(msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(t,i)/a)
        ; locationxi(t)=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if ismissing(locationxi(t)).ne.True  then 
        break
        end if
     end do
end do

;;;judge locationdong
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0))  then 
         a=(msf_mean_annual_dong(t,i+1)-msf_mean_annual_dong(t,i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong(t)=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a)
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
  if mode.eq.0 then 

   if ismissing(locationdong(t)).eq.True  then 
   co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,{230:270}))
   locationdong(t)=msf_mean_annual_dong&lon(co_in)
   delete(co_in)
  end if

end if
end do


locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)

 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

;********************************************************
;draw_msf_hov
;*******************************************************
function draw_msf_hov(wks,LeftString,title)
begin 
print("************************"+LeftString+"**************************")
    diri="/mnt/f/share/rea_msf/"
    fils1=diri+"msf_time_"+LeftString+"_global.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
   
    msf_time_month=clmMonLLT(f1->msf_time(:,:,timeind))
    msf_time_month=(/msf_time_month/10^9/)
    level=msf_time_month&level

    levind=ind(msf_time_month&level.eq.500)
    levind2=ind(msf_time_month&level.le.600.and.msf_time_month&level.ge.400)

    printVarSummary(msf_time_month)
    msf_time_month!2="time"
    
    msf_mean=level_mean(msf_time_month(time|:,lon|:,level|levind2))
    ; msf_mean=dim_avg_n_Wrap(msf_time_month(time|:,lon|:,level|levind2),2)
    ; location=search_location(msf_time_month(time|:,lon|:,level|levind),0)
    location=search_location(msf_mean,0)
    

    locationxi=location[0]
    locationdong=location[1]
    printVarSummary(locationxi)

    ; location_nan=search_location(msf_time_month(time|:,lon|:,level|levind),1)  
    location_nan=search_location(msf_mean,1)  
    locationdong_nan=location_nan[1]

    msf_mean_month=level_mean(msf_time_month(time|:,lon|:,level|:))    ;;msf_mean_month(time,lon)

    ; msf_mean_month=msf_time_month(time|:,level|levind,lon|:)


    printVarSummary(msf_mean_month)
    msf_mean_month_new=msf_mean_month
    ; msf_mean_month_new(0:5,:)=(/msf_mean_month(6:11,:)/)
    ; msf_mean_month_new(6:11,:)=(/msf_mean_month(0:5,:)/)
     
    locationxi_new=locationxi
    ; locationxi_new(0:5)=(/locationxi(6:11)/)
    ; locationxi_new(6:11)=(/locationxi(0:5)/)
    
    locationdong_new=locationdong
    ; locationdong_new(0:5)=(/locationdong(6:11)/)
    ; locationdong_new(6:11)=(/locationdong(0:5)/)
    
    locationdong_nan_new=locationdong_nan
    ; locationdong_nan_new(0:5)=(/locationdong_nan(6:11)/)
    ; locationdong_nan_new(6:11)=(/locationdong_nan(0:5)/)
    

    width_new=tofloat(locationdong_new - locationxi_new)
    ; width_new=tofloat(locationxi_new )
    inten=dim_avg_n_Wrap(msf_mean_month_new(:,{180:240}),1)
    
    




;***********************************************
;draw hov
;***********************************************

    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
  
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame
    res@gsnLeftString=title
   	res@cnLineLabelsOn =False
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnFillPalette="GMT_polar"
  

	res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	res@cnMinLevelValF       = -7.0              ; min level
	res@cnMaxLevelValF       =  7.0               ; max leve
	res@cnLevelSpacingF      =1.0


	res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
	res@tmYLLabelsOn=True
	res@tmYUseLeft=True
	res@tmYLValues    = ispan(0,11,1)
	res@tmYLLabels    = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
    ; res@tmYLLabels    = (/"J","A","S","O","N","D","J","F","M","A","M","J"/)
	res@tmXTOn=False
	res@tmYROn=False

	; res@tmXBValues=(/0,90,180,270,357.5/)
 ;    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
 ;    res@tmXBMinorOn=True
 ;    res@tmXBMinorValues=ispan(0,360,10)

    res@tmXBValues=ispan(90, 280, 30)
    res@tmXBLabels=(/"90~S~o~N~E","120~S~o~N~E","150~S~o~N~E","180","150~S~o~N~W","120~S~o~N~W"\
                    ,"90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(90,280,10)


	res@vpHeightF=0.4
	res@vpWidthF=0.6
    plot = gsn_csm_contour(wks, msf_mean_month_new(:,{90:280}), res)   ; create plot


    ; res@cnLinesOn=True
    ; res@cnFillOn        = False
    ; res@cnLevelSelectionMode="ExplicitLevels"
    ; res@cnLevels=0.0
    ; res@cnLineThicknessF=2.0
 
    ; plot1 = gsn_csm_contour(wks, msf_mean_month_new, res)
    
    polyres=True
    polyres@gsLineColor="black"
    polyres@gsLineDashPattern=2
    polyres@gsLineThicknessF=3.0
    polyres@tfPolyDrawOrder="PreDraw"
    str1=unique_string("dum1")
    str2=unique_string("dum2")
    str3=unique_string("dum3")
    str4=unique_string("dum4")
    plot@$str1$=gsn_add_polyline(wks, plot, locationdong_new,locationdong&time, polyres)
    plot@$str2$=gsn_add_polyline(wks, plot, locationxi_new,locationxi&time, polyres)
   
    polyres@tfPolyDrawOrder="PostDraw"
    polyres@gsLineColor="black"
    polyres@gsLineDashPattern=0
    plot@$str3$=gsn_add_polyline(wks, plot, locationdong_nan_new ,locationdong&time, polyres)
    plot@$str4$=gsn_add_polyline(wks, plot, locationxi_new,locationxi&time, polyres)
    ; overlay(plot, plot1)
    plot@width_new=width_new
    plot@inten=tofloat(inten)
    plot@locationxi=tofloat(locationxi)
    plot@locationdong=tofloat (locationdong)
     r=escorc(width_new, inten)
    plot@r=tofloat( r)
 return plot
 end 

 ;***************************************************************
 ;main program
 ;****************************************************************
  begin 
   ; rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2","mre"/)
   rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2"/)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料沃克环流月气候态_小范围")
    title=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)
    plot=new(6,graphic)
    width_rea=new((/7,12/),float)
    locationxi_rea=new((/7,12/),float)
    locationdong_rea=new((/7,12/),float)
    inten_rea=new((/7,12/),float)
      r_rea=new(7,float)
    do i=0,dimsizes(plot)-1
    plot(i)=draw_msf_hov(wks,rea(i),title(i))
    plot_new=plot(i)
    width_rea(i,:)=plot_new@width_new
    locationxi_rea(i,:)=plot_new@locationxi
    locationdong_rea(i,:)=plot_new@locationdong
    inten_rea(i,:)=plot_new@inten
    r_rea(i)=plot_new@r
    end do
   printVarSummary(width_rea)
    width_rea(6,:)=dim_avg_n(width_rea(0:5,:), 0)
    locationxi_rea(6,:)=dim_avg_n(locationxi_rea(0:5,:), 0)
    locationdong_rea(6,:)=dim_avg_n(locationdong_rea(0:5,:), 0)
    inten_rea(6,:)=dim_avg_n(inten_rea(0:5,:), 0)
    r_rea(6)=dim_avg_n(r_rea(0:5), 0)
    print("r="+r_rea)

    locationxi_jja=avg(locationxi_rea(6,5:7))
    locationxi_djf=(locationxi_rea(6,0)+locationxi_rea(6,1)+locationxi_rea(6,11))/3.0
    locationdong_jja=avg(locationdong_rea(6,5:7))
    locationdong_djf=(locationdong_rea(6,0)+locationdong_rea(6,1)+locationdong_rea(6,11))/3.0
    locationdong_mam=avg(locationdong_rea(6,2:4))

    width_djf=(width_rea(6,0)+width_rea(6,1)+width_rea(6,11))/3.0
    width_jja=avg(width_rea(6,5:7))
    width_son=avg(width_rea(6,8:10))
    print("width_djf="+width_djf)
    print("width_jja="+width_jja)
    print("width_jja="+width_son)
    print("locationxi_djf="+locationxi_djf)
    print("locationxi_jja="+locationxi_jja)

    print("locationdong_djf="+locationdong_djf)
    print("locationdong_jja="+locationdong_jja)

    print("locationdong_mam="+locationdong_mam)
    ; width_rea=(/dim_rmvmean_n(width_rea,1)/)
    ; inten_rea=(/dim_rmvmean_n(inten_rea,1)/)
  
    pres=True
    pres@gsnMaximize=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.6
    pres@pmLabelBarHeightF=0.05
    pres@pmLabelBarOrthogonalPosF =-0.01
    gsn_panel(wks, plot, (/4,2/), pres)

    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料沃克环流位置_年循环")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    
    res@xyDashPattern=(/0,0,0,0,0,0,1/)
    res@xyLineThicknesses=(/2,2,2,2,2,2,4/) 
    res@xyLineColors = (/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black"/)
    res@tmXBMode="Explicit"
    ; res@tmXBLabels    = (/"J","A","S","O","N","D","J","F","M","A","M","J"/)
    res@tmXBLabels    = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
    res@trXMaxF=11
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBValues=ispan(0, 11, 1)
     
    res@trYMaxF=170.0
    res@trYMinF=130.0

    res@vpHeightF=0.3
    res@vpWidthF=0.7
    inten_line=True
    inten_line@gsLineColor="black"
    inten_line@gsLineDashPattern=1
    plot_line=new(3,graphic)
    ; plot_line(0)=gsn_csm_xy(wks1, ispan(0, 11, 1), width_rea,inten_rea(6,:),res,inten_line)
    res@gsnLeftString="(a)Western edge"
    plot_line(0)=gsn_csm_xy(wks1, ispan(0, 11, 1), locationxi_rea,res)
    delete(res@trYMaxF)
    delete(res@trYMinF)
    res@gsnLeftString="(b)Eastern edge"
    plot_line(1)=gsn_csm_xy(wks1, ispan(0, 11, 1), locationdong_rea,res)
    res@gsnLeftString="(c)Width"
    plot_line(2)=gsn_csm_xy(wks1, ispan(0, 11, 1), width_rea,res)


    text=(/"20CR","ERA","JRA-55","MERRA","NCEP1","NCEP2","ENS"/)
    x=fspan(0.5, 10.5, 7)
    ; y=(/127,127,127,127,127,127,127/)
    y=(/168,168,168,168,168,168,168/)
    color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
    m=new(7,graphic)
    textres=True
    textres@txFontHeightF =12.0
    do i=0,6
    textres@txFontColor=color(i)
    m(i)=gsn_add_text(wks1,plot_line(0),text(i),x(i),y(i),textres)
    end do
    inten_line=True
     inten_line2=True
     inten_line@xyLineColor="red"
    ; inten_line@gsLineDashPattern=1
    ; plot_inten=gsn_add_polyline(wks1, plot_line, ispan(0,11,1), inten_rea(6,:), inten_line)
    inten_rea_month=new(4,typeof(inten_rea))
    width_rea_month=new(4,typeof(width_rea))

 
        width_rea_month(0)=avg(width_rea(6,2:4))
        width_rea_month(1)=avg(width_rea(6,5:7))
        width_rea_month(2)=avg(width_rea(6,8:10))
        width_rea_month(3)=(width_rea(6,0)+width_rea(6,1)+width_rea(6,11))/3.0

        inten_rea_month(0)=avg(inten_rea(6,2:4))
        inten_rea_month(1)=avg(inten_rea(6,5:7))
        inten_rea_month(2)=avg(inten_rea(6,8:10))
        inten_rea_month(3)=(inten_rea(6,0)+inten_rea(6,1)+inten_rea(6,11))/3.0
        delete(res)
;;;画强度与宽度季节上的关系

        mnmx_width=ploy(width_rea(0:5,:))            ;ploy(locationxi_rea,intension_rea)
   printVarSummary(mnmx_width)

        res                    = True              ; plot mods desired
        res@gsnDraw            = False             ; don't draw yet
        res@gsnFrame           = False             ; don't advance frame yet

        
        res@tmXTOn=False
        res@xyMonoLineColor    = False             ; want colored lines
        res@xyLineColors       = (/"Blue"/) ; colors chosen
        res@tmXBMode="Explicit"
        res@tmXBLabels    = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
        res@trXMaxF=11
        res@tmXTOn=False
        res@tmXBValues=ispan(0, 11, 1)
        res@trYMaxF=130.0
        res@trYMinF=80.0
        res@xyLineThicknessF=5.0
        res@tiYAxisString="PWCW(~S~o~N~)"
        res@tiYAxisFontHeightF=0.015
        res@vpHeightF=0.3
        res@vpWidthF=0.7
        res@vpXF=0.1
        inten_line=True
        inten_line@xyLineColors="red"
        inten_line@xyLineThicknessF=5.0
        inten_line@trYMaxF=4.5
        inten_line@trYMinF=1.5
        inten_line@tiYAxisString="PWC Intensity(10~S~9~N~Kg s~S~-1)"
        inten_line@tiYAxisFontHeightF=0.015
        mid_plot= gsn_csm_xy2(wks1,ispan(0,11,1),mnmx_width(0,:),inten_rea(6,:),res,inten_line)       ; create line plot

        delete(res@xyLineColors)
        res@gsnXYFillColors = "Blue"
        res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
        res@gsnXYFillOpacities=0.1

        amip_plot=gsn_csm_xy(wks1,ispan(0,11,1),mnmx_width(1:2,:),res)       ; create line plot
        
        ;;添加legend 
        gres = True
        gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
        gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

        lineres = True
        lineres@lgLineColors = (/"blue","red"/) ; line colors
        lineres@lgLineThicknesses = 5.0                       ; line thicknesses
        lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

        textres = True
        textres@lgLabels = (/"PWCW","PWC Intensity"/)
        textres@lgLabelFontHeights = (/0.012,0.012/)                           ; label font heights
  


        plot = simple_legend(wks1,mid_plot,gres,lineres,textres)   



        overlay(mid_plot, amip_plot)
        draw(mid_plot)
        frame(wks1)

        
        ; pres1=True
        ; gsn_panel(wks1, plot_line, (/3,1/), pres1)
           wks2=gsn_open_wks("pdf", "/mnt/f/share/era/data/年循环")
delete(res)
  res=True
  res@xyMarkLineModes       = "MarkLines"  ; markers and lines
  res@xyMonoMarkLineMode=True
  res@xyMarkers            = 16           ; filled dot
  res@xyMarkerSizeF        = 0.0135
plot = gsn_csm_xy(wks2,inten_rea(6,:),width_rea(6,:),res)  

    end
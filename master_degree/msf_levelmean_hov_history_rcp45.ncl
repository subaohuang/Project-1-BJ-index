
function  search_location(msf_mean_annual)
begin

locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:220})
msf_mean_annual_dong=msf_mean_annual(:,{230:300})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-13
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+12).gt.0))  then 
        locationxi(t)=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
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
        locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
end do
locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)

 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

;********************************************************
;draw_msf_hov_history puls RCP45
;1950-2099年
begin
 
 model="CanESM2"
 ; model="ACCESS1-0"
 diri1="/mnt/f/share/CMIP5/msf_history/"
 fils1=systemfunc("ls "+diri1+"msf_time_"+model+"_history_r1i1p1_195001-200512.nc")
  
 diri2="/mnt/f/share/CMIP5/msf_rcp45/"+model+"/"
 fils2=systemfunc("ls "+diri2+"msf_time_"+model+"_rcp45_r1i1p1_200601-209912.nc")
;;;********************************************* 
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    timeind1=ind(date.le.200512.and.date.ge.195001)
    msf_time_history=f1->msf_time({500},:,timeind1)
    delete(f1)
    delete(date)
 
    date=cd_calendar(f2->time, 1)
    timeind2=ind(date.le.209912.and.date.ge.200601)
    msf_time_rcp45=f2->msf_time({500},:,timeind2)
    delete(f2)
    delete(date)
     

    year=(dimsizes(timeind1)+dimsizes(timeind2))/12-1
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    msf_time=new((/dimsizes(msf_time_history&lon),dimsizes(timeind1)+dimsizes(timeind2)/),typeof(msf_time_history))
    msf_time(:,0:dimsizes(msf_time_history&time)-1)=msf_time_history
    msf_time(:,dimsizes(msf_time_history&time):)=msf_time_rcp45
    msf_time=(/msf_time/10^9/)
    msf_time!0="lon"
    msf_time!1="time"
    
    msf_time&time=ispan(0,dimsizes(timeind1)+dimsizes(timeind2)-1,1)
    msf_time&lon=msf_time_history&lon
    delete(msf_time_history)
    delete(msf_time_rcp45)
    timeind=msf_time&time
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time(time|:,lon|:)
    printVarSummary(msf_mean)
 
    msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
    printVarSummary(msf_mean_annual)

    msf_mean_annual=month_to_annual(msf_mean, 1)

    msf_mean_annual!0="time"
    msf_mean_annual!1="lon"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    msf_mean_annual&lon=msf_mean&lon

    printVarSummary(msf_mean_annual)
    msf_mean_annual!0="time"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)


    location=search_location(msf_mean_annual)
    locationxi=location[0]
    locationdong=location[1]

    rc_xi=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),locationxi , 0, 0)
    rc_xi_line=rc_xi*ispan(0,dimsizes(locationxi)-1,1) + rc_xi@yintercept

    rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
    rc_dong_line=rc_dong*ispan(0,dimsizes(locationdong)-1,1) + rc_dong@yintercept

    width=locationdong - locationxi

    
    width=(/dim_rmvmean_n(width,0)/)
    ; width=(/dim_standardize_n(width, 1, 0)/)
    width=(/runave(width, 3, -1)/)
    time=locationxi&time

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
    rc_width_line=rc_width*ispan(0,dimsizes(locationxi)-1,1) + rc_width@yintercept
    rc_decade=(/rc_width*10/)
  


;***********************************************
;draw hov
;***********************************************

   wks=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/msf_hov_history_rcp45/"+model+"_1958-2012_500hpa_msf_history_hov")

    res            = True                         ; plot mods desired
    res@gsnDraw    = False                        ; don't draw
    res@gsnFrame   = False                        ; don't advance frame

    res@cnLineLabelsOn =False
    res@cnFillDrawOrder="PreDraw"
    res@cnLinesOn=False
    ; res@cnFillPalette="BlueDarkRed18"
    res@cnFillPalette="testcmap"
    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = -7.0               ; min level
    res@cnMaxLevelValF       =  7.0               ; max leve
    res@cnLevelSpacingF      =0.5
    res@lbOrientation         ="Vertical"

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@tmYUseLeft=True

    res@tmYLValues    = ispan(0,dimsizes(timeind)/12-1,10)    
  res@tmYLLabels    = ispan(1950,2100,10)
  ; res@tmYLValues    = ispan(2,dimsizes(timeind)/12-1,5)    
  ;   res@tmYLLabels    = ispan(1958,2012,5)
    res@tmXTOn=False
    res@tmYROn=False

    ; res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBValues=(/0,90,180,270/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
    res@tmXBMinorOn=True

    res@vpHeightF=0.25
    res@vpWidthF=0.4
    plot = gsn_csm_contour(wks, msf_mean_annual, res)   ; create plot


    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = 0
    res@cnLineThicknessF = 4.
    res@cnFillOn        = False

    plot1 = gsn_csm_contour(wks, msf_mean_annual, res)
    
    overlay(plot, plot1)

    ;;;;;add trend
    lineres=True
    lineres@gsLineColor="red"
    lineres@gsLineDashPattern=1
    lineres@gsLineThicknessF=3.0
   
    str1 = unique_string("dum1")  ; "unique_string" will return a unique
                                     ; string every time it is called from
                                     ;  within a single NCL session.
    str2 = unique_string("dum2")

    plot@$str1$=gsn_add_polyline(wks, plot, rc_xi_line, ispan(0,year,1), lineres)
    plot@$str2$=gsn_add_polyline(wks, plot, rc_dong_line, ispan(0,year,1), lineres)

    ;;;;draw_attch_plot***************************************************************************
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False


    resxy@tmXTOn=True
    resxy@tmYLOn=False
    resxy@tmYROn=False
    resxy@tmXTOn=False
    resxy@trYMaxF=dimsizes(timeind)/12-1
   
    ; resxy@trXMinF=-2.0
    ; resxy@trXMaxF=3.0
    ; resxy@tmXBValues=ispan(-1,3,1)
    ; resxy@tmXBLabels=ispan(-1,3,1)
    ; resxy@tmXBMinorValues=fspan(-2,3,26 )
    ; resxy@tmXBMinorOn=True
    ; if fils1.ne."era/era40/msf_time_era40_in_global.nc" then 
    resxy@trXMinF=-12.0
    resxy@trXMaxF=12.0
    resxy@tmXBValues=ispan(-12,12,4)
    resxy@tmXBLabels=ispan(-12,12,4)
    resxy@tmXBMinorValues=fspan(-12,12,25)
    ; else 
    ; resxy@trXMinF=-20.0
    ; resxy@trXMaxF=20.0
    ; resxy@tmXBValues=ispan(-20,20,10)
    ; resxy@tmXBLabels=ispan(-20,20,10)
    ; resxy@tmXBMinorValues=fspan(-20,20,21)
    ; end if
    resxy@gsnRightString=""+sprintf("%2.2f",rc_decade)+"~S~o~N~/dacade"
    resxy@gsnRightStringFontHeightF=0.01
    resxy@tmXBMinorOn=True

    resxy@tmBorderThicknessF=3.0     ;调整边框粗细
    resxy@tmYUseLeft=False
    
    
    resxy@xyLineColors      = (/"blue","red"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines"/)
    resxy@xyDashPatterns    =(/0,0/)
    resxy@xyLineThicknesses=(/2.0,2.0/)
    
    resxy@vpHeightF=0.25
    resxy@vpWidthF=0.1
    plot_line=gsn_csm_xy(wks, (/width,rc_width_line/),locationxi&time, resxy)

    pres1=True
    pres2=True
    plot3= gsn_attach_plots(plot,plot_line,pres1,pres2)

    printVarSummary(wks)
draw(plot)
frame(wks)
 
end
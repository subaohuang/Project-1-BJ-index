; function level_mean(msf_time)
; begin
; msf_time=(/msf_time/10^9/)
; plev=msf_time&level
; psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
; ptop=min(plev)
; dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
; msfdp=msf_time(level|:,time|:,lon|:)*dp
; copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

; msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
; copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
; return(msf_mean)
; end 

function level_mean(msf_time)   ;;;;select 850hpa-200hpa 
begin
; msf_time=(/msf_time/10^9/)
plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
; dp=(/50/)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 

undef("search_location")
function  search_location(msf_mean_annual)
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
   if ismissing(locationdong(t)).eq.True  then 
   co_in=closest_val(0,msf_mean_annual_dong(t,:))
   locationdong(t)=msf_mean_annual_dong&lon(co_in)
   delete(co_in)
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
function draw_msf_hov(wks,LeftString,LeftString2 ,season)
begin
diri="/mnt/f/share/rea_msf/"

f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
date=cd_calendar(f1->time, 1)
; timeind=ind(date.ge.195801.and.date.le.201212)
yearstart=1979
yearend=2012
timeind=ind(date.ge.197901.and.date.le.201212)
; msf_time=f1->msf_time({850:500},:,timeind)
; msf_time=f1->msf_time({500},:,timeind)
 msf_time=f1->msf_time(:,:,timeind)

msf_time=(/msf_time/10^9/)
msf_mean=level_mean(msf_time({600:400},:,:))
; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
; msf_mean=msf_time(time|:,lon|:)
delete(msf_mean&time)
msf_mean&time=date(timeind)
printVarSummary(msf_mean)

msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
printVarSummary(msf_mean_annual)

    msf_mean_annual=month_to_annual(msf_mean, 1)
    msf_djf=msf_mean_annual
 
    if season.eq.1 then ;;express winter
    do i=0,dimsizes(msf_mean&lon)-1
    msf_djf(:,i)=(/month_to_season(msf_mean(:,i), "DJF")/)
   ;;注意点
    end do
    end if 

    if season.eq.2 then ;;express summer
    do i=0,dimsizes(msf_mean&lon)-1
    msf_mean_annual(:,i)=month_to_season(msf_mean(:,i), "JJA")
    end do
    end if 
    
    if season.eq.3 then ;;express spring
    do i=0,dimsizes(msf_mean&lon)-1
    msf_mean_annual(:,i)=month_to_season(msf_mean(:,i), "MAM")
    end do
    end if 

   
    if season.eq.4 then ;;express autumn
    do i=0,dimsizes(msf_mean&lon)-1
    msf_mean_annual(:,i)=month_to_season(msf_mean(:,i), "SON")
    end do
    end if 

    
     printVarSummary(msf_djf)
    msf_djf(0:yearend-yearstart+1-2,:)=(/msf_djf(1:,:)/)
    msf_djf(yearend- yearstart+1-1,:)=(/msf_mean({201212},:)/)
    msf_mean_annual=(/msf_djf/)


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
    rc_xi=(/rc_xi*10/)

    rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
    rc_dong_line=rc_dong*ispan(0,dimsizes(locationdong)-1,1) + rc_dong@yintercept
    rc_dong=(/rc_dong*10/)

    width=locationdong - locationxi

  

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
    rc_width_line=rc_width*ispan(0,dimsizes(locationxi)-1,1) + rc_width@yintercept
    rc_decade=(/rc_width*10/)
    resxy=True
    ; resxy@gsnRightString=""+sprintf("%2.2f",rc_decade)+"~S~o~N~/dacade"
   

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

    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
  
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame

	; res@gsnContourZeroLineThicknessF=4.0
 ;    res@gsnContourNegLineDashPattern=1
 ;    res@gsnContourPosLineDashPattern=0


    res@gsnLeftString=LeftString2
   
	res@cnLineLabelsOn =False
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
    ; res@gsnRightString="rc_xi:"+sprintf("%2.2f",rc_xi)+"~S~o~N~/dacade"+"&rc_dong:"+sprintf("%2.2f",rc_dong)+"~S~o~N~/dacade"

    ; res@cnFillPalette="testcmap"
     res@cnFillPalette  = "BlueDarkRed18"

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

	; res@tmYLValues    = ispan(2,dimsizes(timeind)/12-1,5)    
	; res@tmYLLabels    = ispan(1960,2010,5)
    res@tmYLValues    = ispan(1,dimsizes(timeind)/12-1,5)
    res@tmYLLabels=ispan(1980, 2010, 5)

	res@tmXTOn=False
	res@tmYROn=False


    ; res@tmXBValues=(/0,90,180,270/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(0,360,10)
    res@tmXBValues=ispan(90, 290, 30)
    res@tmXBLabels=(/"90~S~o~N~E","120~S~o~N~E","150~S~o~N~E","180","150~S~o~N~W","120~S~o~N~W"\
                    ,"90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(90,280,10)

	res@tmXBMinorOn=True

	res@vpHeightF=0.4
    res@vpWidthF=0.6
    plot = gsn_csm_contour(wks, msf_mean_annual(:,{90:290}), res)   ; create plot
    

    ; res@cnLevelSelectionMode = "ExplicitLevels"
    ; res@cnLinesOn=True
    ; res@cnLevels = (/0/)
    ; res@cnLineThicknessF = 4.
    ; res@cnFillOn        = False
    ; res@cnLineColor="green"

    ; plot1 = gsn_csm_contour(wks, msf_mean_annual(:,{10:180}), res)
    
    ; overlay(plot, plot1)

    ;;;;;add trend
    lineres=True
    lineres@gsLineColor="red"
    lineres@gsLineDashPattern=1
    lineres@gsLineThicknessF=3.0
    lineres@tfPolyDrawOrder="PreDraw"
    str1 = unique_string("dum1")  ; "unique_string" will return a unique
                                     ; string every time it is called from                                ;  within a single NCL session.
    str2 = unique_string("dum2")

    plot@$str1$=gsn_add_polyline(wks, plot, rc_xi_line, ispan(0,54,1), lineres)
    plot@$str2$=gsn_add_polyline(wks, plot, rc_dong_line, ispan(0,54,1), lineres)
    ;;;

    ;;add 黑线
    polyres=True
    polyres@gsLineColor="black"
    polyres@gsLineDashPattern=0
    polyres@gsLineThicknessF=4.0
    polyres@tfPolyDrawOrder="PostDraw"

    str3=unique_string("dum3")
    str4=unique_string("dum4")
    plot@$str3$=gsn_add_polyline(wks, plot, locationdong,ispan(0,54,1), polyres)
    plot@$str4$=gsn_add_polyline(wks, plot, locationxi,ispan(0,54,1), polyres)
    ;;
   


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

    ; if  LeftString.eq."MERRA" then 
    ; resxy@trXMinF=-12.0
    ; resxy@trXMaxF=12.0
    ; resxy@tmXBValues=ispan(-12,12,4)
    ; resxy@tmXBLabels=ispan(-12,12,4)
    ; resxy@tmXBMinorValues=fspan(-12,12,25)
    ; else 
    resxy@trXMinF=-15.0
    resxy@trXMaxF=15.0
    resxy@tmXBValues=ispan(-15,15,5)
    resxy@tmXBLabels=ispan(-15,15,5)
    resxy@tmXBMinorValues=fspan(-16,16,33)
    ; end if

    ; resxy@gsnRightString=""+sprintf("%2.2f",rc_decade)+"~S~o~N~/dacade"
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

    ; pres1=True
    ; pres2=True
    ; plot3= gsn_attach_plots(plot,plot_line,pres1,pres2)

    printVarSummary(wks)

 ; return[/plot,rc_xi_line,rc_dong_line/]
 return(plot)
 end 

 ;***************************************************************
 ;main program
 ;****************************************************************
  begin 
   
    ; wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个再分析资料(1979-2012)计算的msf_hov")
    
    rea= (/"20cr","era","jra55","MERRA","ncep1","ncep2"/)
 LeftString2= (/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/全球6个再分析资料(1979-2012)计算的msf_hov,没有趋势线")
    plot=new(6,graphic)
    season=1
    do i=0,5
    plot(i)=draw_msf_hov(wks,rea(i),LeftString2(i),season)   ;1 means winter 0 means aunnal 2 means summer
   end do
    pres=True
    pres@gsnMaximize=True
    pres@gsnPanelLabelBar=True
    pres@lbLabelStride=2
    pres@pmLabelBarWidthF=0.8
    pres@pmLabelBarHeightF=0.05
    pres@pmLabelBarOrthogonalPosF =0.005
    pres@gsnPanelDebug=True
    gsn_panel(wks, plot, (/3,2/), pres)
 
    end

    
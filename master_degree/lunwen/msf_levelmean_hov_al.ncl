
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
msf_mean_annual_xi=msf_mean_annual(:,{120:250})        ;;120:250加30
msf_mean_annual_dong=msf_mean_annual(:,{230:290})      ;; 
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
;;; 
undef("search_location_alt") ;;;altlantic
function  search_location_alt(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{270:330})   ;;;
msf_mean_annual_dong=msf_mean_annual(:,{320:370})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+2).gt.0))  then 
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
   if ismissing(locationxi(t)).eq.True  then 
   co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
   test=msf_mean_annual_xi(t,{290:330})
   co_in=closest_val(0,msf_mean_annual_xi(t,{290:330}))
   locationxi(t)=test&lon(co_in)
   delete(co_in)
  end if

end do
; locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

;;; 
undef("search_location_id") ;;;altlantic
function  search_location_id(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{35:55})
msf_mean_annual_dong=msf_mean_annual(:,{100:220})
msf_mean_annual_xi=(/msf_mean_annual_xi+0.5/)
msf_mean_annual_dong=(/msf_mean_annual_dong+0.5/)

;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
   co_in=closest_val(0,msf_mean_annual_xi(t,:))
   locationxi(t)=msf_mean_annual_xi&lon(co_in)
   delete(co_in)
end do

;;;judge locationdong
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(t,i+1).ge.0).and.((msf_mean_annual_dong(t,i).le.0))  then 
         a=(msf_mean_annual_dong(t,i+1)-msf_mean_annual_dong(t,i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong(t)=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a)
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
  ;  if ismissing(locationxi(t)).eq.True  then 
  ;  co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
  ;  test=msf_mean_annual_xi(t,{-70:-30})
  ;  co_in=closest_val(0,msf_mean_annual_xi(t,{-70:-30}))
  ;  locationxi(t)=test&lon(co_in)
  ;  delete(co_in)
  ; end if
end do
; locationxi=(/linmsg(locationxi, -1)/)
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
    year1=1979
    year2=2012
    timeind=ind(date.ge.197901.and.date.le.201212)

    msf_time_r=f1->msf_time(:,:,timeind)
    msf_time_f=lonPivot(msf_time_r(time|:,level|:,lon|:), 30)
    printVarSummary(msf_time_f)
    msf_time=msf_time_f(level|:,lon|:,time|:)



    msf_time=(/msf_time/10^9/)
    msf_mean=msf_time
    delete(msf_mean&time)
    msf_mean&time=date(timeind)
    printVarSummary(msf_mean)
    delete(date)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_jja=month_to_season(msf_mean(time|:,level|:,lon|:), "JJA")
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    msf_mam=month_to_season(msf_mean(time|:,level|:,lon|:), "MAM")
    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

        ; ;;注意点
        ; msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
        ; msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{201212})/)
        ; ;;;
     
    ; msf_djf_mean=level_mean(msf_djf(:,{600:400},:))

    msf_djf_mean=msf_djf(:,{500},:)
    msf_djf_mean=(/(msf_djf(:,{500},:)+msf_mam(:,{500},:))/2.0/)
    msf_mean_annual=(/msf_djf_mean(:,:)/)*1.0
    delete(msf_djf_mean)



    msf_mean_annual!0="time"
    msf_mean_annual!1="lon"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    msf_mean_annual&lon=msf_mean&lon
    
    ; msf_mean_annual=(/runave_n(msf_mean_annual, 3, 1, 0)/)

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
 
   

    width=(/dim_rmvmean_n(width,0)/)
    ; width=(/dim_standardize_n(width, 1, 0)/)
    width=(/runave(width, 7, -1)/)
    time=locationxi&time
    
    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
    rc_width_line=rc_width*ispan(0,dimsizes(locationxi)-1,1) + rc_width@yintercept
    rc_decade=(/rc_width*10/)
;;;;

;;;;alt
    location_alt=search_location_alt(msf_mean_annual)
    locationxi_alt=location_alt[0]
    locationdong_alt=location_alt[1]
    rc_xi_alt=regCoef_n(ispan(0,dimsizes(locationxi_alt)-1,1),locationxi_alt , 0, 0)
    rc_xi_alt_line=rc_xi_alt*ispan(0,dimsizes(locationxi_alt)-1,1) + rc_xi_alt@yintercept
    rc_xi_alt=(/rc_xi_alt*10/)

    rc_dong_alt=regCoef_n(ispan(0,dimsizes(locationdong_alt)-1,1),locationdong_alt, 0, 0)
    rc_dong_alt_line=rc_dong_alt*ispan(0,dimsizes(locationdong_alt)-1,1) + rc_dong_alt@yintercept
    rc_dong_alt=(/rc_dong_alt*10/)

    width_alt=locationdong_alt - locationxi_alt

    rc_width_alt=regCoef_n(ispan(0,dimsizes(locationxi_alt)-1,1),width_alt , 0, 0)
    rc_width_alt_line=rc_width_alt*ispan(0,dimsizes(locationxi_alt)-1,1) + rc_width_alt@yintercept
    rc_decade_alt=(/rc_width_alt*10/)
    print("alt的增宽速度="+rc_decade_alt)
;;;

    location_id=search_location_id(msf_mean_annual)
    locationxi_id=location_id[0]
    locationdong_id=location_id[1]



     

;***********************************************
;draw hov
;***********************************************

    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
  
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame

    res@gsnLeftString=LeftString2
    res@gsnLeftStringFontHeightF=30.0

    res@cnLineLabelsOn =False
    res@cnFillDrawOrder="PreDraw"
    res@cnLinesOn=False

    res@cnFillPalette  = "BlueDarkRed18"
    ; res@cnFillPalette="MPL_jet"

    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = -7.0               ; min level
    res@cnMaxLevelValF       =  7.0               ; max leve
    res@cnLevelSpacingF      =0.5
    res@lbOrientation         ="Vertical"

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXBMode       ="Explicit"
    res@tmYLLabelsOn=True
    res@tmYUseLeft=True
    res@tmYLValues    = ispan(1,dimsizes(timeind)/12-1,5)
    res@tmYLLabels=ispan(1980, 2010, 5)

    res@tmXTOn=False
    res@tmYROn=False


    ; res@tmXBValues=ispan(30, 390, 60)
    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    ; res@tmXBLabels=(/"30E","90E","150E","150W"\
    ;                 ,"90W","30W","30E"/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
   
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)
    res@tmYLMinorOn=True
    res@tmYLMinorValues=ispan(1,34,1)

     
    res@tmXBLabelFontHeightF=0.02
    res@tmXBLabelFontThicknessF=5.0
    res@tmYLLabelFontThicknessF=5.0
    res@tmYLLabelFontHeightF=0.02
    res@tmYLMajorThicknessF=4.0
    res@tmXBMajorThicknessF=4.0
    res@tmBorderThicknessF=4.0


    res@vpHeightF=0.4
    res@vpWidthF=0.7
    plot = gsn_csm_contour(wks, msf_mean_annual(:,:), res)   ; create plot

    ;;;;;add trend
    lineres=True
    lineres@gsLineColor="red"
    lineres@gsLineDashPattern=1
    lineres@gsLineThicknessF=3.0
    lineres@tfPolyDrawOrder="PreDraw"
    str1 = unique_string("dum1")  ; "unique_string" will return a unique
                                     ; string every time it is called from                                ;  within a single NCL session.
    str2 = unique_string("dum2")
    
    str11 = unique_string("dum1")  ; "unique_string" will return a unique
                                     ; string every time it is called from                                ;  within a single NCL session.
    str22 = unique_string("dum2")

    plot@$str1$=gsn_add_polyline(wks, plot, rc_xi_line, ispan(0,34,1), lineres)
    plot@$str2$=gsn_add_polyline(wks, plot, rc_dong_line, ispan(0,34,1), lineres)
    plot@$str11$=gsn_add_polyline(wks, plot, rc_xi_alt_line, ispan(0,34,1), lineres)
    plot@$str22$=gsn_add_polyline(wks, plot, rc_dong_alt_line, ispan(0,34,1), lineres)
    ;;;

    ;;add 黑线
    polyres=True
    polyres@gsLineColor="black"
    polyres@gsLineDashPattern=0
    polyres@gsLineThicknessF=4.0
    polyres@tfPolyDrawOrder="PostDraw"

    str3=unique_string("dum3")
    str4=unique_string("dum4")
    str5=unique_string("dum5")
    str6=unique_string("dum6")
    str7=unique_string("dum7")
     str8=unique_string("dum8")

    plot@$str3$=gsn_add_polyline(wks, plot, locationdong,ispan(0,34,1), polyres)
    plot@$str4$=gsn_add_polyline(wks, plot, locationxi,ispan(0,34,1), polyres)

    plot@$str5$=gsn_add_polyline(wks, plot, locationdong_alt,ispan(0,34,1), polyres)
    plot@$str6$=gsn_add_polyline(wks, plot, locationxi_alt,ispan(0,34,1), polyres)
    plot@$str7$=gsn_add_polyline(wks, plot, locationxi_id,ispan(0,34,1), polyres)
    ; plot@$str8$=gsn_add_polyline(wks, plot, locationdong_id,ispan(0,34,1), polyres)
    ;;
   
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

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/全球6个再分析资料(1979-2012)计算的msf_hov_冬季")
    plot=new(6,graphic)
    season=1
    do i=0,5
    plot(i)=draw_msf_hov(wks,rea(i),LeftString2(i),season)   ;1 means winter 0 means aunnal 2 means summer
   end do
    pres=True
    pres@gsnMaximize=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.8
    pres@pmLabelBarHeightF=0.05
    pres@pmLabelBarOrthogonalPosF =0.005
    pres@gsnPanelDebug=True
    gsn_panel(wks, plot, (/3,2/), pres)
  ; draw(plot(0))
  ; frame(wks)
    end

    
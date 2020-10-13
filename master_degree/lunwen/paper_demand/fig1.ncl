
function  search_location(msf_mean_annual)
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:300})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+8).gt.0))  then 
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
; locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)

 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;;;;二维的位置查找
undef("search_location_2d")
function  search_location_2d(msf_mean_annual)
begin
msf_mean_annual_xi=msf_mean_annual({120:250})
msf_mean_annual_dong=msf_mean_annual({230:290})
locationxi=msf_mean_annual@_FillValue
locationdong=locationxi
;;;judge locationxi
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(i+1).ge.0).and.((msf_mean_annual_xi(i).le.0).and.(msf_mean_annual_xi(i+8).gt.0))  then 
         a=(msf_mean_annual_xi(i+1)-msf_mean_annual_xi(i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi=msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(i)/a
        ; locationxi=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if  (msf_mean_annual_xi(i).eq.0)  then 
        locationxi=msf_mean_annual_dong&xi(i)
        end if  
        if ismissing(locationxi.ne.True)  then 
        break
        end if
     end do
;;;judge locationdong

    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(i+1).le.0).and.((msf_mean_annual_dong(i).ge.0))  then 
         a=(msf_mean_annual_dong(i+1)-msf_mean_annual_dong(i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong=msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(i)/a    
        ; locationdong=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
         if  (msf_mean_annual_dong(i).eq.0)  then 
        locationdong=msf_mean_annual_dong&lon(i)
        end if    
        if ismissing(locationdong.ne.True)  then 
        break
        end if
     end do
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

;;;;二维的位置查找
undef("search_location_2d_alt")
function  search_location_2d_alt(msf_mean_annual)
begin
msf_mean_annual_xi=msf_mean_annual({270:330})
msf_mean_annual_dong=msf_mean_annual({320:370})
locationxi=msf_mean_annual@_FillValue
locationdong=locationxi
;;;judge locationxi
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(i+1).ge.0).and.((msf_mean_annual_xi(i).le.0).and.(msf_mean_annual_xi(i+8).gt.0))  then 
         a=(msf_mean_annual_xi(i+1)-msf_mean_annual_xi(i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi=msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(i)/a
        ; locationxi=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if  (msf_mean_annual_xi(i).eq.0)  then 
        locationxi=msf_mean_annual_dong&xi(i)
        end if  
        if ismissing(locationxi.ne.True)  then 
        break
        end if
     end do
;;;judge locationdong

    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(i+1).le.0).and.((msf_mean_annual_dong(i).ge.0))  then 
         a=(msf_mean_annual_dong(i+1)-msf_mean_annual_dong(i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong=msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(i)/a    
        ; locationdong=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
         if  (msf_mean_annual_dong(i).eq.0)  then 
        locationdong=msf_mean_annual_dong&lon(i)
        end if    
        if ismissing(locationdong.ne.True)  then 
        break
        end if
     end do
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;;;
undef("level_mean")
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
;************************************************
; 画20cr,ncep,era的msf气候态
;************************************************
function draw_all_msf(LeftString,wks,title)
begin
    diri="/mnt/f/share/rea_msf/"
    fils1=systemfunc("ls "+diri+"msf_time_"+LeftString+"_global.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    level=f1->level
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind1=ind(f1->level.eq.400)
    levind2=ind(f1->level.eq.600)
    levind3=ind(level.ge.400.and.f1->level.le.600)
    levind=ind(f1->level.eq.500)
    msf_time_20cr_r=f1->msf_time(:,:,timeind)
    ud_time_20cr_r=f1->ud_time(:,:,timeind)
    w_time_20cr_r=f1->w_time(:,:,timeind)
    
    msf_time_20cr_f=lonPivot(msf_time_20cr_r(time|:,level|:,lon|:), 30)
    ud_time_20cr_f=lonPivot(ud_time_20cr_r(time|:,level|:,lon|:), 30)
    w_time_20cr_f=lonPivot(w_time_20cr_r(time|:,level|:,lon|:), 30)
    
    msf_time_20cr=msf_time_20cr_f(level|:,lon|:,time|:)
    ud_time_20cr=ud_time_20cr_f(level|:,lon|:,time|:)
    w_time_20cr=w_time_20cr_f(level|:,lon|:,time|:)


    msf_time_20cr=(/msf_time_20cr/10^9/)
    
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)

    printVarSummary(ud_time_20cr)
    printVarSummary(msf_time_20cr)
   
    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)
;;;
; if LeftString.eq."MERRA"  then
;   msf_clm_20cr=dim_avg_n_Wrap(msf_time_20cr(time|:,level|:,lon|:), 0)
; end if
;;;选择冬季


msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
ud_clm_20cr=dim_avg_n_Wrap(month_to_season(ud_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
w_clm_20cr=dim_avg_n_Wrap(month_to_season(w_time_20cr(time|:,level|:,lon|:),"DJF"), 0)

; msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"JJA"), 0)
; ud_clm_20cr=dim_avg_n_Wrap(month_to_season(ud_time_20cr(time|:,level|:,lon|:),"JJA"), 0)
; w_clm_20cr=dim_avg_n_Wrap(month_to_season(w_time_20cr(time|:,level|:,lon|:),"JJA"), 0)

; msf_clm_20cr=dim_avg_n_Wrap(msf_time_20cr(time|:,level|:,lon|:), 0)
; ud_clm_20cr=dim_avg_n_Wrap(ud_time_20cr(time|:,level|:,lon|:), 0)
; w_clm_20cr=dim_avg_n_Wrap(w_time_20cr(time|:,level|:,lon|:), 0)



msf_clm_20cr_levelmean=dim_avg_n_Wrap(msf_clm_20cr(level|:,lon|:),0)

intension_pa=dim_avg_n_Wrap(msf_clm_20cr_levelmean({150:210}), 0)
intension_id=dim_avg_n_Wrap(msf_clm_20cr_levelmean({70:100}), 0)
intension_al=dim_avg_n_Wrap(msf_clm_20cr_levelmean({320:340}), 0)

print("********************"+LeftString+"*********************")
print("intension_id="+intension_id)
print("intension_pa="+intension_pa)
print("intension_al="+intension_al)




;;;
;;;location

;     location_r=search_location(msf_time_20cr(time|:,level|levind1,lon|:))
;     locationxi_r=location_r[0]
;     locationdong_r=location_r[1]

;     locationxi=new(12,typeof(location_r[0]))
;     locationdong=new(12,typeof(location_r[0]))

;     do i=0,11
;     locationxi(i)=dim_avg_n(locationxi_r(i:dimsizes(locationxi_r)-1:12),0)
;     locationdong(i)=dim_avg_n(locationdong_r(i:dimsizes(locationxi_r)-1:12),0)
;     end do  

;     location_r1=search_location(msf_time_20cr(time|:,level|levind2,lon|:))
;     locationxi_r1=location_r1[0]
;     locationdong_r1=location_r1[1]

;     locationxi1=new(12,typeof(location_r1[0]))
;     locationdong1=new(12,typeof(location_r1[0]))

;     do i=0,11
;     locationxi1(i)=dim_avg_n(locationxi_r1(i:dimsizes(locationxi_r1)-1:12),0)
;     locationdong1(i)=dim_avg_n(locationdong_r1(i:dimsizes(locationxi_r1)-1:12),0)
;     end do  


; xi=(locationxi(10)+locationxi(11)+locationxi(1))/3.0
; dong=(locationdong(10)+locationdong(11)+locationdong(1))/3.0

; xi1=(locationxi1(10)+locationxi1(11)+locationxi1(1))/3.0
; dong1=(locationdong1(10)+locationdong1(11)+locationdong1(1))/3.0


  location_pa=search_location_2d(msf_clm_20cr(level|levind,lon|:))
   xi_pa=location_pa[0]
   dong_pa=location_pa[1]

  location_alt=search_location_2d_alt(msf_clm_20cr(level|levind,lon|:))
   xi_alt=location_alt[0]
   dong_alt=location_alt[1]
   
   width_alt=dong_alt - xi_alt
   dong_id=xi_pa
   xi_id=40.0
; location2=search_location_2d(msf_clm_20cr_levelmean)
; width=location2[1]-location2[0]
print("**********"+title+"****************")
print("locationxi="+xi_alt)
print("locationdong="+dong_alt)
print("width="+width_alt)


;************************************************
; 画图气候态
;************************************************

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    ; res@gsnLeftString=title
   
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    ; res@cnFillPalette  = "BlueDarkRed18"
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 0.5
    res@cnFillOn = True 
    ; res@gsnYRefLine=500
    ; res@gsnYRefLineColor="black"
    ; res@gsnContourNegLineDashPattern=2
    ; res@gsnContourZeroLineThicknessF=4.0

    res@cnLineLabelsOn=False
    res@cnLineLabelInterval=2
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False
    res@cnLineThicknessF=2.0

    ; res@cnExplicitLineLabelsOn=True
    ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
    ; res@cnLabelScalingMode="ScaleFactor"
    ; res@cnLabelScaleValueF=10^9

    res@lbLabelBarOn=False
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"

    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    ; res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
    ;                 ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
    
      res@tmXBLabels=(/"30~S~o","90~S~o","150~S~o","150~S~o"\
                    ,"90~S~o","30~S~o","30~S~o~N~E"/)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)
     
    res@tiYAxisString=""
    
    res@tmXBLabelFontHeightF=0.05
    res@tmYLLabelFontHeightF=0.02
    res@tmYLMajorThicknessF=1.0
    res@tmXBMajorThicknessF=4.0

    res@tmBorderThicknessF=4.0


    res@vcRefMagnitudeF = 1.5                ; define vector ref mag
    res@vcRefLengthF    = 0.03        ; define length of vec ref
    res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.03       ; thin out vectors
    res@vcMapDirection  = False

    res@vcLineArrowThicknessF=2.0
    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15
    
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    if title.eq."(a)20CR".or.title.eq."(c)JRA55" then
    res@vcRefAnnoOrthogonalPosF  = -1.02
    end if
   if title.eq."(b)ERAIM".or.title.eq."(d)MERRA" then
    res@vcRefAnnoOrthogonalPosF  = -1.
    end if

 if title.eq."(e)NCEP1".or.title.eq."(f)NCEP2" then
    res@vcRefAnnoOrthogonalPosF  = -1.16
    end if
    
    



    res@vcRefAnnoPerimOn=False
    res@vcRefAnnoString2On=False
 

    res@vpHeightF=0.3
    res@vpWidthF=0.7

    
    no_ba=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA"/)
    no_ya=(/"(b)ERAIM","(d)MERRA","(f)NCEP2"/)

    if any(no_ba.eq.title) then
    res@tmXBOn=False
    end if

    if any(no_ya.eq.title) then
    res@tmYLOn=False
    res@tmYUseLeft=False
    end if
    
    plot=gsn_csm_pres_hgt_vector(wks,msf_clm_20cr(:,:),ud_clm_20cr(:,:),w_clm_20cr(:,:),res)
  ; plot=gsn_csm_pres_hgt(wks, msf_clm_20cr(:,：), res)
  
    res@cnLevelSelectionMode ="ExplicitLevels"
    res@cnLevels=(/-0.5/)

    res@cnFillOn=False
    res@cnMonoLineColor=True
    res@cnLineColor="blue"
    res@cnLineDashPattern=0
    res@cnLineThicknessF =4.0

   plot1=gsn_csm_pres_hgt(wks,msf_clm_20cr(:,{31:180}),res)
    res@cnLineColor="black"
    res@cnLineDashPattern=0
     res@cnLevels=(/0/)
   plot2=gsn_csm_pres_hgt(wks,msf_clm_20cr(:,{120:}),res)

   overlay(plot, plot1)
   overlay(plot, plot2)

 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 30, 25)+360.0
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon6=fspan(30, 39, 25)

      ploylon4=fspan(xi_pa, dong_pa, 100)
      ploylon5=fspan(xi_alt, dong_alt, 100)
      ploylon7=fspan(xi_id, dong_id, 100)

      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel6=conform_dims(dimsizes(ploylon6), ploylevel, 0)
      
      ploylevel4=conform_dims(dimsizes(ploylon4), 500, 0)
      ploylevel5=conform_dims(dimsizes(ploylon5), 500, 0)
      ploylevel7=conform_dims(dimsizes(ploylon7), 500, 0)



    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str31 = unique_string("dum31")
    str4 = unique_string("dum4")
    str5 = unique_string("dum5")
    str51 = unique_string("dum51")
    str6 = unique_string("dum6")
    str7 = unique_string("dum7")
    str8 = unique_string("dum8")
    str9 = unique_string("dum9")     
    str10 = unique_string("dum10")
    str11 = unique_string("dum11") 
    ployres@gsLineDashPattern=0
     
    plot@$str1$= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)

    ployres@gsLineDashPattern=16

    plot@$str2$= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)

    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str31$= gsn_add_polyline(wks,plot,  ploylon6,ploylevel6,  ployres)

     ployres@gsLineDashPattern=0
     ployres@gsLineThicknessF=4.0
     plot@$str4$= gsn_add_polyline(wks,plot,  ploylon4,ploylevel4,  ployres)
     plot@$str5$= gsn_add_polyline(wks,plot,  ploylon5,ploylevel5,  ployres)
     plot@$str51$= gsn_add_polyline(wks,plot,  ploylon7,ploylevel7,  ployres)

     ploymarker=True
     ploymarker@gsMarkerColor="black"
     ploymarker@gsMarkerThicknessF=1.5
     ploymarker@gsMarkerIndex=10


     plot@$str6$= gsn_add_polymarker(wks,plot,xi_pa+4,500,ploymarker)
     ploymarker@gsMarkerIndex=11     
     plot@$str7$= gsn_add_polymarker(wks,plot,dong_pa-4,500,ploymarker)

     ploymarker@gsMarkerIndex=10 
     plot@$str8$= gsn_add_polymarker(wks,plot,xi_alt+4,500,ploymarker)
     ploymarker@gsMarkerIndex=11     
     plot@$str9$= gsn_add_polymarker(wks,plot,dong_alt-4,500,ploymarker)
     
    ploymarker@gsMarkerIndex=10 
     plot@$str10$= gsn_add_polymarker(wks,plot,xi_id+4,500,ploymarker)
     ploymarker@gsMarkerIndex=11     
     plot@$str11$= gsn_add_polymarker(wks,plot,dong_id-4,500,ploymarker)

delete(res)
 return (plot)
end 
;;;;
;;主程序
;;;

begin
  rea= (/"20cr","era","jra55","MERRA","ncep1","ncep2"/)
    wks_type = "pdf"
    ; wks_type@wkWidth = 1800
    ; wks_type@wkHeight = 1600
    ; wks_type@wkWidth = 8000
    ; wks_type@wkHeight = 8000
    wks=gsn_open_wks(wks_type, "/mnt/f/share/era/data/再分析资料DJF平均")

   plot=new(6,graphic)
   title=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)
   do i=0,5
    plot(i)=draw_all_msf(rea(i),wks,title(i))
  end do



    pres=True
    pres@gsnPanelLabelBar=True
    pres@gsnMaximize=True
    pres@pmLabelBarWidthF=0.6
    pres@pmLabelBarHeightF=0.07
    ; pres@pmLabelBarOrthogonalPosF=1.2
    pres@gsnPanelFigureStrings=(/"a","b","c","d","e","f"/)
    pres@gsnPanelFigureStringsJust="TopLeft"
    pres@gsnPanelFigureStringsBackgroundFillColor="Transparent"
    pres@gsnPanelFigureStringsPerimOn=False
    pres@gsnPanelXWhiteSpacePercent=0.0001
    pres@gsnPanelRight=0.98


    ; pres@lbTitleString="10~S~9~N~Kg s~S~-1"
    pres@lbTitleFontThicknessF=5.0
    pres@lbTitlePosition  = "Bottom"              ; title position
    pres@lbTitleFontHeightF= .02               ; make title smaller
    pres@gsnPanelDebug=True
             ; title direction
    gsn_panel(wks, plot, (/3,2/), pres)
    
 end
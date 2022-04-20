
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
; locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)

 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
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
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind1=ind(f1->level.eq.400)
    levind2=ind(f1->level.eq.600)
    msf_time_20cr=f1->msf_time(:,:,timeind)
    ud_time_20cr=f1->ud_time(:,:,timeind)
    w_time_20cr=f1->w_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)

    printVarSummary(ud_time_20cr)
    printVarSummary(msf_time_20cr)

    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)

;;;选择冬季
    msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
    ud_clm_20cr=dim_avg_n_Wrap(month_to_season(ud_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
    w_clm_20cr=dim_avg_n_Wrap(month_to_season(w_time_20cr(time|:,level|:,lon|:),"DJF"), 0)

    ;;;
; if LeftString.eq."mre"  then
  ; msf_clm_20cr=dim_avg_n_Wrap(msf_time_20cr(time|:,level|:,lon|:), 0)
  ;   ud_clm_20cr=dim_avg_n_Wrap(ud_time_20cr(time|:,level|:,lon|:), 0)
  ;   w_clm_20cr=dim_avg_n_Wrap(w_time_20cr(time|:,level|:,lon|:), 0)
; end if
;;;
;;;location

    location_r=search_location(msf_time_20cr(time|:,level|levind1,lon|:))
    locationxi_r=location_r[0]
    locationdong_r=location_r[1]

    locationxi=new(12,typeof(location_r[0]))
    locationdong=new(12,typeof(location_r[0]))

    do i=0,11
    locationxi(i)=dim_avg_n(locationxi_r(i:dimsizes(locationxi_r)-1:12),0)
    locationdong(i)=dim_avg_n(locationdong_r(i:dimsizes(locationxi_r)-1:12),0)
    end do  

    location_r1=search_location(msf_time_20cr(time|:,level|levind2,lon|:))
    locationxi_r1=location_r1[0]
    locationdong_r1=location_r1[1]

    locationxi1=new(12,typeof(location_r1[0]))
    locationdong1=new(12,typeof(location_r1[0]))

    do i=0,11
    locationxi1(i)=dim_avg_n(locationxi_r1(i:dimsizes(locationxi_r1)-1:12),0)
    locationdong1(i)=dim_avg_n(locationdong_r1(i:dimsizes(locationxi_r1)-1:12),0)
    end do  


xi=(locationxi(10)+locationxi(11)+locationxi(1))/3.0
dong=(locationdong(10)+locationdong(11)+locationdong(1))/3.0

xi1=(locationxi1(10)+locationxi1(11)+locationxi1(1))/3.0
dong1=(locationdong1(10)+locationdong1(11)+locationdong1(1))/3.0


;************************************************
; 画图气候态
;************************************************

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString=title
    res@gsnRightString="10~S~9~N~Kg s~S~-1"
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1.0
    res@cnFillOn = True 
    ; res@gsnYRefLine=500
    ; res@gsnYRefLineColor="black"

    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=2.5

    res@cnLineLabelsOn=False
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

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
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
     
    res@tiYAxisString=""

    res@vcRefMagnitudeF = 2.0                ; define vector ref mag
    res@vcRefLengthF    = 0.03        ; define length of vec ref
    res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.025       ; thin out vectors
    res@vcMapDirection  = False

    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15

    res@vcRefAnnoOrthogonalPosF  = -1.14

    res@vcRefAnnoPerimOn=False
    res@vcRefAnnoString2On=False


    res@vpHeightF=0.3
    res@vpWidthF=0.7

   plot=gsn_csm_pres_hgt_vector(wks,msf_clm_20cr,ud_clm_20cr,w_clm_20cr,res)
 
 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 39, 25)
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon4=fspan(xi, dong, 100)
      ploylon5=fspan(xi1, dong1, 100)

      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      
      ploylevel4=conform_dims(dimsizes(ploylon4), 400, 0)
      ploylevel5=conform_dims(dimsizes(ploylon5), 600, 0)



    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str4 = unique_string("dum4")
    str5 = unique_string("dum5")
    str6 = unique_string("dum6")
    str7 = unique_string("dum7")
    str8 = unique_string("dum8")
    str9 = unique_string("dum9")     
    ployres@gsLineDashPattern=0
     
    plot@$str1$= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)

    ployres@gsLineDashPattern=16

    plot@$str2$= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)

    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)
    
     ployres@gsLineDashPattern=0
     ployres@gsLineThicknessF=2.0
     plot@$str4$= gsn_add_polyline(wks,plot,  ploylon4,ploylevel4,  ployres)
      plot@$str5$= gsn_add_polyline(wks,plot,  ploylon5,ploylevel5,  ployres)

     ploymarker=True
     ploymarker@gsMarkerColor="black"
     ploymarker@gsMarkerThicknessF=1.0
     ploymarker@gsMarkerIndex=10

     plot@$str6$= gsn_add_polymarker(wks,plot,xi+10,400,ploymarker)
     ploymarker@gsMarkerIndex=11     
     plot@$str7$= gsn_add_polymarker(wks,plot,dong-10,400,ploymarker)
     
     plot@$str8$= gsn_add_polymarker(wks,plot,xi1+10,600,ploymarker)
     ploymarker@gsMarkerIndex=11     
     plot@$str9$= gsn_add_polymarker(wks,plot,dong1-10,600,ploymarker)

 return (plot)
end 
;;;;
;;主程序
;;;

begin
  rea= (/"20cr","era","jra55","MERRA","ncep1","mre"/)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料DJF平均")

   plot=new(6,graphic)
   ; title=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)
   title=(/"(a)20CR","(","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)
   do i=0,5
    plot(i)=draw_all_msf(rea(i),wks,title(i))
  end do



    pres=True
    pres@gsnPanelLabelBar=True
    pres@gsnMaximize=False

    pres@pmLabelBarWidthF=0.6
    pres@pmLabelBarHeightF=0.07
    pres@lbTitleString="10~S~9~N~Kg s~S~-1"
    pres@lbTitlePosition  = "Bottom"              ; title position
    pres@lbTitleFontHeightF= .008               ; make title smaller
             ; title direction
    ; gsn_panel(wks, plot, (/3,2/), pres)

    draw(plot(1))
    frame(wks)
    
 end
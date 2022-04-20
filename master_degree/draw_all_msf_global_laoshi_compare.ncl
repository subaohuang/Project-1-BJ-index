
;************************************************
; 画20cr,ncep,era的msf气候态
;************************************************
begin
    diri="/mnt/f/share/rea_msf/"
    fils1=systemfunc("ls "+diri+"msf_time_mre_global.nc")
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

    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)

    delete(w_time_20cr&time)
    w_time_20cr&time=date(timeind)

    delete(ud_time_20cr&time)
    ud_time_20cr&time=date(timeind)


    msf_djf=month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF")
    ud_djf=month_to_season(ud_time_20cr(time|:,level|:,lon|:),"DJF")
    w_djf=month_to_season(w_time_20cr(time|:,level|:,lon|:),"DJF")

    msf_djf&time=ispan(1979,2012,1)  
    ud_djf&time=ispan(1979,2012,1)  
    w_djf&time=ispan(1979,2012,1)  
 
   
    msf_djf_1=dim_avg_n_Wrap(msf_djf({1979:1994},:,:), 0)
    msf_djf_2=dim_avg_n_Wrap(msf_djf({1995:2012},:,:),0)

    w_djf_1=dim_avg_n_Wrap(w_djf({1979:1994},:,:), 0)
    w_djf_2=dim_avg_n_Wrap(w_djf({1995:2012},:,:),0)

    ud_djf_1=dim_avg_n_Wrap(ud_djf({1979:1994},:,:), 0)
    ud_djf_2=dim_avg_n_Wrap(ud_djf({1995:2012},:,:),0)
   
    msf_diff=msf_djf_1
    w_diff=w_djf_1
    ud_diff=ud_djf_1

    msf_diff=(/msf_djf_2- msf_djf_1/)
    ud_diff=(/ud_djf_2- ud_djf_1/)
    w_diff=(/w_djf_2- w_djf_1/)

;************************************************
; 画图气候态
;************************************************
 wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/ENS的MSF两个时间段的平均")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnRightString="10~S~9~N~Kg s~S~-1"
    res@gsnMaximize=False
    res@gsnLeftString="1979-1994 mean"
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
    plot=new(3,graphic)

    plot(0)=gsn_csm_pres_hgt_vector(wks,msf_djf_1,ud_djf_1,w_djf_1,res)
    res@gsnLeftString="1995-2012 mean"
    plot(1)=gsn_csm_pres_hgt_vector(wks,msf_djf_2,ud_djf_2,w_djf_2,res)

    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -3
    res@cnMaxLevelValF = 3
    res@cnLevelSpacingF = 0.5
    
    
    res@vcRefMagnitudeF = 1.0                ; define vector ref mag
      res@gsnLeftString="(b)-(a)"
    plot(2)=gsn_csm_pres_hgt_vector(wks,msf_diff,ud_diff,w_diff,res)
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


      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)




    ployres@gsLineDashPattern=0
     
    do i=0,1
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")

    plot@$str1$= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)

    ployres@gsLineDashPattern=16

    plot@$str2$= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel2,  ployres)

    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot(i),  ploylon3,ploylevel3,  ployres)
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
    gsn_panel(wks, plot, (/3,1/), pres)
    
 end
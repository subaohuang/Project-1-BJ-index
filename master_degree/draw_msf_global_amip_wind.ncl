
;;****************************************************************
;;main program
;;*************************************************************

begin
 ; diri="/mnt/f/share/cmip/history/"
 ; model="CanESM2"
  model="CanAM4"
 ; model="GFDL-CM3"
 ; diri="/mnt/f/share/CMIP5/msf_history_rcp45/"
 diri="/mnt/f/share/CMIP5/msf_AMIP/"

 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")
 f1=addfile(fils1, "r")
 date=cd_calendar(f1->time, 1)
 
 yearstart=197901
 yearend=200812

 timeind=ind(date.le.yearend.and.date.ge.yearstart)
 msf_time=f1->msf_time(:,:,timeind)
 msf_time=(/msf_time/10^9/)
 msf_time&level@units="hpa"
 msf_time&time=ispan(0, dimsizes(timeind)-1, 1)

 msf_clm=dim_avg_n_Wrap(msf_time, 2)

    ud_time=f1->ud_time(:,:,timeind)
    w_time=f1->w_time(:,:,timeind)
    
    ud_clm=dim_avg_n_Wrap(ud_time, 2)
    w_clm=dim_avg_n_Wrap(w_time, 2)

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/"+model+"amip_"+yearstart+"-"+yearend+"")
;********************************************************
;;draw
;********************************************************

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="20CR"
    ; res@gsnRightString="10~S~9~N~Kg s~S~-1"
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1.0
    res@cnFillOn = True 
     
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
    res@vcRefLengthF    = 0.04        ; define length of vec ref
    res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.02       ; thin out vectors
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
    plot=new(1,graphic)
    plot1=gsn_csm_pres_hgt_vector(wks,msf_clm,ud_clm,w_clm,res)
    
;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 8.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 10)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
 
    dum1=new(1,graphic)
    dum2=dum1
    dum3=dum1

     i=0
    ployres@gsLineDashPattern=0
    dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
     

        overlay(plot,plot1)
        overlay(plot,iplot)
        draw(plot)
        frame(wks)
  end 
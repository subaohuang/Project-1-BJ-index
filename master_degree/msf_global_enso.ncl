;*******************************************
;enino时期以及1994年IOD事件全纬度的walker环流变化情况
;*********************************************
begin
 diri="/mnt/f/share/era/data/"
 fils1=systemfunc("ls "+diri+"msf_time_era_global.nc")
 f1=addfile(fils1, "r")
 msf_time=f1->msf_time(:,:,:)
 ud_time=f1->ud_time(:,:,:)
 w_time=f1->w_time
 level=f1->level(:)
 lon=f1->lon(:)

 msf_time&level@units="hpa"
 copy_VarCoords(msf_time, ud_time)
 copy_VarCoords(msf_time, w_time)
 msf_time=(/msf_time/10^9/)
 printVarSummary(msf_time)

;;;读入enino指数
    nino3_4=asciiread(diri+"era_nino3_4.ascii", 456, "float")  ;;已经去趋势了
    nino3_4_nor=dim_standardize_n(nino3_4, 1, 0)

    enino=ind(nino3_4_nor.gt.1.0)
    lanila=ind(nino3_4_nor.lt.-1.0)
    mid=ind(nino3_4_nor.le.1.0.and.nino3_4_nor.ge.-1.0)
;;;

;;;选择msf，ENSO成分
    msf_enino=dim_avg_n_Wrap(msf_time(:,:,enino),2)
    msf_lanina=dim_avg_n_Wrap(msf_time(:,:,lanila),2)
    msf_mid=dim_avg_n_Wrap(msf_time(:,:,mid),2)

    w_enino=dim_avg_n_Wrap(w_time(:,:,enino),2)
    w_lanina=dim_avg_n_Wrap( w_time(:,:,lanila),2)
    w_mid=dim_avg_n_Wrap( w_time(:,:,mid),2)

    ud_enino=dim_avg_n_Wrap(ud_time(:,:,enino),2)
    ud_lanina=dim_avg_n_Wrap( ud_time(:,:,lanila),2)
    ud_mid=dim_avg_n_Wrap(ud_time(:,:,mid),2)

;;;选择1994年IOD事件，选择秋季（SON）
    msf_iod=dim_avg_n_Wrap(msf_time(:,:,{188:190}),2)
     
    w_iod=dim_avg_n_Wrap(w_time(:,:,{188:190}),2)

    ud_iod=dim_avg_n_Wrap(ud_time(:,:,{188:190}),2)

; ;;;1997年极端enino事件
    ; msf_str_enino=dim_avg_n_Wrap(msf_time(:,:,{228:239}),2)
         
    ; w_str_enino=dim_avg_n_Wrap(w_time(:,:,{228:239}),2)

    ; ud_str_enino=dim_avg_n_Wrap(ud_time(:,:,{228:239}),2)
;;

;;; 1983年
    msf_str_enino=dim_avg_n_Wrap(msf_time(:,:,{48:59}),2)
         
    w_str_enino=dim_avg_n_Wrap(w_time(:,:,{48:59}),2)

    ud_str_enino=dim_avg_n_Wrap(ud_time(:,:,{48:59}),2)


;;;---------
 plot=new(3,graphic)
 
 
wks=gsn_open_wks("pdf", ""+diri+"msf_global_enso")
wks1=gsn_open_wks("pdf", ""+diri+"msf_global_iod1994&enino1997")
wks3=gsn_open_wks("pdf", ""+diri+"msf_global_1995")
;;;画图---------------------------------------
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="enino"
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    ; res@cnFillPalette  = "BlueDarkRed18"
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -4
    res@cnMaxLevelValF = 4
    res@cnLevelSpacingF = 0.5
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
    res@cnLinePalette="BlueDarkRed18"
    res@cnLineLabelsOn=True
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

    ; res@cnExplicitLineLabelsOn=True
    ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
    ; res@cnLabelScalingMode="ScaleFactor"
    ; res@cnLabelScaleValueF=10^9

    res@lbLabelBarOn=False
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)

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

    res@vcRefAnnoOrthogonalPosF  = -1.05
    res@vcRefAnnoPerimOn=False

     res@vpHeightF=0.3
     res@vpWidthF=0.7
    
   plot(0)=gsn_csm_pres_hgt_vector(wks,w_enino,ud_enino,w_enino,res)
   res@gsnRightString="lanina"
   plot(1)=gsn_csm_pres_hgt_vector(wks,w_lanina,ud_lanina,w_lanina,res)
   res@gsnRightString="normal"
   plot(2)=gsn_csm_pres_hgt_vector(wks,w_mid,ud_mid,w_mid,res)

   ;;;IOD1994年事件
   res@gsnRightString="1994 iod modoki"
   plot_iod=gsn_csm_pres_hgt_vector(wks1,w_iod,ud_iod,w_iod,res) 
  
   ;;;1997年极端enino事件
   delete(res@gsnRightString)
    res@gsnRightString="1997enino"
   plot_enino=gsn_csm_pres_hgt_vector(wks1,w_str_enino,ud_str_enino,w_str_enino,res) 
   
   ;;; 1995年
   plot_1995=gsn_csm_pres_hgt_vector(wks3,dim_avg_n_Wrap(msf_time(:,:,192:203),2),dim_avg_n_Wrap(ud_time(:,:,192:203),2), \
                                     dim_avg_n_Wrap(w_time(:,:,192:203),2),res) 

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
      printVarSummary(ploylevel1)
  

     
     dum1=new(3,graphic)
     dum2=new(3,graphic) 
     dum3=new(3,graphic) 

    ployres@gsLineDashPattern=0
    do i=0,2
    dum1(i)= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
    end do
    ployres@gsLineDashPattern=16
    do i=0,2
    dum2(i)= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    end do

    do i=0,2
    dum3(i)= gsn_add_polyline(wks,plot(i),  ploylon3,ploylevel3,  ployres)
    end do

;-------------------------------
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = 0
    res@cnLineThicknessF = 4.
    res@cnFillOn        = False
    res@gsnRightString=""
  plot1=new(3,graphic)

    plot1(0)= gsn_csm_pres_hgt(wks, msf_enino, res)
    plot1(1)= gsn_csm_pres_hgt(wks, msf_lanina, res)
    plot1(2)= gsn_csm_pres_hgt(wks, msf_mid, res)

    do i=0,2
    overlay(plot(i),plot1(i))
    end do
;;;

pres=True
pres@gsnPanelLabelBar=True
pres@gsnMaximize=False
gsn_panel(wks, plot, (/3,1/), pres)


;;;画1994年秋季的IOD事件和1997年的enino事件
   
    dum_iod1= gsn_add_polyline(wks1,plot_iod,  ploylon1,ploylevel1,  ployres)
    dum_iod2= gsn_add_polyline(wks1,plot_iod,  ploylon2,ploylevel2,  ployres)
    dum_iod3= gsn_add_polyline(wks1,plot_iod,  ploylon3,ploylevel3,  ployres)

    dum_enino1= gsn_add_polyline(wks1,plot_enino,  ploylon1,ploylevel1,  ployres)
    dum_enino2= gsn_add_polyline(wks1,plot_enino,  ploylon2,ploylevel2,  ployres)
    dum_enino3= gsn_add_polyline(wks1,plot_enino,  ploylon3,ploylevel3,  ployres)

    plot1_iod= gsn_csm_pres_hgt(wks1, msf_iod, res)
    plot1_enino= gsn_csm_pres_hgt(wks1, msf_str_enino, res)

    overlay(plot_iod, plot1_iod)
    overlay(plot_enino, plot1_enino)
    gsn_panel(wks1, (/plot_iod,plot_enino/), (/2,1/), pres)

;;;
draw(plot_1995)
frame(wks3)
end
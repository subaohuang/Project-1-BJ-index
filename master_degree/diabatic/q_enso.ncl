begin
diri="/mnt/f/share/ersst/"
fils0=systemfunc("ls "+diri+"nino3_4.ascii")
nino3_4=asciiread(fils0, 457, "float")
nino3_4_de=dtrend_n(nino3_4(0:407), True,0)
printVarSummary(nino3_4_de)
nino3_4_de!0="time"
nino3_4_de&time=ispan(1,408,1)

nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)


enino=ind(nino3_4_nor.gt.1.0)
lanila=ind(nino3_4_nor.lt.-1.0)
Neutral=ind(nino3_4_nor.le.1.0.and.nino3_4_nor.ge.-1.0)

diri1="/mnt/f/share/era/"
fils1=systemfunc("ls "+diri1+"Q_time_era.nc")
f1=addfile(fils1, "r")
Q_time=f1->Q_time(:,:,0:455)
Q_time_enino=Q_time(:,:,enino)
Q_time_lanina=Q_time(:,:,lanila)
Q_time_neutral=Q_time(:,:,Neutral)
Q_enino_com=dim_avg_n_Wrap(Q_time_enino,2)
Q_lanina_com=dim_avg_n_Wrap(Q_time_lanina,2)
Q_neutral_com=dim_avg_n_Wrap(Q_time_neutral, 2)
Q_clm=dim_avg_n_Wrap(Q_time, 2)

Q_time_de=dtrend_n(rmMonAnnCycLLT(Q_time),True,2)
Q_time_enino_de=dim_avg_n(Q_time_de(:,:,enino),2)
Q_time_lanina_de=dim_avg_n(Q_time_de(:,:,lanila),2)
Q_time_neutral_de=dim_avg_n(Q_time_de(:,:,Neutral),2)
Q_time_sum_de=Q_time_enino_de+Q_time_lanina_de

copy_VarCoords(Q_time_enino(:,:,1), Q_time_enino_de)
copy_VarCoords(Q_time_lanina(:,:,1), Q_time_lanina_de)
copy_VarCoords(Q_time_lanina(:,:,1), Q_time_neutral_de)
copy_VarCoords(Q_time_lanina(:,:,1), Q_time_sum_de)


plot=new(4,graphic)
 wks=gsn_open_wks("png", ""+diri1+"Q_era_enso")
  gsn_define_colormap(wks,"hotcolr_19lev")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="1979-2017"
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -3
    res@cnMaxLevelValF = 3
    res@cnLevelSpacingF = 0.5
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    ; res@cnMonoFillColor=False
    
    ; res@cnFillColors=(/8,11,14,18,48,50,67,74,84,93/)
    res@cnLineLabelsOn=False
    ; res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res@cnLabelScalingMode="ScaleFactor"
       ; res@cnLabelScaleValueF=10^9
    res@tmYRMode="Automatic"
    res@lbOrientation="Vertical"
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(100,290,40)
    res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    res@tmXBMinorOn=False

    plot(0)=gsn_csm_pres_hgt(wks,Q_enino_com,res)
    plot(1)=gsn_csm_pres_hgt(wks,Q_lanina_com,res)
    plot(2)=gsn_csm_pres_hgt(wks,Q_neutral_com,res)
    plot(3)=gsn_csm_pres_hgt(wks,Q_clm,res)
    pres=True
    pres@gsnMaximize=True

   ;---Attach the polyline
  ploylevel=1000
  ploylon1=fspan(100, 149, 25)
  ploylon2=fspan(280, 290, 11)
  ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
  ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
 
  dum2=new(4,graphic)  
  dum1=new(4,graphic)
 ;---Resources for the polyline  添加下垫面的陆地
  ployres                  = True          ; polyline mods desired
  ployres@gsLineThicknessF = 15.0            ; line thickness   
  ployres@gsLineColor      = "black"         ; color of lines
  ployres@gsLineDashPattern=15

   do i=0,3
 dum1(i)= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
  end do
ployres@gsLineDashPattern=0
do i=0,3
 dum2(i)= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel1,  ployres)
  end do
    gsn_panel(wks, plot, (/2,2/), pres)
   
delete(dum1)
delete(dum2)

 plot2=new(4,graphic)
 wks2=gsn_open_wks("png", ""+diri1+"Q_era_enso_ano")
    delete( res@cnLevelSelectionMode )
    delete(res@cnMinLevelValF )
    delete( res@cnMaxLevelValF ) 
    delete( res@cnLevelSpacingF )
   res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -2.0
    res@cnMaxLevelValF = 2.0
    res@cnLevelSpacingF = 0.5
    plot2(0)=gsn_csm_pres_hgt(wks2,Q_time_enino_de,res)
    plot2(1)=gsn_csm_pres_hgt(wks2,Q_time_lanina_de,res)

    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.3
    res@cnMaxLevelValF = 0.3
    res@cnLevelSpacingF = 0.05


    plot2(2)=gsn_csm_pres_hgt(wks2,Q_time_neutral_de,res)
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -2.0
    res@cnMaxLevelValF = 2.0
    res@cnLevelSpacingF = 0.5
 dum2=new(4,graphic)  
  dum1=new(4,graphic)

    plot2(3)=gsn_csm_pres_hgt(wks2,Q_time_sum_de,res)
    
    pres2=True
    pres2@gsnMaximize=True
 ;---Attach the polyline
   dum2=new(4,graphic)  
   dum1=new(4,graphic)
ployres@gsLineDashPattern=15

   do i=0,3
 dum1(i)= gsn_add_polyline(wks2,plot2(i),  ploylon1,ploylevel1,  ployres)
  end do
ployres@gsLineDashPattern=0
do i=0,3
 dum2(i)= gsn_add_polyline(wks2,plot2(i),  ploylon2,ploylevel1,  ployres)
  end do

      gsn_panel(wks2, plot2, (/2,2/), pres2)

;     res@cnLevelSelectionMode = "ExplicitLevels"
;     res@cnLinesOn=True
;     ; res@cnLevels = ispan(-50,400,50)
;     res@cnLineLabelsOn=True
;     res@cnLineLabelBackgroundColor="Transparent"
;     res@cnLineThicknessF = 3.
;     res@cnFillOn        = False
;     res@gsnContourNegLineDashPattern=1
;     plot1 = gsn_csm_pres_hgt(wks, Q_new, res)

; overlay(plot,plot1)
; draw(plot)
; frame(wks)
end
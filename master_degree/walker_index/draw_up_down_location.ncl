function gain_data(data1,data2)
 begin
 
 rc1=regCoef(data1&time,data1)

 rc2=regCoef(data2&time,data2)
  print(rc1)
  print(rc2)
 data_fuc=new((/2,4,dimsizes(data1&time)/),double)


 data_fuc(0,0,:) = data1(:)
 data_fuc(0,1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(0,2,:) =avg(data1(:))
 data_fuc(0,3,:)= rc1
 data_fuc(1,0,:) = data2(:)
 data_fuc(1,1,:) = rc2*data2&time+rc2@yintercept
 data_fuc(1,2,:) =avg(data2(:))
 data_fuc(1,3,:) =rc2
return(data_fuc)
end 

function draw_panel(wks1,wks2,time,data1,data2,title,rc1,rc2)
 begin
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnLeftString=title
    resxy@gsnRightString="upward location"
    resxy@tmXBMode="Explicit"
    ; resxy@trXMaxF=456
    resxy@tmXBValues=ispan(100,240,20)
    resxy@tmXBLabels=ispan(100,240,20)+"~S~o~N~E"
    resxy@tmXBMinorOn=False
    
    
    resxy@tmYLMode="Explicit"
    
    resxy@trYMaxF=dimsizes(data1&time)
    resxy@tmYLValues=ispan(13,dimsizes(data1&time),60)
    resxy@tmYLLabels=ispan(1980,2015,5)
    ; resxy@tmXMajorGrid=True
    ; resxy@tmXMajorGridLineDashPattern=1


    resxy@vpHeightF=0.7
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9

    resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
    resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    resxy@xyDashPatterns    =(/0,0,2/)
    

    resxy@tmXBLabelFontAspectF=2.5
    resxy@tmYLLabelFontAspectF=2.5
    resxy@tmXTOn=False
    resxy@tmYROn=False
    ; resxy@tiYAxisString = title

    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy@pmLegendSide           = "Top"               ; Change location of 
    resxy@pmLegendParallelPosF   = 0.85                  ; move units right
    resxy@pmLegendOrthogonalPosF = -0.18                ; move units down
    resxy@pmLegendWidthF         = 0.10                ; Change width and
    resxy@pmLegendHeightF        = 0.1               ; height of legend.
    resxy@lgPerimOn              = False               ; turn off box around
    resxy@lgLabelFontHeightF     = .008                 ; label font height
    resxy@xyExplicitLegendLabels = (/"Anual change", "linear trend"+"("+sprintf("%5.3f",rc1(1)*120)+"~S~o~N~/decade"+")","average"/)         ; create explicit labels
    




    ; resxy@tmLabelAutoStride = True                    ; nice stride on labels
    contour=new(2,graphic)
    resxy2=resxy
    delete(resxy2@tmXBValues)
    delete(resxy2@tmXBLabels)
    delete(resxy2@xyExplicitLegendLabels)

    delete(resxy2@pmLegendDisplayMode)
    delete(resxy2@pmLegendSide)
    delete(resxy2@pmLegendParallelPosF   )
    delete(resxy2@pmLegendOrthogonalPosF )
    delete(resxy2@pmLegendWidthF        )
    delete(resxy2@pmLegendHeightF       )
    delete(resxy2@lgPerimOn              )
    delete(resxy2@lgLabelFontHeightF    )
    delete(resxy2@gsnRightString)


    resxy2@gsnRightString="downward location"
    resxy2@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy2@pmLegendSide           = "Top"               ; Change location of 
    resxy2@pmLegendParallelPosF   = 0.15                  ; move units right
    resxy2@pmLegendOrthogonalPosF = -0.18                ; move units down
    resxy2@pmLegendWidthF         = 0.10                ; Change width and
    resxy2@pmLegendHeightF        = 0.1               ; height of legend.
    resxy2@lgPerimOn              = False               ; turn off box around
    resxy2@lgLabelFontHeightF     = .008                 ; label font height
    resxy2@xyExplicitLegendLabels = (/"Anual change", "linear trend"+"("+sprintf("%5.3f",rc2(1)*120)+"~S~o~N~/decade"+")","average"/)         ; create explicit labels
    

    resxy2@tmXBValues=ispan(180,280,20)
    resxy2@tmXBLabels=""+ispan(180,280,20)+"~S~o~N~E"
    resxy2@trXMaxF=280

    contour(0)=gsn_csm_xy(wks1,data1(:,:),time,resxy)
    contour(1)=gsn_csm_xy(wks2,data2(:,:),time,resxy2)
    
 
    ; contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    ; contour(1)=gsn_csm_xy(wks, time,data2,resxy)

return contour
end 

begin 
 diri="/mnt/f/share/era/"
 shuju="era"
 fils1=systemfunc("ls "+diri+"msf_era_level_location1979-2017.nc")
 f1=addfile(fils1, "r")
 locationxi_a=f1->locationxi_a
 locationdong_a=f1->locationdong_a
 locationxi=f1->locationxi
 printVarSummary(locationxi)
  print("locationxi=" +locationxi_a)

 printVarSummary(locationxi_a)

 wks1=gsn_open_wks("png", ""+diri+"locationxi_"+shuju+"_chazhi")
 wks2=gsn_open_wks("png", ""+diri+"locationdong_"+shuju+"_chazhi")

 
 data= gain_data(locationxi_a,locationdong_a)
    pres1=True
    pres2=True
    pres1@gsnMaximize=True
    contour=draw_panel(wks1,wks2,locationxi_a&time,data(0,0:2,:),data(1,0:2,:),shuju,data(0,3,:),data(1,3,:) )

    draw(contour(0))
    draw(contour(1))
    frame(wks1)
    frame(wks2)
end
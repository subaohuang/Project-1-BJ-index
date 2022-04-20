function gain_data(data1)
 begin
 
 rc1=regCoef(data1&time,data1)


 data_fuc=new((/3,dimsizes(data1&time)/),float)


 data_fuc(0,:) = data1(:)
 data_fuc(1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(2,:) =avg(data1)

return(data_fuc)
end 
begin
diri1="/mnt/f/share/era/"
fils1=systemfunc("ls "+diri1+"Q_time_era.nc")
f1=addfile(fils1, "r")
Q_time=f1->Q_time(:,:,0:455)
location=new(dimsizes(Q_time&time),float)

 do t=0,455
	lon=maxind(ndtooned(Q_time(:,:,t)))
   m= ind_resolve(lon,dimsizes(Q_time(:,:,t)) )
   ind_lo=m(0,1)
   location(t)=Q_time&lon(ind_lo)
   delete(lon)
   delete(m)
   delete(ind_lo)
end do
; end do
time=ispan(1,456,1)
 print(location)
location!0="time"
 location&time=ispan(1,456,1)
 printVarSummary(location)
 data=gain_data(location)
 wks=gsn_open_wks("png", ""+diri1+"Q_location")
 res=True
 plot=gsn_csm_xy(wks, data, time, res)
 diri2="/mnt/f/share/era/ceshi/"

 wks2=gsn_open_wks("png", ""+diri2+"Q_location")

delete(res)

    res = True
   
    res@gsnLeftString="era"
    res@gsnRightString="1979-2017"
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    ; res@cnMinLevelValF = -3
    ; res@cnMaxLevelValF = 3
    ; res@cnLevelSpacingF = 0.5
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
do i=0,455
plot=gsn_csm_pres_hgt(wks2,Q_time(:,:,i),res)
end do 

end

;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
undef("draw_msf_trend")
function draw_msf_trend(LeftString)
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
     levind=ind(f1->level.eq.500)

    msf_time_20cr=f1->msf_time(:,{100:310},timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    delete(f1) 
    ; delete(date)
    ; delete(timeind) 
;;;去趋势
   ; msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
   msf_time_20cr_de=msf_time_20cr
   copy_VarCoords(msf_time_20cr, msf_time_20cr_de)
;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
 msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")

 msf_son&time=ispan(1979,2012,1)
 msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr_de(:,:,{201212})/)
 ;;;

 msf_clm=dim_avg_n((/msf_djf(:,:,:)/)*1.0,0)
 copy_VarCoords(msf_son(1,:,:), msf_clm)
 msf_cold=msf_djf
 msf_cold=(/msf_djf(:,:,:)/) 

do i=0,33
    year=1979+i
    wks=gsn_open_wks("png", "/mnt/f/share/era/data/djf_every_pi/"+LeftString+"_djf_"+year)
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString=LeftString
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

    res@vpHeightF=0.3
    res@vpWidthF=0.7

   plot=gsn_csm_pres_hgt(wks,msf_cold(i,:,:),res)
   draw(plot)
   frame(wks)
end do

;;;画图
return [/plot/]
end

;***********************************************************
;主程序
;***********************************************************         
begin





;;;画图
lo1=draw_msf_trend("20cr")
lo2=draw_msf_trend("era")
lo3=draw_msf_trend("jra55")
lo4=draw_msf_trend("MERRA")
lo5=draw_msf_trend("ncep1")
lo6=draw_msf_trend("ncep2")
; 
; lo5=draw_msf_trend(fils5,wks,"mre",pa)



end







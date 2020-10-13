

;************************************************
; 画20cr,ncep,era的msf气候态
;************************************************
begin
    diri="/mnt/f/share/rea_msf/"
    fils1=systemfunc("ls "+diri+"msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind=ind(f1->level.eq.500)
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
;;;

;;;选择冬季
    msf_clm_20cr=dim_avg_n_Wrap( month_to_season(msf_time_20cr(time|:,level|:,lon|:),"JJA"), 0)
;;;location


;************************************************
; 画图气候态
;************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/MERRA再分析资料质量流函数")

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    ; res@gsnLeftString=LeftString
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

    res@vpHeightF=0.3
    res@vpWidthF=0.7

   plot=gsn_csm_pres_hgt(wks,msf_clm_20cr,res)

draw(plot)
frame(wks)
end




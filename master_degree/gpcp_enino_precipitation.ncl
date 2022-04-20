
;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
begin
    
    f1=addfile("/mnt/f/share/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    
    pre=f1->precip(timeind,:,:)
    
    pre=(/rmMonAnnCycTLL(pre)/)

    delete(pre&time)
    pre&time=date(timeind)
    delete(f1) 
    delete(date)
    delete(timeind)
    
    


    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    sst=lonFlip(f2->sst(timeind,{20:-20},:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    ssta_son=month_to_season(ssta_de, "SON")
    ssta_djf=month_to_season(ssta_de, "DJF")
    
     delete(ssta_de&time)
     ssta_de&time=date(timeind)

    ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:)=(/ssta_djf(1:,:,:)/)
    ssta_djf(dimsizes(ssta_djf&time)-1,:,:)=(/ssta_de({201212},:,:)/)

    ; ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
    ssta_s=(/ssta_djf(:,:,:)/)*1.0
    copy_VarCoords(ssta_djf, ssta_s)

    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
   
    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)

    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
;;;
;;;去趋势
   pre_de=dtrend_msg_n(ispan(1,dimsizes(pre&time),1),pre,False,True,0)
   copy_VarCoords(pre, pre_de)
;;;   
  
;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 
 pre_djf=month_to_season(pre_de(time|:,lat|:,lon|:), "DJF")

 
 pre_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    pre_djf(0:dimsizes(pre_djf&time)-2,:,:)=(/pre_djf(1:,:,:)/)
    pre_djf(dimsizes(pre_djf&time)-1,:,:)=(/pre_de({201212},:,:)/)
 ;;;
;;;
  pre_cold=(/pre_djf(:,:,:)*1.0/) 

;;;挑选年份 modoki year
 modokiind=ind(emi_std.ge.std_emi)
 pre_modoki=dim_avg_n((/pre_djf(modokiind,:,:)/)*1.0,0)
 ; pre_modoki=dim_avg_n((/pre_djf(modokiind+1,:,:)/),0)   ;;;;djf取的要加1
 copy_VarCoords(pre_djf(1,:,:), pre_modoki)
;;;

;;;挑选年份 canino  enino year
 caind=ind(nino3_std.ge.std_nino3)
 pre_canonical=dim_avg_n((/pre_djf(caind,:,:)/)*1.0,0)
 ; pre_canonical=dim_avg_n((/pre_djf(caind+1,:,:)/),0)
 copy_VarCoords(pre_djf(1,:,:), pre_canonical)
;;;

;;;挑选年份 canino  lanina year
 laind=ind(nino3_std.le.(/-1*std_nino3/))
 pre_lanina=dim_avg_n((/pre_djf(laind,:,:)/)*1.0,0)
 ; pre_canonical=dim_avg_n((/pre_djf(caind+1,:,:)/),0)
 copy_VarCoords(pre_djf(1,:,:), pre_lanina)
;;; normal time
 noind=ind(nino3_std.ge.(/-0.5*std_nino3/).and.nino3_std.le.(/0.5*std_nino3/))
 pre_noind=dim_avg_n((/pre_djf(noind,:,:)/)*1.0,0)
 ; pre_canonical=dim_avg_n((/pre_djf(caind+1,:,:)/),0)
 copy_VarCoords(pre_djf(1,:,:), pre_noind)

;;;挑选年份 lanina modoki year
 lamodokiind=ind(emi_std.le.(/-1*std_emi/))
 pre_lamodoki=dim_avg_n((/pre_djf(lamodokiind,:,:)/)*1.0,0)
 ; pre_canonical=dim_avg_n((/pre_djf(caind+1,:,:)/),0)
 copy_VarCoords(pre_djf(1,:,:), pre_lamodoki)
;;
 
;;;;-----------------------------------------
;;画图
;;;;-------------------------------------------
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/不同厄尔尼诺事件时期降水的区别")
  gsn_define_colormap(wks,"MPL_bwr")
     res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
    ; res@cnFillDrawOrder="PostDraw"
    ; res@cnFillPalette=cmap
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -6
    res@cnMaxLevelValF       =   6
    res@cnLevelSpacingF      =  0.25

    res@tmXTOn=False
    res@tmYROn=False
    
    ; res@mpFillDrawOrder="PreDraw"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -50.0
    res@mpMaxLatF             = 50.0
    res@mpMinLonF             = 90
    res@mpMaxLonF             = 330
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

   plot=new(4,graphic)
    plot(0) = gsn_csm_contour_map(wks,pre_canonical,res)
    plot(1) = gsn_csm_contour_map(wks,pre_lanina,res)

    
    delete(res@cnMinLevelValF)
    delete(res@cnMaxLevelValF)    
 
    ; res@cnMinLevelValF       =   -4
    ; res@cnMaxLevelValF       =   4
    ; res@cnLevelSpacingF      =  1
    plot(2) = gsn_csm_contour_map(wks,pre_modoki,res)
    plot(3) = gsn_csm_contour_map(wks,pre_lamodoki,res)
    pres=True
    pres@gsnPanelLabelBar=True

gsn_panel(wks, plot, (/2,2/), pres)

end

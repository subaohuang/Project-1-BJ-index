;;;************************************************************
;; investigate and draw the annual cycle of WC movemont
;；画不同季节的气候态的PWC的空间结构的变化情况
;; linshuheng 2019.10.16
;;;*********************************************************

;;;===================read data ==========================
;;;
begin

   fils1 = systemfunc("ls /home/ys17-19/lsh/data/rea_msf/msf_time_era_global.nc ")
    f1 = addfile(fils1, "r")

    time = f1->time
    date = cd_calendar(time, 1)

    tind = ind(date.le.201212.and.date.ge.197901)
    msf_time = f1->msf_time(:,:,tind)
    msf_time=(/msf_time/10^9/)

   ;;;去趋势
    msf_time = (/dtrend_msg_n(ispan(1,dimsizes(msf_time&time),1),msf_time,False,True,2)/)
   ;;;

    msf_time&time := date(tind)

   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_time(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_time(time|:,level|:,lon|:), "DJF")
    msf_mam=month_to_season(msf_time(time|:,level|:,lon|:), "MAM")
    msf_jja=month_to_season(msf_time(time|:,level|:,lon|:), "JJA")

    msf_son&time=ispan(1979,2012,1)
    msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time(:,:,{201212})/)

    ;;;取气候态
    msf_mam_clm = dim_avg_n_Wrap(msf_mam, 0)
    msf_jja_clm = dim_avg_n_Wrap(msf_jja, 0)
    msf_son_clm = dim_avg_n_Wrap(msf_son, 0)
    msf_djf_clm = dim_avg_n_Wrap(msf_djf, 0)

;;;===================draw==========================
    wks = gsn_open_wks("pdf", "沃克环流在不同季节的空间结构变化")

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False

    res2@cnLinesOn=True
    res2@cnLevelSelectionMode ="ExplicitLevels"
    res2@cnLevels=0.2
    res2@cnFillOn=False
    res2@cnLineColor="red"
    res2@cnMonoLineColor=True
    res2@cnLineThicknessF =4.0
    
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    

    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    ; res2@tmXBValues=(/0,90,180,270,357.5/)
    ; res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    ; res2@tmXBMinorOn=True
    ; res2@tmXBMinorValues=ispan(0,360,10)

    res2@tmXBValues=ispan(100,310,40)
    res2@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","60~S~o~N~W"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(100,310,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.35
    res2@vpWidthF=0.65
   

    plot=new(1,graphic)
    
    plot=gsn_csm_pres_hgt(wks,msf_jja_clm(:,{120:285}),res2)

    res2@cnLineColor="blue"
    plot1=gsn_csm_pres_hgt(wks,msf_djf_clm(:,{120:285}),res2)
    ; res2@cnLineColor="black"
    ; plot2=gsn_csm_pres_hgt(wks,msf_clm,res2)    
    res2@cnLineColor="pink"
    
    plot2=gsn_csm_pres_hgt(wks,msf_son_clm(:,{120:260}),res2)  
    res2@cnLineColor="green"
    plot3=gsn_csm_pres_hgt(wks,msf_mam_clm(:,{120:285}),res2)  
      
    overlay(plot,plot1)
    overlay(plot,plot2)
    overlay(plot,plot3)   
    
    draw(plot)
    frame(wks)



    end
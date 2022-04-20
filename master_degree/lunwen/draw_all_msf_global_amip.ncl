;draw_msf_hov
begin
;画时间序列图的时候，把BCC-CSM1-1给去掉
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)
; models=(/"BNU-ESM"/)
dim=dimsizes(models)


msf_total=new((/dim,12,144/), float)



yearstart=197901
yearend=200812  
year1=1979
year2=2008
; p  = 0.975
p=0.95
year=year2-year1+1
locationxi_rea=new((/dim,year/), float)
;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")
 diri="/mnt/f/share/CMIP5/msf_AMIP/"
 
 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")

;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind=get1Dindex(f1->level,z0)
   
    msf_time_r= f1->msf_time(levind,:,timeind)
    delete(f1)

    msf_time_r=(/msf_time_r/10^9/)
    msf_time=tofloat(msf_time_r)
    copy_VarMeta(msf_time_r, msf_time)
    printVarSummary(msf_time)
    delete(msf_time_r)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")

    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点


    ;;;
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ;;;


    ; msf_djf=msf_son
    
    
    msf_djf_clm = dim_avg_n_Wrap(msf_djf, 0)
    
    msf_djf_clm = dim_avg_n_Wrap(msf_mean(time|:,level|:,lon|:), 0)
    delete(msf_mean)

    delete(msf_djf)
    delete(msf_son)
    msf_total(i,:,:)=msf_djf_clm(:,:)
    delete(msf_djf_clm)
end do

msf=lonPivot(dim_avg_n_Wrap(msf_total, 0),30)
printVarSummary(msf)

;;;读入MRE的数据
 fils1="/mnt/f/share/rea_msf/msf_time_mre_global.nc"

    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    msf_time_r= f1->msf_time(levind,:,timeind)
    delete(f1)

    msf_time_r=(/msf_time_r/10^9/)
    msf_time=tofloat(msf_time_r)
    copy_VarMeta(msf_time_r, msf_time)
    printVarSummary(msf_time)
    delete(msf_time_r)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF

    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ;;;
    ; msf_djf=msf_son

    msf_djf_clm=dim_avg_n_Wrap(msf_djf, 0)
      msf_djf_clm = dim_avg_n_Wrap(msf_mean(time|:,level|:,lon|:), 0)

    delete(msf_mean)
    delete(msf_djf)
    msf_rea=lonPivot(msf_djf_clm,30)
    printVarSummary(msf_rea)
;;;;

  msf_diff=msf-msf_rea
  copy_VarCoords(msf, msf_diff)



;************************************************
; 画图气候态
;************************************************






wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/amip资料和再分析资料DJF平均2")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="(a)ENS"
   
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    ; res@cnFillPalette  = "BlueDarkRed18"
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 0.5
    res@cnFillOn = True 
    ; res@gsnYRefLine=500
    ; res@gsnYRefLineColor="black"

    ; res@gsnContourNegLineDashPattern=2
    ; res@gsnContourZeroLineThicknessF=4.0

    res@cnLineLabelsOn=False
    res@cnLineLabelInterval=2
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False


    res@lbLabelBarOn=False
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"

    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)
     
    res@tiYAxisString=""

    res@vpHeightF=0.3
    res@vpWidthF=0.7
    
    plot=new(3, graphic)
    plot(0)=gsn_csm_pres_hgt(wks, msf_rea(:,:), res)
    res@gsnLeftString="(b)MME"
    plot(1)=gsn_csm_pres_hgt(wks, msf(:,:), res)

    res@gsnLeftString=""
    plot(2)=gsn_csm_pres_hgt(wks, msf_diff(:,:), res)
    

    res@cnLevelSelectionMode ="ExplicitLevels"
    res@cnLevels=(/-0.5/)

    res@cnFillOn=False
    res@cnMonoLineColor=True
    res@cnLineColor="blue"
    res@cnLineDashPattern=0
    res@cnLineThicknessF =4.0
 

    plot1=gsn_csm_pres_hgt(wks,msf_rea(:,{31:180}),res)
    res@cnLineColor="black"
    res@cnLineDashPattern=0
    res@cnLevels=(/0/)
    plot2=gsn_csm_pres_hgt(wks,msf_rea(:,{120:}),res)

   overlay(plot(0), plot1)
   overlay(plot(0), plot2)

    res@cnLevelSelectionMode ="ExplicitLevels"
    res@cnLevels=(/-0.5/)
    res@cnFillOn=False
    res@cnMonoLineColor=True
    res@cnLineColor="blue"
    res@cnLineDashPattern=0
    res@cnLineThicknessF =4.0
    plot3=gsn_csm_pres_hgt(wks,msf(:,{31:180}),res)
    res@cnLineColor="black"
    res@cnLineDashPattern=0
    res@cnLevels=(/0/)
    plot4=gsn_csm_pres_hgt(wks,msf(:,{120:}),res)
   overlay(plot(1), plot3)
   overlay(plot(1), plot4)


 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 30, 25)+360.0
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon6=fspan(30, 39, 25)

    

      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel6=conform_dims(dimsizes(ploylon6), ploylevel, 0)
      

    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str31 = unique_string("dum31")
   
   
    str4 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str5 = unique_string("dum2")
    str6 = unique_string("dum3")
    str61 = unique_string("dum31")
    
    ployres@gsLineDashPattern=0
    
    i=0
 
    plot@$str1$= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    plot@$str2$= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot(i),  ploylon3,ploylevel3,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str31$= gsn_add_polyline(wks,plot(i),  ploylon6,ploylevel6,  ployres)

        i=1
 
    plot@$str4$= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    plot@$str5$= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str6$= gsn_add_polyline(wks,plot(i),  ploylon3,ploylevel3,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str61$= gsn_add_polyline(wks,plot(i),  ploylon6,ploylevel6,  ployres)


pres=True
pres@gsnPanelLabelBar=True
pres@lbLabelStride=2
gsn_panel(wks, plot, (/2,1/), pres)

end 
  
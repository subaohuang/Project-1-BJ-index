
;;****************************************************************
;;main program
;;*************************************************************

begin

    model="mre"
    diri="/mnt/f/share/rea_msf/"
    fils1=systemfunc("ls "+diri+"msf_time_"+model+"_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    yearstart=197901
    yearend=201212
    year1=1979
    year2=2012
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    msf_time_r=f1->msf_time(:,:,timeind)
    ud_time = f1->ud_time(:,:,timeind)

    ; msf_time_f=lonPivot((msf_time_r(time|:,level|:,lon|:)),180)
    msf_time_f=msf_time_r
    ; msf_time_f&lon=fspan(-180,177.5,143)
    printVarSummary(msf_time_f)
    msf_time=msf_time_f(level|:,lon|:,time|:)
    msf_time=(/msf_time/10^9/)
    msf_time&level@units="hpa"
    msf_time&time=ispan(0, dimsizes(timeind)-1, 1)

    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    delete(f1)
    delete(fils1)

 msf_clm=dim_avg_n_Wrap(msf_time, 2)

 ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 msf_son=month_to_season(msf_time(time|:,level|:,lon|:), "SON")
 msf_djf=month_to_season(msf_time(time|:,level|:,lon|:), "DJF")

 msf_son&time=ispan(year1,year2,1)
 msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time(:,:,{201212})/)
 ;;;

 msf_clm=dim_avg_n((/msf_djf(:,:,:)/)*1,0)
 copy_VarCoords(msf_son(1,:,:), msf_clm)

;;
; msf_djf_1984=msf_djf({2012},:,:)
 msf_djf_1984 =dim_avg_n_Wrap(msf_djf(:,:,:),0)

 ; fils1="/mnt/f/share/rea_msf/location_rea_r1i1p1_djf.nc"
 ; f1=addfile(fils1, "r")
 ; locationdong=f1->locationdong

 ; printVarSummary(locationdong)
 ; locationdong_merra=locationdong(1,:)
 ; print(locationdong_merra)

 ;;;求相关系数
  ud_jja = month_to_season(ud_time(time|:,level|:,lon|:), "MAM")
  printVarSummary(ud_jja)  

  ud_pa = dim_avg_n_Wrap(ud_jja(:,{850},{150:210}),1)
  
  ud_id = dim_avg_n_Wrap(ud_jja(:,{850},{60:90}), 1)

  corr = escorc(ud_pa, ud_id)

  print("corr="+corr)

  exit()




wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/"+model+"_"+yearstart+"-"+yearend+"")
;********************************************************
;;draw
;********************************************************

    cmap=read_colormap_file("MPL_PuOr")
    cmap=(/cmap(::-1,0:3)/)

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString=model
    res2@gsnRightString=""+yearstart+"-"+yearend
   
    res2@cnLinesOn=False

    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -6
    res2@cnMaxLevelValF = 6
    res2@cnLevelSpacingF = 1
    res2@cnFillPalette=cmap

    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    res2@lbLabelBarOn=True

    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    res2@tmXBValues=(/0,90,180,270,357.5/)
    res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(0,360,10)
    res2@tmYRMode="Automatic"

    res2@vpHeightF=0.3
    res2@vpWidthF=0.7
   
    plot=new(1,graphic)
    
    plot=gsn_csm_pres_hgt(wks,msf_djf_1984,res2)
  

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
     

        draw(plot)
        frame(wks)
  end 
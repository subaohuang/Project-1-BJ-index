

;************************************************
; 画20cr,ncep,era的msf气候态
;************************************************
function draw_all_msf(LeftString,wks,title)
begin
    diri="/mnt/f/share/rea_msf/"

    fils1=systemfunc("ls "+diri+"msf_time_"+LeftString+"_global.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind1=ind(f1->level.eq.400)
    levind2=ind(f1->level.eq.600)
    msf_time_20cr=f1->msf_time(:,:,timeind)
    ud_time_20cr=f1->ud_time(:,:,timeind)
    w_time_20cr=f1->w_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
    
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)

    printVarSummary(ud_time_20cr)
    printVarSummary(msf_time_20cr)
   
    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)
;;;
      ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
      msf_son=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "SON")
      msf_jja=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "JJA")
      msf_djf=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "DJF")
      msf_year=month_to_annual(msf_time_20cr(time|:,level|:,lon|:), 1)
      
      msf_son&time=ispan(1979,2012,1)
      msf_year&year=ispan(1979,2012,1)
      msf_jja&time=ispan(1979,2012,1)
      msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
      ;;注意点
      msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
      msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr(:,:,{201212})/)
      ;;;
      printVarSummary(msf_djf)
      msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
      msf_djf_1=dim_avg_n_Wrap(msf_djf({1979:1994},:,:), 0)
      msf_djf_2=dim_avg_n_Wrap(msf_djf({1995:2012},:,:),0)

      msf_jja_1=dim_avg_n_Wrap(msf_jja({1979:1994},:,:), 0)
      msf_jja_2=dim_avg_n_Wrap(msf_jja({1995:2012},:,:),0)

      msf_year_1=dim_avg_n_Wrap(msf_year({1979:1994},:,:), 0)
      msf_year_2=dim_avg_n_Wrap(msf_year({1995:2012},:,:),0)


    ud_clm_20cr=dim_avg_n_Wrap(month_to_season(ud_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
    w_clm_20cr=dim_avg_n_Wrap(month_to_season(w_time_20cr(time|:,level|:,lon|:),"DJF"), 0)

    
;************************************************
; 画图气候态
;************************************************

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    
   
    res@gsnMaximize=False
    res@gsnRightString=title
    res@cnLinesOn=True
    res@cnInfoLabelOn=False
    res@lbLabelBarOn=False
    res@cnLineLabelsOn=False
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

    res@cnLevelSelectionMode ="ExplicitLevels"
    res@cnLevels=(/1.0,-1.0/)

    res@cnFillOn=False
    res@cnLineColor="blue"
    res@cnMonoLineColor=True
    res@cnLineThicknessF =4.0


    res@vpHeightF=0.4
    res@vpWidthF=0.6
    plot=gsn_csm_pres_hgt(wks, msf_djf_1, res)
    res@cnLineColor="red"
    plot1=gsn_csm_pres_hgt(wks, msf_djf_2, res)
   
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
      




    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
  
    ployres@gsLineDashPattern=0
     
    plot@$str1$= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)

    ployres@gsLineDashPattern=16

    plot@$str2$= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)

    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)
    
    overlay(plot, plot1)
return(plot)
end
;;;;
;;主程序
;;;

    begin
    rea= (/"mre","era","jra55","MERRA","ncep1","ncep2"/)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/6个再分析资料DJF平均_全球2个时间段的比较_1.0")

    plot=new(6,graphic)
    title=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)
    do i=0,5
    plot(i)=draw_all_msf(rea(i),wks,title(i))
    end do
    pres=True
    gsn_panel(wks, plot, (/3,2/), pres)
 end
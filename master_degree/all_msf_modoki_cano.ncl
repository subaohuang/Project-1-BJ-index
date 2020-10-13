function decade_trend(data)
    local data_a, data_a_new, rc
    begin

    ; data_new=doubletofloat(dim_standardize_n_Wrap(data(::-1,:,:),1,(/0,1/)))       ;标准化
    ; print(data_new)
    ; copy_VarMeta(data(::-1,:,:), data_new)
    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)
    ;-----------------------------------------------
    rc = regCoef_n(data&time,data,0,2)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data&level
    rc&lon=data&lon
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
    ; print(rc)
return(rc)
end
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end

;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
undef("draw_msf_trend")
function draw_msf_trend(fils1,wks,LeftString)
begin
    diri="/mnt/f/share/"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    
    msf_time_20cr=f1->msf_time(:,{100:310},timeind)
    ; ud_time_20cr=f1->ud_time(:,:,timeind)
    ; w_time_20cr=f1->w_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
   
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    delete(f1) 

;;;

  
;;;去趋势
   msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
   ; msf_time_20cr_de=msf_time_20cr
;;;
   copy_VarCoords(msf_time_20cr, msf_time_20cr_de)
  printVarSummary(msf_time_20cr_de)
  
  

;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
 msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")
 

 msf_son&time=ispan(1958,2012,1)
 msf_djf&time=ispan(1958,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

 printVarSummary(msf_son)
 printVarSummary(msf_djf)

 msf_clm=dim_avg_n((msf_son(:,:,:)+msf_djf(:,:,:))*0.5, 0)
 copy_VarCoords(msf_son(1,:,:), msf_clm)
;;

;;;挑选年份 modoki year
 modoki_yr=(/1963,1968,1977,1979,1987,1990,1991,1992,1994,2002,2004,2009/)
 modokiind=get1Dindex(msf_son&time,modoki_yr)
 msf_modoki=dim_avg_n((/(msf_son(modokiind,:,:)+msf_djf(modokiind+1,:,:))*0.5/),0)

 ; msf_modoki=dim_avg_n((/msf_djf(modokiind+1,:,:)/),0)   ;;;;djf取的要加1
 copy_VarCoords(msf_son(1,:,:), msf_modoki)
;;;

;;;挑选年份 canino year
 canonical_yr=(/1965,1972,1976,1982,1997/)
 caind=get1Dindex(msf_son&time,canonical_yr)
 msf_canonical=dim_avg_n((/(msf_son(caind,:,:)+msf_djf(caind+1,:,:))*0.5/),0)
 
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_son(1,:,:), msf_canonical)
;;;



;;;画图

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString=LeftString
    res2@cnLinesOn=True
    res2@cnLevelSelectionMode ="ExplicitLevels"
    res2@cnLevels=0.0
    res2@cnFillOn=False
    res2@cnLineColor="red"
    res2@cnMonoLineColor=True
    res2@cnLineThicknesses=3.0
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
     
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7
   
    plot=new(1,graphic)
    
    plot(0)=gsn_csm_pres_hgt(wks,msf_canonical,res2)
    
    res2@cnLineColor="green"

    plot1=gsn_csm_pres_hgt(wks,msf_modoki,res2)
   
   res2@cnLineColor="black"

    plot2=gsn_csm_pres_hgt(wks,msf_clm,res2)    
    overlay(plot,plot1)
   
   overlay(plot,plot2)

return(plot)
end

;***********************************************************
;主程序
;***********************************************************         
begin

fils1="ncep/20CR/msf_time_20cr_global_6--6.nc"
fils2="ncep/ncep1/msf_time_ncep1_global.nc"
fils3="era/era40/msf_time_era40_in_global.nc"
fils4="JRA-55/msf_time_jra55_global.nc"

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个再分析资料_modoki和canomsf的位置秋冬")

plot=new(4,graphic)
plot(0)=draw_msf_trend(fils1,wks,"20CR")
plot(1)=draw_msf_trend(fils2,wks,"NCPE1")
plot(2)=draw_msf_trend(fils3,wks,"ERA")
plot(3)=draw_msf_trend(fils4,wks,"JRA-55")

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
 
    dum1=new(4,graphic)
    dum2=dum1
    dum3=dum1

    do i=0,3
    ployres@gsLineDashPattern=0
    dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
    end do 

    pres=True
    pres@gsnPanelLabelBar=False

    gsn_panel(wks, plot, (/2,2/), pres)
end




;计算3D Q 非绝热加热项Q的情况，era数据
function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin


     data_a=data
    ; data_a2=dim_standardize_n_Wrap(data_a,1,0)        ;标准化
    ; data_a1=doubletofloat(data_a2)

    ; data_a1&time=ispan(1,456,1)
    data_a&time=ispan(1,dimsizes(data&time),1)
    data_a@units=units
    data_a@long_name=long_name
    data_a_new = data_a(level|:,longitude|:,time|:)

    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)

    ;-----------------------------------------------
    rc = regCoef(data_a_new&time,data_a_new)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data_a_new&level
    rc&lon=data_a_new&longitude
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
    rc=(/rc*10.0/)
return(rc)
end
begin
 diri="/mnt/f/share/era/data/data_series/"
 fils1=systemfunc("ls "+diri+"era_w.nc")
 fils2=systemfunc("ls "+diri+"era.T.nc")
 f1=addfile(fils1, "r")
 f2=addfile(fils2, "r")
 time=f1->time
 lon=f1->longitude
 lat=f1->latitude
 level=f1->level
 date=cd_calendar(time, 1)

 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201201)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)

 levind1=ind(level.eq.100)
 levind2=ind(level.eq.1000)
 
 omega=f1->w(tind1:tind2,levind1:levind2,latind1:latind2,:)
 w=dble2flt(omega*omega@scale_factor+omega@add_offset)
 copy_VarMeta(omega, w)
 delete(omega)

 tem=f2->t(tind1:tind2,levind1:levind2,latind1:latind2,:)
 t=dble2flt((/tem*tem@scale_factor+tem@add_offset/))
 copy_VarMeta(tem, t)

 delete(tem)
 p=level(levind1:levind2)*100.0

 pt=pot_temp(p, t, 1, False)
 printMinMax(pt, 1)
                                                          
 dtdp=center_finite_diff_n(pt, p,False , 0, 1)                                                         ;计算稳定度
 copy_VarCoords(pt, dtdp)

 printMinMax(w, 1)
 printVarSummary(dtdp)
 Q=new((/dimsizes(pt&time),dimsizes(pt&level),dimsizes(pt&latitude),dimsizes(pt&longitude)/),float)
 do i=0,dimsizes(p)-1
 Q(:,i,:,:)=((p(i)/100000.0)^0.286)*w(:,i,:,:)*dtdp(:,i,:,:)*3600.0*24.0
 end do
 copy_VarCoords(pt, Q)
 Q_new=dim_avg_n_Wrap(Q(time|:,latitude|:,level|:,longitude|:),(/0,1/))
 Q_time=dim_avg_n_Wrap(Q(time|:,latitude|:,level|:,longitude|:), 1)
 rc_Q=decade_trend(Q_time,"k/day","dabatic heating")
 ; print(Q_new)
print(Q_new)
 Q_time_w=Q_time(level|:,longitude|:,time|:)
 wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/Q_era")
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="1979-2017"
    ; res@lbLabelStrings=ispan(-6,6,1)
  
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.1
    res@cnMaxLevelValF = 0.1
    res@cnLevelSpacingF = 0.01
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
    ; res@tmXBMode="Explicit"
    ; res@tmXBValues=ispan(100,290,40)
    ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    ; res@tmXBMinorOn=False

    plot=gsn_csm_pres_hgt(wks,rc_Q,res)


    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = ispan(-50,400,50)
    res@cnLineLabelsOn=True
    res@cnLineLabelBackgroundColor="Transparent"
    res@cnLineThicknessF = 3.
    res@cnFillOn        = False
    res@gsnContourNegLineDashPattern=1
    plot1 = gsn_csm_pres_hgt(wks, Q_new, res)

overlay(plot,plot1)
draw(plot)
frame(wks)
; nlev=dimsizes(Q_time_w&level)
; nlon=dimsizes(Q_time_w&longitude)
; ntime=dimsizes(Q_time_w&time)
; ;****************************************************************************
; ;随时间变化的量
;     diro=diri
;     filo1="Q_time_era.nc"
;     system("/bin/rm -f "+diro + filo1+"")
;     fout=addfile(diro + filo1,"c")
;     setfileoption(fout,"DefineMode",True)
;    ;======================================文件名
;     fAtt               = True            ; assign file attributes
;     fAtt@title         = "Q的值，取5到-5平均了"    
;     fAtt@creation_date = "2017.8.25"
;     fileattdef( fout, fAtt)            ; copy file attributes    
;     ;----------------------------坐标变量的定义
;     dimnames=(/"level","lon","time"/)
;     dimsize=(/nlev,nlon,ntime/)
;     dimunlim=(/False,False,False/)
;     filedimdef(fout, dimnames, dimsize, dimunlim)
;  ;-----------------------------定义变量名，类型和维数
;     filevardef(fout, "level",typeof(level), (/"level"/))
;     filevardef(fout, "lon", typeof(lon),(/"lon"/))
;     filevardef(fout, "time", typeof(time),(/"time"/))
;     filevardef(fout, "Q_time", typeof(Q_time_w),(/"level","lon","time"/))    
 
    
  
;     filevarattdef(fout, "level", Q_time_w&level)
;     filevarattdef(fout, "lon", Q_time_w&longitude)
   
;     setfileoption(fout,"DefineMode",False)

;     fout->level=(/Q_time_w&level/)
;     fout->lon=(/Q_time_w&longitude/)
;     fout->time=(/Q_time_w&time/)
;     fout->Q_time=(/Q_time_w/)

 end
;******************************************************************
 ;IOD & IOD modoki sst
;***************************************************************
;;选择月份的function
function iod_select(iod_year,month)                                                          ;
begin
 iod_month=new((/dimsizes(iod_year)/),string)
   if month.le.9 then
  do i=0,dimsizes(iod_year)-1
  iod_month(i)=""+iod_year(i)+""+0+""+month+""
  end do
  end if
    if month.gt.9.and.month.le.12 then
  do i=0,dimsizes(iod_year)-1
  iod_month(i)=""+iod_year(i)+""+""+month+""
  end do
  end if

   if month.gt.12.and.month.le.21 then
  do i=0,dimsizes(iod_year)-1
  iod_month(i)=""+iod_year(i)+""+0+""+(/month-12/)+""
 
  end do
  end if

  if month.gt.21.and.month.le.24 then
  do i=0,dimsizes(iod_year)-1
  iod_month(i)=""+iod_year(i)+""+""+(/month-12/)+""
  end do
  end if
  iod_month_new= stringtointeger(iod_month)

return iod_month_new
end 
;;----

begin
;;;读入sst
	diri="/mnt/f/share/ersst/"
	diri_pi="/mnt/f/share/era/"
	f0=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
	date=f0->time
	lon=f0->longitude
	lat=f0->latitude

	time=cd_calendar(date, 1)

	timeind=ind(time.le.201312.and.time.ge.195001)
  time_new=time(timeind)    ;;;所选取时间的日期
	latind=ind(lat.le.30.5.and.lat.ge.-30.5)
	sst_d=lonFlip(f0->sst(timeind,latind,:))
    ; sst_r=sst_d(:,:,{30:300})
    sst_r=sst_d(:,:,:)
    delete(sst_d)
;;;
    printVarSummary(sst_r)
;;;海温取年循环，去趋势，做3个月平滑处理
   ssta=dtrend_msg_n(ispan(0,dimsizes(sst_r&time)-1,1),rmMonAnnCycLLT(sst_r(latitude|:,longitude|:,time|:)),\
                     False,False,2)
   copy_VarCoords(sst_r(latitude|:,longitude|:,time|:), ssta)
 
   ssta=(/runave_n_Wrap(ssta, 3, -1, 2 )/)
   printVarSummary(ssta)
;;;----读入u，v, noaa 20cr数据
    diri_20c="/mnt/f/share/ncep/20CR/"   
    fils1=systemfunc("ls "+diri_20c+"uwnd.mon.mean.nc")
    fils2=systemfunc("ls "+diri_20c+"vwnd.mon.mean.nc")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date_20c=f1->time
    time_20c=cd_calendar(date_20c, 1)

    timeind_20c=ind(time_20c.le.201212.and.time_20c.ge.195001)
    time_20c_new=time_20c(timeind_20c)
    ; print(time_20c_new)
    u=f1->uwnd(timeind_20c,{850},{30:-30},:)
    v=f2->vwnd(timeind_20c,{850},{30:-30},:)
;;;

;;;速度去循环和趋势，3个月平滑处理
    ua=dtrend_msg_n(ispan(0,dimsizes(u&time)-1,1),rmMonAnnCycLLT(u(lat|:,lon|:,time|:)),\
                         False,False,2)
    copy_VarCoords(u(lat|:,lon|:,time|:), ua)
    ua=(/runave_n_Wrap(ua, 3, -1, 2 )/)
     va=dtrend_msg_n(ispan(0,dimsizes(v&time)-1,1),rmMonAnnCycLLT(v(lat|:,lon|:,time|:)),\
                         False,False,2)
    copy_VarCoords(v(lat|:,lon|:,time|:), va)
    va=(/runave_n_Wrap(va, 3, -1, 2 )/)
;;; 


;;; 挑选IOD事件 传统IOD：1961，1972,1982,1987,1997,2011,2012 正事件 ; IOD modoki :1963,1967,1976,1991,1994,2003,2008
    
    ; iod_year=(/1961,1972,1982,1987,1997,2011,2012/)
    iod_year=(/1961,1972,1982,1987,1997,2011/)
    iod_year1=iod_year+1

    iod_modoki_year=(/1963,1967,1976,1991,1994,2003,2008/)
    iod_modoki_year1=iod_modoki_year+1
    month=(/5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12/)

    iod_yr_mon=new((/dimsizes(iod_year)*dimsizes(month)/),string)
    iod_modoki_yr_mon=new((/dimsizes(iod_year)*dimsizes(month)/),string)

   ; do i=0,dimsizes(iod_year)-1
   ;  do j=0,dimsizes(month)-1
   ;    if (j+5).le.9 then 
   ;     iod_yr_mon(i*dimsizes(month)+j)=""+iod_year(i)+""+0+""+month(j)+""
   ;     iod_modoki_yr_mon(i*dimsizes(month)+j)=""+iod_modoki_year(i)+""+0+""+month(j)+""
   ;   else 
   ;    iod_yr_mon(i*dimsizes(month)+j)=""+iod_year(i)+""+""+month(j)+""
   ;     iod_modoki_yr_mon(i*dimsizes(month)+j)=""+iod_modoki_year(i)+""+month(j)+""
   ;    end if
   ;   end do
   ; end do
;;;

;;;选择相应的iod事件对应的月份 这里取出5-12月

    iod_time=new((/dimsizes(iod_year),dimsizes(month)/),integer)
    iod_modoki_time=new((/dimsizes(iod_modoki_year),dimsizes(month)/),integer)

    do i=0,dimsizes(month)-1
      if i.le.7 then 
        iod_time(:,i)=iod_select(iod_year,i+5)
        iod_modoki_time(:,i)=iod_select(iod_modoki_year,i+5)  ;;;而维数组，后面维是表示从5月到来年4月
      else
        iod_time(:,i)=iod_select(iod_year1,i+5)
        iod_modoki_time(:,i)=iod_select(iod_modoki_year1,i+5)  ;;;而维数组，后面维是表示从5月到来年4月
      end if
    end do
;;;


;;提取温度场成分，并平均
  iod_ssta=new((/dimsizes(ssta&latitude),dimsizes(ssta&longitude),dimsizes(month)/),typeof(ssta))
  iod_modoki_ssta=new((/dimsizes(ssta&latitude),dimsizes(ssta&longitude),dimsizes(month)/),typeof(ssta))
  do i=0,dimsizes(month)-1
   iod_ssta(:,:,i)=dim_avg_n_Wrap(ssta(:,:,get1Dindex(time_new,iod_time(:,i))),2)
   iod_modoki_ssta(:,:,i)=dim_avg_n_Wrap(ssta(:,:,get1Dindex(time_new,iod_modoki_time(:,i))),2)
  end do
  
  printVarSummary(iod_ssta)
  printVarSummary(iod_modoki_ssta)
;;;

;;;提取速度场成分
    iod_ua=new((/dimsizes(ua&lat),dimsizes(ua&lon),dimsizes(month)/),typeof(ua))
    iod_modoki_ua=new((/dimsizes(ua&lat),dimsizes(ua&lon),dimsizes(month)/),typeof(ua))
    do i=0,dimsizes(month)-1
    iod_ua(:,:,i)=dim_avg_n_Wrap(ua(:,:,get1Dindex(time_20c_new,iod_time(:,i))),2)
    iod_modoki_ua(:,:,i)=dim_avg_n_Wrap(ua(:,:,get1Dindex(time_20c_new,iod_modoki_time(:,i))),2)
    end do

    iod_va=new((/dimsizes(va&lat),dimsizes(va&lon),dimsizes(month)/),typeof(va))
    iod_modoki_va=new((/dimsizes(va&lat),dimsizes(va&lon),dimsizes(month)/),typeof(va))
    do i=0,dimsizes(month)-1
    iod_va(:,:,i)=dim_avg_n_Wrap(va(:,:,get1Dindex(time_20c_new,iod_time(:,i))),2)
    iod_modoki_va(:,:,i)=dim_avg_n_Wrap(va(:,:,get1Dindex(time_20c_new,iod_modoki_time(:,i))),2)
    end do
  



wks=gsn_open_wks("pdf", diri_pi+"data/"+"sst_iod_05_04")
wks2=gsn_open_wks("pdf", diri_pi+"data/"+"sst_iodmodoki05_04")
wks3=gsn_open_wks("pdf", diri_pi+"data/"+"sst_iodmodoki下一年05_12")

;;;画海温场
 gsn_define_colormap(wks, "MPL_seismic")
 gsn_define_colormap(wks2, "MPL_seismic")
 gsn_define_colormap(wks3, "BlRe")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillOn=True
    res@cnLinesOn=False
    res@gsnAddCyclic=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -1.5
    res@cnMaxLevelValF = 1.5
    res@cnLevelSpacingF = 0.25
    res@gsnScalarContour=True

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(iod_ssta&latitude)
    res@mpMaxLatF             = max(iod_ssta&latitude)
    res@mpMinLonF             = min(iod_ssta&longitude)
    res@mpMaxLonF             = max(iod_ssta&longitude)
    res@mpCenterLonF          = avg(iod_ssta&longitude)
    res@lbLabelBarOn=False

    res@vcRefMagnitudeF          = 3.0           ; define vector ref mag
    res@vcRefLengthF             = 0.02        ; define length of vec ref
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcMinDistanceF           = 0.01         ; thin vectors
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角
    res@vcMinMagnitudeF =0.5

   
    plot1=new(12,graphic)
    plot2=new(12,graphic)
    plot3=new(8,graphic)





  do i=0,11
   plot1(i)=gsn_csm_vector_scalar_map(wks, iod_ua(:,:,i),iod_va(:,:,i), iod_ssta(:,:,i), res)
   plot2(i)=gsn_csm_vector_scalar_map(wks2, iod_modoki_ua(:,:,i) , iod_modoki_va(:,:,i) ,iod_modoki_ssta(:,:,i), res)
end do

  pres=True
  pres@gsnMaximize=True
  pres@gsnPanelLabelBar=True
  pres@gsnPanelFigureStrings=(/"5","11","6","12","7","1","8","2","9","3","10","4"/)
  pres@gsnPanelFigureStringsJust="RightBottom"
  gsn_panel(wks,(/plot1(0),plot1(6),plot1(1),plot1(7),plot1(2),plot1(8),plot1(3),plot1(9),plot1(4),plot1(10),plot1(5),\
                  plot1(11)/), (/6,2/), pres)
   
  gsn_panel(wks2, (/plot2(0),plot2(6),plot2(1),plot2(7),plot2(2),plot2(8),plot2(3),plot2(9),plot2(4),plot2(10),\
                   plot2(5),plot2(11)/), (/6,2/), pres)
  
  
  do i=0,7
   plot3(i)=gsn_csm_vector_scalar_map(wks3, iod_modoki_ua(:,:,i+12) , iod_modoki_va(:,:,i+12) ,iod_modoki_ssta(:,:,i+12), res)
end do
  delete(pres@gsnPanelFigureStrings)
   pres@gsnPanelFigureStrings=(/"5","9","6","10","7","11","8","12"/)
  gsn_panel(wks3,(/plot3(0),plot3(4),plot3(1),plot3(5),plot3(2),plot3(6),plot3(3),plot3(7)/),(/4,2/),pres)



 end

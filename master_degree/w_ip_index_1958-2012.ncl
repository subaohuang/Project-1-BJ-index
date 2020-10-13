;;二维的情况
function malon2(X,case)
 begin
 lon_min=new(dimsizes(X&time),float)
 if case.eq."min" then 
     do t=0,dimsizes(X&time)-1
        dims = dimsizes(X(:,t))
        x1d = ndtooned(X(:,t))      ; convert 2D array to 1D for use in minind
        inds = ind_resolve(minind (x1d), dims)    ; convert 1D array back to 2D 
      
        ilon =inds(0,0)        ; select the longitude index where the X array is at its' minimum  
        lon_min(t)= X&lon(ilon)   ; insert the longitude index into the lon coordinate variable
     end do
 end if 

 if case.eq."max" then 
      do t=0,dimsizes(X&time)-1
      dims = dimsizes(X(:,t))
      x1d = ndtooned(X(:,t))      ; convert 2D array to 1D for use in minind
      inds = ind_resolve(maxind (x1d), dims)    ; convert 1D array back to 2D 
      ilon = inds(0,0)        ; select the longitude index where the X array is at its' minimum  
      lon_min(t) = X&lon(ilon)   ; insert the longitude index into the lon coordinate variable
    end do
 end if 
 return(lon_min)
 end
;;;计算每个月的气候态
function clm_month(x[*]:numeric)  
   local aveMonth
   begin
   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check; Compute all 12 monthly averages first. 
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do
   copy_VarAtts (x,aveMonth)
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   return (aveMonth)
end

function clm_ano(x[*]:numeric)  ;去除一维时间序列的年循环
begin

   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check; Compute all 12 monthly averages first. 
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do
   copy_VarAtts (x,aveMonth)
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   xAnom=x
   do yr=0,ntim-1,nmos
      xAnom(yr:yr+nmos-1) = (/ x(yr:yr+nmos-1)- aveMonth /)
   end do
   return (xAnom)
end

undef("rmtrend") ;;;滑动趋势
function rmtrend(x,nave,year)
begin
length=dimsizes(x&time)
rc=new((/length-nave+1/),typeof(x))
rc_test=new((/length-nave+1/),typeof(x))

if nave.le.21 then 
  siglv=0.1
  else 
  siglv=0.2
end if 

do t=nave,length  
   ssta_rc=regCoef_n(ispan(1,nave,1),x(t-nave:t-1),0,0)
   
   tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
 
   if prob.le.siglv
    rc_test(t-nave)=ssta_rc
   end if
   rc(t-nave)=ssta_rc 
end do 
rc!0="time"
rc&time=ispan(year+(nave-1)/2, (year+(nave-1)/2+length-nave), 1)
rc@rc_test=rc_test
return (rc)
end

;;;
begin

diri="/mnt/f/share/era/data/"
fils1=systemfunc("ls "+diri+"msf_time_era_global.nc")
f1=addfile(fils1, "r")
fils2=systemfunc("ls /mnt/f/share/ncep/20CR/omega.mon.mean.nc") 
; fils2=systemfunc("ls /mnt/f/share/era/data/era_20cm_pre_globle.nc") 
f2=addfile(fils2, "r")
date=cd_calendar(f2->time, 1)


timeind=ind(date.le.201212.and.date.ge.195801)
w=f2->omega(timeind,:,{10:-10},:)
; w_r=f2->w(timeind,:,{5:-5},:)
; w=w_r*w_r@scale_factor+w_r@add_offset
; copy_VarCoords(w_r, w)
; w!3="lon"
; w&lon=w_r&longitude
; delete(w_r)
w_r=dim_avg_n_Wrap(w,2)
printVarSummary(w_r)


;;;去除hadley的影响
  ; w_had=dim_avg_n_Wrap(w,3)

  ; w_ano=w

  ; do i=0,dimsizes(w&lon)-1
  ; w_ano(:,:,:,i)=(/w(:,:,:,i)-w_had(:,:,:)/)
  ; end do

  ; printVarSummary(w_ano)

  ; w_r=dim_avg_n_Wrap(w_ano, 2)
;;;


w_time=w_r(level|:,lon|:,time|:)
w_time=(/-100*w_time/)
delete(w_r)
printVarSummary(w_time)


;;根据w定义的指数
 ; w_indian=w_time({500},{40:120},:)
 ; w_pacific=w_time({500},{125:230},:)

w_indian=dim_avg_n_Wrap(w_time({850:200},{40:120},:),0)
w_pacific=dim_avg_n_Wrap(w_time({850:200},{125:230},:),0)



;;;




;;;;
 in_lon_w=malon2(w_indian,"max")
 pa_lon_w=malon2(w_pacific,"max")
 
 ipindex_w=pa_lon_w - in_lon_w

 
 opt=diri+"ip_index_1958_2012.ascii"
 opt1=diri+"pa_index_1958_2012.ascii"
 opt2=diri+"id_index_1958_2012.ascii"
 asciiwrite(opt, sprintf ("%6.2f,", ipindex_w))
 asciiwrite(opt1, sprintf ("%6.2f,", pa_lon_w))
 asciiwrite(opt2, sprintf ("%6.2f,", in_lon_w))

  ; ipindex_w=(/dtrend(ipindex_w, True)/)

;;; 月份变化
  ipindex_season=clm_month(ipindex_w)
  ipindex_season!0="time"
  ipindex_season&time=ispan(1, 12, 1)
;;;



;;  
    ipindex_anu=month_to_annual(ipindex_w, 1)
    ipindex_anu_std=dim_standardize(ipindex_anu, 1)
    ipindex_anu_std!0="time" 
    ipindex_anu_std&time=ispan(1, 55, 1)
;; 




;;;读入nino3.4 index
    ; nino3_4_r=asciiread(diri+"Hadley_nino3_4_1950_2010_nodtrend.ascii", 732, "float")
    
    ; ; nino3_4=(/dtrend(nino3_4, True)/)
    
    ; nino3_4=nino3_4_r()

    ; nino3_4_anu=month_to_annual(nino3_4, 1)
    ; nino3_4_anu_std=dim_standardize(nino3_4_anu, 1)
    ; nino3_4_anu_std!0="time" 
    ; nino3_4_anu_std&time=ispan(1, 55, 1)
    
;;;读入pdo index
delete(fils1)
delete(f1)
delete(date)
delete(timeind)
  fils1=systemfunc("ls /mnt/f/share/era/data/pdo.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->TSAXIS, 1)
  timeind=ind(date.ge.195801.and.date.le.201212)
  pdo=f1->PDO(timeind)
  pdo!0="time"
  printVarSummary(pdo)
  pdo_anu=month_to_annual(pdo, 1)
  pdo_anu&year=ispan(1, 55, 1)


    
; ;;;相关性
    ; r1=escorc_n(nino3_4_anu_std, ipindex_anu_std, 0, 0)
     
    ; r2=escorc_n(runave_n_Wrap(nino3_4_anu_std,9,0,0), runave_n_Wrap(ipindex_anu_std, 9, 0, 0), 0, 0)

    ; mxlag    = 120
    ; x_Lead_y = esccr(ipindex,nino3_4,mxlag)
    ; y_Lead_x = esccr(nino3_4,ipindex,mxlag)    ; switch the order of the series

    ; ccr = new ( 2*mxlag+1, double)    
    ; ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
    ; ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
    ; x=ispan(-120,120,1)     
    ; wks=gsn_open_wks("png", ""+diri+"indopa")
    ; res=True
    ; plot=gsn_csm_xy(wks,x,ccr,res)
    ; print("r1="+r1)
    ; print("r2="+r2)

; ;;;

;;;趋势 1950-1979 1979-2010
  ipindex_anu_std_rc=regCoef(ipindex_anu_std&time,ipindex_anu_std)
  ip_rc =ipindex_anu_std_rc*ipindex_anu_std&time+ipindex_anu_std_rc@yintercept
   
   
  ; ipindex_anu_std_mk=trend_manken(ipindex_anu_std,False,0)  ;;;;mk检验
  ; print("mk1="+ipindex_anu_std_mk(0))
  ; print("mk2="+ipindex_anu_std_mk(1))

  ipindex_anu_std_rc1=regCoef(ispan(1,23,1),ipindex_anu_std(0:22))
  ip_rc1 =ipindex_anu_std_rc1*ispan(1, 23, 1)+ipindex_anu_std_rc1@yintercept 
  
  ipindex_anu_std_rc2=regCoef(ispan(24,55,1),ipindex_anu_std(23:54))
  ip_rc2=ipindex_anu_std_rc2*ispan(24, 55, 1)+ipindex_anu_std_rc2@yintercept   
;;;
   
   
;;;滑动趋势
  iprc31 = rmtrend(ipindex_anu_std,31,1958)
  iprc21 = rmtrend(ipindex_anu_std,21,1958)
  iprc11 = rmtrend(ipindex_anu_std,11,1958)

;;;
 



;;;画bar图
   wks=gsn_open_wks("pdf", ""+diri+"w定义的印太指数_不去趋势")
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    ; resxy@gsnLeftString="(c)"
    resxy@gsnXYBarChart = True            ; Create bar plot
    ; resxy@tmYLMode="Explicit"
    ; resxy@trYMaxF=2.0
    ; resxy@trYMinF=-2.0
    
    ; resxy@tmYLValues=fspan(-2, 2, 9)
    ; resxy@tmYLLabels=fspan(-2, 2, 9)
    resxy@gsnYRefLine  = 0.   
    resxy@gsnBelowYRefLineBarColors="blue"
    resxy@gsnAboveYRefLineBarColors="red"
    resxy@tiYAxisString="index"
    resxy@tmYLMinorOn=False
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(ipindex_anu_std&time)
    resxy@tmXBValues=ispan(1,55,5)
    resxy@tmXBLabels=ispan(1958,2012,5)
    resxy@tmXBMinorValues=ispan(1,55,1)
    
    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9  
    plot = gsn_csm_xy (wks,ipindex_anu_std&time , ipindex_anu_std,resxy)

;;;将nino3_4指数画上去
 
   res=True
   res@gsLineColor="black"
   res@gsLineThicknessF=2.0
   ; Polyline=gsn_add_polyline(wks, plot, ipindex_anu_std&time, nino3_4_anu_std, res)
     res@gsLineColor="red"
   Polyline_0=gsn_add_polyline(wks, plot, ipindex_anu_std&time, pdo_anu, res)


   res@gsLineColor="blue"
   Polyline1=gsn_add_polyline(wks, plot, ipindex_anu_std&time, ip_rc, res)
    res@gsLineColor="green"
   Polyline2=gsn_add_polyline(wks, plot, ispan(1, 23, 1), ip_rc1, res)
    res@gsLineColor="green"
   Polyline3=gsn_add_polyline(wks, plot, ispan(24, 55, 1), ip_rc2, res)

   draw(plot)
   frame(wks)

exit()

;;画月份变化的情况
  wks2=gsn_open_wks("pdf", ""+diri+"印太指数的月份变化")
  res2=True
  ; res2@xyLabelMode="MarkLines"
  plot2=gsn_csm_xy(wks2,ispan(0, 11, 1), ipindex_season, res2)


;;;画滑动趋势图
  wks3=gsn_open_wks("pdf", ""+diri+"印太指数的滑动趋势")
  res3=True
  res3@gsnFrame=False
  res3@gsnDraw=False
  res3@tmXTOn=False
  res3@tmYROn=False
  res3@gsnYRefLine=0
  res3@gsnYRefLineColor="blue"
  res3@tmXMajorGrid=True
  res3@tmYMajorGrid=True
  res3@tmXMajorGridLineDashPattern=2
  res3@tmXMajorGridLineColor="gray"
  res3@tmXMajorGridThicknessF=0.5
  res3@tmYMajorGridLineDashPattern=2
  res3@tmYMajorGridLineColor="gray"
  res3@tmYMajorGridThicknessF=0.5
  res3@gsnRightString="move window 31yr"
  res3@vpHeightF=0.3
  res3@vpWidthF= 0.7
  plot_rm=new(3,graphic)

  plot_rm(0) = gsn_csm_xy (wks3,iprc31&time , iprc31,res3)
  res3@gsnRightString="move window 21yr"
  plot_rm(1) = gsn_csm_xy (wks3,iprc21&time , iprc21,res3)
  res3@gsnRightString="move window 11yr"
  plot_rm(2) = gsn_csm_xy (wks3,iprc11&time , iprc11,res3)
  
;;;给显著性检验通过的点加标记

  resmak=True
  resmak@gsMarkerColor="red"
  resmak@gsMarkerIndex=3
  
  Polyline3=gsn_add_polymarker(wks3, plot_rm(0), iprc31&time, iprc31@rc_test, resmak)

  Polyline4=gsn_add_polymarker(wks3, plot_rm(1), iprc21&time, iprc21@rc_test, resmak)

  Polyline5=gsn_add_polymarker(wks3, plot_rm(2), iprc11&time, iprc11@rc_test, resmak)



  pres=True
  gsn_panel(wks3, plot_rm, (/3,1/), pres)
 end
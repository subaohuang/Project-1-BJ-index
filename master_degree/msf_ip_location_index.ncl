
 ;;;挑出最大的位置，三维的情况
 ; function malon(X,case)
 ; begin
 ; lon_min=new(dimsizes(X&time),float)
 ; if case.eq."min" then 
 ;     do t=0,dimsizes(X&time)-1
 ;        dims = dimsizes(X(:,:,t))
 ;        x1d = ndtooned(X(:,:,t))      ; convert 2D array to 1D for use in minind
 ;        inds = ind_resolve(minind (x1d), dims)    ; convert 1D array back to 2D 
 ;        ilon =inds(0,1)        ; select the longitude index where the X array is at its' minimum  
 ;        lon_min(t)= X&lon(ilon)   ; insert the longitude index into the lon coordinate variable
 ;     end do
 ; end if 

 ; if case.eq."max" then 
 ;      do t=0,dimsizes(X&time)-1
 ;      dims = dimsizes(X(:,:,t))
 ;      x1d = ndtooned(X(:,:,t))      ; convert 2D array to 1D for use in minind
 ;      inds = ind_resolve(maxind (x1d), dims)    ; convert 1D array back to 2D 
 ;      ilon = inds(0,1)        ; select the longitude index where the X array is at its' minimum  
 ;      lon_min(t) = X&lon(ilon)   ; insert the longitude index into the lon coordinate variable
 ;    end do
 ; end if 
 ; return(lon_min)
 ; end
;;;

function malon(X,case)
 begin
 lon_min=new(dimsizes(X&time),float)
 if case.eq."min" then 
     do t=0,dimsizes(X&time)-1
        ilon =minind(X(:,t))       ; select the longitude index where the X array is at its' minimum  
        lon_min(t)= X&lon(ilon)   ; insert the longitude index into the lon coordinate variable
     end do
 end if 

 if case.eq."max" then 
      do t=0,dimsizes(X&time)-1
        ilon =maxind(X(:,t))       ; select the longitude index where the X array is at its' minimum  
        lon_min(t)= X&lon(ilon)   ; insert the longitude index into the lon coordinate variable
    end do
 end if 
 return(lon_min)
 end

begin
diri="/mnt/f/share/ncep/ncep1/"
diri_pi="/mnt/f/share/era/data/"
fils1=systemfunc("ls "+diri+"msf_time_necp1.nc")
f1=addfile(fils1, "r")

msf_time=f1->msf_time(:,:,0:731)
printVarSummary(msf_time)

; msf_indian=msf_time(:,{40:110},:)
; msf_pacific=msf_time(:,{150:260},:)


plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
printVarSummary(dp)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)
printVarSummary(msfdp)
msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
printVarSummary(msf_mean)
 msf_mean_new=msf_mean(lon|:,time|:)
 msf_indian=msf_mean_new({40:110},:)
 msf_pacific=msf_mean_new({150:260},:)



in_lon=malon(msf_indian,"min")
pa_lon=malon(msf_pacific,"max")

printVarSummary(in_lon)
printVarSummary(pa_lon)

ipindex=pa_lon - in_lon

 opt=diri_pi+"msf_ip_index_1950_2010.ascii"
 opt1=diri_pi+"msf_pa_location_index_1950_2010.ascii"
 opt2=diri_pi+"msf_id_location_index_1950_2010.ascii"

 asciiwrite(opt, sprintf ("%6.2f,", ipindex))
 asciiwrite(opt1, sprintf ("%6.2f,", pa_lon))
 asciiwrite(opt2, sprintf ("%6.2f,", in_lon))



;;;1997年SON的位置与1994年纯IOD事件SON的位置
    ; ipindex_iod=dim_avg_n_Wrap(ipindex(188:190), 0)
    ; ipindex_enino=dim_avg_n_Wrap(ipindex(224:226), 0)

    ipindex_iod=dim_avg_n_Wrap(ipindex(180:191), 0)
    ipindex_enino=dim_avg_n_Wrap(ipindex(216:227), 0)
    
    print("iod="+ipindex_iod)
    print("enino="+ipindex_enino)
;;;---    



;;;计算指数
    ipindex_anu=month_to_annual(ipindex, 1)
    
    ; ipindex_anu=(/dtrend(ipindex_anu, True)/)
    ipindex_anu=ipindex_anu
    ipindex_anu_std=dim_standardize(ipindex_anu, 1)
    ipindex_anu_std!0="time" 
    ipindex_anu_std&time=ispan(1, 61, 1)
    



; ;;;找到最大量值的位置的量值是多少
;   msf_min=dim_min_n_Wrap(msf_indian, (/0,1/))
;   msf_max=dim_min_n_Wrap(msf_pacific, (/0,1/))
;   msfindex=msf_max-msf_min

;   msfindex_anu=month_to_annual(msfindex, 1)
;   msfindex_anu=(/dtrend(msfindex_anu, True)/)

;   msfindex_anu_std=dim_standardize(msfindex_anu, 1)
;   msfindex_anu_std!0="time" 
;   msfindex_anu_std&time=ispan(1, 38, 1)

;   print("msfindex16="+msfindex_anu_std(9))

; ;;;强度除位置的指数
;    index=(/msfindex_anu_std/ipindex_anu_std/)
;    index!0="time"
;    index&time=ispan(1,38,1)
;    ; print(index)
; ;;;   


; ;;;质量加权平均的指数
;     printVarSummary(msf_time)
;     plev=msf_time&level
;     psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
;     ptop=min(plev)
;     dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
;     printVarSummary(dp)
;     msfdp=msf_time(level|:,time|:,lon|:)*dp
;     copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)
;     printVarSummary(msfdp)
;     msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
;     copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
;     printVarSummary(msf_mean)

;     msf_mean_iod=dim_avg_n_Wrap(msf_mean(:,{40:140}), 1)
;     msf_mean_pac=dim_avg_n_Wrap(msf_mean(:,{140:260}), 1)

;     msfindex2=msf_mean_pac-msf_mean_iod 
;     msfindex2_anu=month_to_annual(msfindex2, 1)
;     msfindex2_anu=(/dtrend(msfindex2_anu, True)/)
;     msfindex2_anu_std=dim_standardize(msfindex2_anu, 1)
;     msfindex2_anu_std!0="time" 
;     msfindex2_anu_std&time=ispan(1, 38, 1)
;     print("msfindex2="+msfindex2_anu_std(16))
; ;;;

; ;;;强度除位置的指数
;    index2=(/msfindex2_anu_std/ipindex_anu_std/)
;    index2!0="time"
;    index2&time=ispan(1,38,1)
;     print(index2)
; ;;;   
  

;;;读入nino3.4 index
    nino3_4=asciiread(diri_pi+"Hadley_nino3_4_1950_2010_nodtrend.ascii", 732, "float")
    
   ; nino3_4=(/dtrend(nino3_4, True)/)

    nino3_4_anu=month_to_annual(nino3_4, 1)
    nino3_4_anu_std=dim_standardize(nino3_4_anu, 1)
    nino3_4_anu_std!0="time" 
    nino3_4_anu_std&time=ispan(1, 61, 1)
    
;;;读入pdo index
  fils1=systemfunc("ls /mnt/f/share/era/data/pdo.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->TSAXIS, 1)
  timeind=ind(date.ge.195001.and.date.le.201012)
  pdo=f1->PDO(timeind)
  pdo!0="time"
  printVarSummary(pdo)
  pdo_anu=month_to_annual(pdo, 1)
  pdo_anu&year=ispan(1, 61, 1)


;;;


; ;;;相关性
;     r1=escorc_n(nino3_4_anu_std, msfindex2_anu_std, 0, 0)

;     mxlag    = 120
;     x_Lead_y = esccr(ipindex,nino3_4,mxlag)
;     y_Lead_x = esccr(nino3_4,ipindex,mxlag)    ; switch the order of the series

;     ccr = new ( 2*mxlag+1, double)    
;     ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
;     ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
;     x=ispan(-120,120,1)     
;     wks=gsn_open_wks("png", ""+diri+"indopa")
;     res=True
;     plot=gsn_csm_xy(wks,x,ccr,res)
;     print("r1="+r1)
; ;;;

;;;计算指数随时间的回归
  
  ipindex_anu_std_rc=regCoef(ipindex_anu_std&time,ipindex_anu_std)
  ip_rc =ipindex_anu_std_rc*ipindex_anu_std&time+ipindex_anu_std_rc@yintercept

  ipindex_anu_std_rc1=regCoef(ispan(1,30,1),ipindex_anu_std(0:29))
  ip_rc1 =ipindex_anu_std_rc1*ispan(1, 30, 1)+ipindex_anu_std_rc1@yintercept 

  ipindex_anu_std_rc2=regCoef(ispan(1,30,1),ipindex_anu_std(30:59))
  ip_rc2=ipindex_anu_std_rc2*ispan(1, 30, 1)+ipindex_anu_std_rc2@yintercept   
;;;

;;;画bar图
  wks=gsn_open_wks("pdf", ""+diri_pi+"msf定义的印太指数index")
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
    resxy@tmXBValues=ispan(1,61,5)
    resxy@tmXBLabels=ispan(1950,2010,5)
    resxy@tmXBMinorValues=ispan(1,61,1)
    
    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9  
    plot = gsn_csm_xy (wks,ipindex_anu_std&time , ipindex_anu_std,resxy)



;;;将nino3_4指数画上去
 
   res=True
   res@gsLineColor="black"
   res@gsLineThicknessF=2.0

   Polyline=gsn_add_polyline(wks, plot, ipindex_anu_std&time, nino3_4_anu_std, res)
   res@gsLineColor="blue"
   Polyline1=gsn_add_polyline(wks, plot, ipindex_anu_std&time, ip_rc, res)
   res@gsLineColor="red"

   Polyline_0=gsn_add_polyline(wks, plot, ipindex_anu_std&time, pdo_anu, res)
  



   res@gsLineColor="green"
   Polyline2=gsn_add_polyline(wks, plot, ispan(1, 30, 1), ip_rc1, res)
   
   Polyline3=gsn_add_polyline(wks, plot, ispan(31, 61, 1), ip_rc2, res)


   draw(plot)
   frame(wks)
   
  end
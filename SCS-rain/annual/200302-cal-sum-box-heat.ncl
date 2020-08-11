;;********************************************************
;;计算SST倾向以及SST值与净长波辐射，短波辐射，感热，潜热以及平流和对流,
;;其中热通量的方向都是向上为正，向下为负。。
;;;的回归场，热通量数据是NCEP-NCAR数据，海洋数据是SODA2.2.4,SST是哈哈得来SST数据
;;由于数据之间的精度不同，需要进行插值预处理,SODA精度是0.5×0.5,HADISST是1×1
;;NCEP-NCAR的精度为1.87×1.9
;;s
;;                                         lnshuheng 2020.02.18
;;
 ;; u_lateral_east is wrong ,modifyed in 2020.03.09
;**********************************************************
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
function  wgt_annual(total_month)
begin
    weight1 = ispan(1,12,1)*1.0
    weight1 = (2*weight1-1)/24

    weight2 = ispan(1,12,1)*1.0
    weight2 = (25-2*weight2)/24
    total_annual_test = new(34,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,32,1
      total_annual_test(i+1) = sum(total_month(i*12:(i+1)*12-1)*weight1) + sum(total_month((i+1)*12:(i+2)*12-1)*weight2)
    end do    
return total_annual_test
end    
    


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta_annual = month_to_annual(ssta(:,{lat1:lat2},{lon1:lon2}), 1)
    ; ssta_annual =(/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)

    ;;;sst 倾向值
    dsst = center_finite_diff_n(ssta, 1, False, 0, 0) 
    copy_VarCoords(ssta, dsst)
    dsst_annual = month_to_annual(dsst(:,{lat1:lat2},{lon1:lon2}), 1)
    ; dsst_annual =(/dtrend_msg_n(ispan(1,yearnum,1),dsst_annual,False,True,0)/)
    dsst_annual_mean = dim_avg_n_Wrap(dsst_annual(:,{0:15},{110:155}), (/1,2/))
    ;;;
    delete(date)
    delete(timeind)

    sst_mean = dim_avg_n_Wrap(sst(:,{0:15},{110:155}), (/1,2/))



;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_mlt_mn_1980-2015_chazhi.nc", "r") 
    f7=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_net_heating_mn_1980-2015_chazhi.nc", "r")     


    indyear = (2013-1980+1)*12-1
    u = f2->u(0:indyear,:,{lat1:lat2},{lon1:lon2})
    v = f3->v(0:indyear,:,{lat1:lat2},{lon1:lon2})
    wt = f4->wt(0:indyear,:,{lat1:lat2},{lon1:lon2})
    T = f5->temp(0:indyear,:,{lat1:lat2},{lon1:lon2})
    mlt = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})
    net_heating = f7->net_heating(0:indyear,{lat1:lat2},{lon1:lon2})
   
    depth = T&depth(0:4)
printVarSummary(wt)



;;;calculate lateral  and bottom boundaries of the domain Heat transport, bottom is 45m , lateral is 
  ;;; four horizatal boundaries

    ;;parameter 
     pi = get_pi("float")
     dx = u&lat
     dx = 6378388.*cos(0.0174533*u&lat)
     dx!0 = "lat"
     copy_VarCoords(u(1,1,:,1), dx)

     dy = (pi/180)*1*6378388.
     dz = dpres_plevel(depth, 45.58, 0, 0)
     print(dz)
     print(depth)
     dz= dz(::-1)
   

    
     dx_3d_test = conform_dims(dimsizes(T(:,{0:45.6},{0:15},{110:155})) ,dx({0:15}) , 2)
     dz_3d_test = conform_dims(dimsizes(T(:,{0:45.6},{0:15},{110:155})) ,dz , 1)


     dx_3d_test  = where(ismissing(T(:,{0:45.6},{0:15},{110:155})), 0.0, dx_3d_test)
     dz_3d_test  = where(ismissing(T(:,{0:45.6},{0:15},{110:155})), 0.0, dz_3d_test)
     
     volume = dim_sum_n(dx_3d_test*dz_3d_test*dy, (/1,2,3/))


     Tbox =dim_sum_n_Wrap(T(:,{0:45.6},{0:15},{110:155})*dy*dx_3d_test*dz_3d_test, (/1,2,3/))/volume
     copy_VarCoords(T(:,1,1,1), Tbox)

     ; volume = sum(45*dy*45.58*dx({0:15})) 

     dT = center_finite_diff_n(Tbox, 1, False, 0, 0)
     copy_VarCoords(Tbox, dT)
     dT_annual  = month_to_annual(dT, 1) 
    
     Tbox_annual = month_to_annual(Tbox, 1)


  
     dT_test = center_finite_diff_n(Tbox_annual, 1, False, 0, 0)
     dT_wgt = dT_test
     dT_wgt(1:33) = Tbox_annual(1:33) - Tbox_annual(0:32)
   
     

    ;; lateral heat transport north and south v*(T -Tbox)dxdz,south: 0N ,110-155,north: 15N ,110-155
       x1 = 0
       x2 = 15
       y1  = 110.5
       y2 = 155
       z1  = 45.6
    
       dim =dimsizes(v(:,{0:z1},{x1},{y1:y2}))
     

       dz_3d =conform_dims(dim ,dz , 1)
       Tbox_3d = conform_dims(dim ,Tbox , 0)
       printVarSummary(dz_3d)  
       
       
       v_lateral_south = dim_sum_n_Wrap(v(:,{0:z1},{x1},{y1:y2})*(T(:,{0:z1},{x1},{y1:y2}) - Tbox_3d)*dx({0})*dz_3d, (/1,2/))/volume
       v_lateral_north = dim_sum_n_Wrap(v(:,{0:z1},{x2},{y1:y2})*(T(:,{0:z1},{x2},{y1:y2}) - Tbox_3d)*dx({15})*dz_3d, (/1,2/))/volume

    ;; lateral heat transport west and east u*(T -Tbox)dydz,west: 0N ,110-155,north: 15N ,110-155
       dim := dimsizes(v(:,{0:z1},{x1:x2},{y1}))
       ; dx_3d =conform_dims(dim ,dx({0:15}) , 2)
       dz_3d := conform_dims(dim ,dz , 1)
       Tbox_3d :=conform_dims(dim, Tbox, 0)

       u_lateral_west = dim_sum_n_Wrap(u(:,{0:z1},{x1:x2},{y1})*(T(:,{0:z1},{x1:x2},{y1}) - Tbox_3d)*dy*dz_3d, (/1,2/))/volume
       u_lateral_east = dim_sum_n_Wrap(u(:,{0:z1},{x1:x2},{y2})*(T(:,{0:z1},{x1:x2},{y2}) - Tbox_3d)*dy*dz_3d, (/1,2/))/volume
    ;;;
    

    ;;; lateral heat transport vertical w*(T -Tbox)dxdy 
       
       dim:= dimsizes(wt(:,4,{x1:x2},{y1:y2}))
       dx_3d :=conform_dims(dim ,dx({0:15}) , 1)
       Tbox_3d :=conform_dims(dim, Tbox, 0)       
  
       w_bottom = dim_sum_n_Wrap(wt(:,4,{x1:x2},{y1:y2})*(T(:,4,{x1:x2},{y1:y2}) - Tbox_3d)*dx_3d*dy, (/1,2/))/volume
    ;;;
    
    
    ;;;;v_lateral_south west and east to define effect from ITF     
        
       dz_3d1 =conform_dims(dimsizes(v(:,{0:z1},{x1},{y1:130})), dz , 1)
       Tbox_3d1 = conform_dims(dimsizes(v(:,{0:z1},{x1},{y1:130})) ,Tbox , 0)
         
    
       v_lateral_south_west = dim_sum_n_Wrap(v(:,{0:z1},{x1},{y1:130})*(T(:,{0:z1},{x1},{y1:130}) - Tbox_3d1)*dx({0})*dz_3d1, (/1,2/))/volume

       dz_3d2 =conform_dims(dimsizes(v(:,{0:z1},{x1},{130:y2})), dz , 1)
       Tbox_3d2 = conform_dims(dimsizes(v(:,{0:z1},{x1},{130:y2})) ,Tbox , 0)
         
       v_lateral_south_east = dim_sum_n_Wrap(v(:,{0:z1},{x1},{130:y2})*(T(:,{0:z1},{x1},{130:y2}) - Tbox_3d2)*dx({0})*dz_3d2, (/1,2/))/volume
    

        dz_3d3 =conform_dims(dimsizes(v(:,{0:z1},{7.5:x2},{y1})), dz , 1)
        Tbox_3d3 = conform_dims(dimsizes(v(:,{0:z1},{7.5:x2},{y1})) ,Tbox , 0)
       
        u_lateral_west_north = dim_sum_n_Wrap(u(:,{0:z1},{7.5:x2},{y1})*(T(:,{0:z1},{7.5:x2},{y1}) - Tbox_3d3)*dy*dz_3d3, (/1,2/))/volume
        
        dz_3d4 =conform_dims(dimsizes(v(:,{0:z1},{x1:7.5},{y1})), dz , 1)
        Tbox_3d4 = conform_dims(dimsizes(v(:,{0:z1},{x1:7.5},{y1})) ,Tbox , 0)

        u_lateral_west_south = dim_sum_n_Wrap(u(:,{0:z1},{x1:7.5},{y1})*(T(:,{0:z1},{x1:7.5},{y1}) - Tbox_3d4)*dy*dz_3d4, (/1,2/))/volume

    ;;;




    
    ;;;平流项
       hadv_mean = v_lateral_north +v_lateral_south + u_lateral_west + u_lateral_east
       vadv_mean = w_bottom

       copy_VarCoords(Tbox, hadv_mean)
       copy_VarCoords(Tbox, vadv_mean)

       ;;convert to ℃/month
       hadv_mean = 30*86400*hadv_mean
       vadv_mean = 30*86400*vadv_mean
       v_lateral_north = 30*86400*v_lateral_north
       v_lateral_south = 30*86400*v_lateral_south
       u_lateral_east  = 30*86400*u_lateral_east
       u_lateral_west = 30*86400*u_lateral_west

       v_lateral_south_west = 30*86400*v_lateral_south_west
       v_lateral_south_east = 30*86400*v_lateral_south_east

       u_lateral_west_north = 30*86400*u_lateral_west_north
       u_lateral_west_south = 30*86400*u_lateral_west_south
       

    ;;;

    


    ;;;年平均
       ; hadv_mean = (/rmMonAnnCycTLL(hadv_mean)/)
       ; vadv_mean = (/rmMonAnnCycTLL(vadv_mean)/)
       ; hadv_mean_annual = month_to_annual(hadv_mean, 0)
       ; vadv_mean_annual = month_to_annual(vadv_mean, 0)

    ;;;权重平均 from paper warming of the upper Equatorial indian ocean and changes
       
        hadv_mean_annual = wgt_annual(hadv_mean)
        vadv_mean_annual = wgt_annual(vadv_mean)

        v_lateral_south_annual = wgt_annual(v_lateral_south)
        v_lateral_north_annual = wgt_annual(v_lateral_north)

        u_lateral_east_annual = wgt_annual(u_lateral_east)
        u_lateral_west_annual = wgt_annual(u_lateral_west)

        v_lateral_south_west_annual = wgt_annual(v_lateral_south_west)
        v_lateral_south_east_annual = wgt_annual(v_lateral_south_east)

        u_lateral_west_south_annual = wgt_annual(u_lateral_west_south)
        u_lateral_west_north_annual = wgt_annual(u_lateral_west_north)
    ;;;

    ;;;
        total_ocean_annual = hadv_mean_annual
        total_ocean_annual = hadv_mean_annual + vadv_mean_annual
    ;;;

;;;read NCEP-NCAR heat flux data
  
    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/NCEP-NCAR-heatflux-1degree.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)
  
    
    lhtfl = f3->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f3->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f3->nlwrs(timeind,{lat1:lat2},{lon1:lon2})


    ;;;nswrs_pen 
     R = 0.58
     L1 = 0.35
     L2 = 23
     nswrs_pen = nswrs
     ; nswrs_pen =  nswrs*(R*exp(-1*mlt/L1)+ (1-R)*exp(-1*mlt/L2))
     mo = (R*exp(-1*45.6/L1)+ (1-R)*exp(-1*45.6/L2))
     nswrs_pen = mo*nswrs

  
    ;;;

    ;;;convert to units degree ℃/year
      ; density = 1025  ;;kg/m3
      ; cp  = 3850         ;;j/kg℃
      ; coef  = density*cp
      coef = 4.088*10^6

      lhtfl = 30*86400*lhtfl/coef
      shtfl = 30*86400*shtfl/coef
      nswrs = 30*86400*nswrs/coef
      nlwrs = 30*86400*nlwrs/coef
      nswrs_pen = 30*86400*nswrs_pen/coef
    ;;;

    ; lhtfl  = (/rmMonAnnCycTLL(lhtfl)/)
    ; shtfl  = (/rmMonAnnCycTLL(shtfl)/)
    ; nswrs  = (/rmMonAnnCycTLL(nswrs)/)
    ; nlwrs  = (/rmMonAnnCycTLL(nlwrs)/)
    ; nswrs_pen  = (/rmMonAnnCycTLL(nswrs_pen)/)

    ;;;转化成方向向下

    lhtfl = (/-1*lhtfl/)
    shtfl = (/-1*shtfl/)
    nswrs = (/-1*nswrs/)
    nlwrs = (/-1*nlwrs/)
    nswrs_pen = (/-1*nswrs_pen/)

    Q_net = lhtfl
    Q_net = lhtfl + shtfl + nswrs + nlwrs - nswrs_pen
    copy_VarCoords(lhtfl, Q_net)
    ;;;

    ;;; Heat transport 

    printVarSummary(Q_net)
    Q_net&lat@units = "degrees_north"
    Q_net&lon@units = "degrees_east"


    
    Q_net(:,{x1:x2},{y1:y2}) = where(ismissing(sst(:,{0:15},{110:155})),Q_net@missing_value, Q_net(:,{x1:x2},{y1:y2}))  ;;;去掉陆地上的辐射 
   
    Q_net_tr = dim_sum_n_Wrap(Q_net(:,{x1:x2},{y1:y2})*dx_3d*dy, (/1,2/))/volume
    copy_VarCoords(Tbox, Q_net_tr)

    ;;;
    
    ; Q_net_tr_annual = month_to_annual(Q_net_tr, 0)   ;;;这里是累加
; 
     Q_net_tr_annual = wgt_annual(Q_net_tr)

;;ocena net_heating 
        net_heating = 30*86400*net_heating/coef
        net_heating_tr = dim_sum_n_Wrap(net_heating(:,{x1:x2},{y1:y2})*dx_3d*dy, (/1,2/))/volume
        copy_VarCoords(Tbox, net_heating_tr)


        net_heating_annual = wgt_annual(net_heating_tr)

  
;;;total 

       total_month = Q_net_tr
       total_month = Q_net_tr + hadv_mean + vadv_mean
       
       total_annual  = Q_net_tr_annual + hadv_mean_annual + vadv_mean_annual
;;;


;;;trend

   rc_total  = regCoef_n(ispan(1, yearnum, 1),total_annual, 0, 0)
   rc_atm = regCoef_n(ispan(1, yearnum, 1),Q_net_tr_annual, 0, 0)
   rc_ocean = regCoef_n(ispan(1, yearnum, 1),total_ocean_annual, 0, 0)
   rc_dt = regCoef_n(ispan(1, yearnum, 1),dT_wgt, 0, 0)

   rc_hadv = regCoef_n(ispan(1, yearnum, 1),hadv_mean_annual, 0, 0)
   rc_vadv = regCoef_n(ispan(1, yearnum, 1),vadv_mean_annual, 0, 0)

   rc_u_lateral_west = regCoef_n(ispan(1, yearnum, 1),u_lateral_west_annual, 0, 0)
   rc_u_lateral_east = regCoef_n(ispan(1, yearnum, 1),u_lateral_east_annual, 0, 0)

   rc_v_lateral_south = regCoef_n(ispan(1, yearnum, 1),v_lateral_south_annual, 0, 0)
   rc_v_lateral_north = regCoef_n(ispan(1, yearnum, 1),v_lateral_north_annual, 0, 0)

  

   rc_v_lateral_south_west = regCoef_n(ispan(1, yearnum, 1),v_lateral_south_west_annual, 0, 0)
   rc_v_lateral_south_east = regCoef_n(ispan(1, yearnum, 1),v_lateral_south_east_annual, 0, 0)


   rc_u_lateral_west_south = regCoef_n(ispan(1, yearnum, 1),u_lateral_west_south_annual, 0, 0)
   rc_u_lateral_west_north = regCoef_n(ispan(1, yearnum, 1),u_lateral_west_north_annual, 0, 0)





   total_trend = rc_total*ispan(1, yearnum, 1)+ rc_total@yintercept 
   atm_trend = rc_atm*ispan(1, yearnum, 1)+ rc_atm@yintercept 
   ocean_trend = rc_ocean*ispan(1, yearnum, 1)+ rc_ocean@yintercept 
   
   rc_total = rc_total*yearnum
   rc_atm = rc_atm*yearnum
   rc_ocean = rc_ocean*yearnum
   rc_dt = rc_dt*yearnum

   rc_hadv  = rc_hadv*yearnum
   rc_vadv  = rc_vadv*yearnum

   rc_u_lateral_west = rc_u_lateral_west*yearnum
   rc_u_lateral_east = rc_u_lateral_east*yearnum
    
   rc_v_lateral_south = rc_v_lateral_south*yearnum
   rc_v_lateral_north = rc_v_lateral_north*yearnum

    rc_v_lateral_south_west = rc_v_lateral_south_west*yearnum
    rc_v_lateral_south_east = rc_v_lateral_south_east*yearnum

    rc_u_lateral_west_south = rc_u_lateral_west_south*yearnum
    rc_u_lateral_west_north = rc_u_lateral_west_north*yearnum
   
   



   print("total_rc="+rc_total)
   print("atm_rc="+rc_atm)
   print("ocean_rc="+rc_ocean)
   print("dt="+rc_dt)

   print("hadv_rc="+rc_hadv)
   print("vadv_rc="+rc_vadv)

   print("rc_u_lateral_east="+rc_u_lateral_east)
   print("rc_u_lateral_west="+rc_u_lateral_west)

   print("rc_v_lateral_south="+rc_v_lateral_south)
   print("rc_v_lateral_north="+rc_v_lateral_north)

   print("rc_v_lateral_south_west="+rc_v_lateral_south_west)
   print("rc_v_lateral_south_east="+rc_v_lateral_south_east)

   print("rc_u_lateral_west_south="+rc_u_lateral_west_south)
   print("rc_u_lateral_west_north="+rc_u_lateral_west_north)
   
   
  
   
   rc_Tbox = regCoef_n(ispan(1, yearnum, 1),Tbox_annual, 0, 0)

   print("rc_Tbox="+rc_Tbox)



;;;

;;;;两段平均 80-96 97-13
time = ispan(1980, 2013, 1)*1.0
m = new(34,"float")

Q_net_tr_annual_1 = time(0:16)
Q_net_tr_annual_1= dim_avg_n_Wrap(Q_net_tr_annual(0:16), 0)
Q_net_tr_annual_2  =time(17:)
Q_net_tr_annual_2= dim_avg_n_Wrap(Q_net_tr_annual(17:), 0)

total_ocean_annual_1 = time(0:16)
total_ocean_annual_1= dim_avg_n_Wrap(total_ocean_annual(0:16), 0)
total_ocean_annual_2 = time(17:)
total_ocean_annual_2= dim_avg_n_Wrap(total_ocean_annual(17:), 0)

total_annual_1 = time(0:16)
total_annual_1=  dim_avg_n_Wrap(total_annual(0:16), 0)
total_annual_2 = time(17:)
total_annual_2=  dim_avg_n_Wrap(total_annual(17:), 0)

dT_wgt_mean_1 = time(0:16)
dT_wgt_mean_1  =dim_avg_n_Wrap(dT_wgt(0:16), 0)
dT_wgt_mean_2 = time(17:)
dT_wgt_mean_2  =dim_avg_n_Wrap(dT_wgt(17:), 0)

print("dT_wgt_mean_1="+dim_avg_n_Wrap(dT_wgt(0:16), 0))
print("dT_wgt_mean_2="+dim_avg_n_Wrap(dT_wgt(17:), 0))

Q_net_diff  = dim_avg_n_Wrap(Q_net_tr_annual(17:), 0) - dim_avg_n_Wrap(Q_net_tr_annual(0:16), 0)
ocean_diff = dim_avg_n_Wrap(total_ocean_annual(17:), 0) - dim_avg_n_Wrap(total_ocean_annual(0:16), 0)
dT_wgt_diff = dim_avg_n_Wrap(dT_wgt(17:), 0) - dim_avg_n_Wrap(dT_wgt(0:16), 0)
ocean_vadv_diff = dim_avg_n_Wrap(vadv_mean_annual(17:), 0) - dim_avg_n_Wrap(vadv_mean_annual(0:16), 0)
; ocean_vadv_diff = dim_avg_n_Wrap(vadv_mean(12*17:), 0) - dim_avg_n_Wrap(vadv_mean(0:12*17-1), 0)
; ocean_hadv_diff = dim_avg_n_Wrap(hadv_mean(12*17:), 0) - dim_avg_n_Wrap(hadv_mean(0:12*17-1), 0)
ocean_hadv_diff = dim_avg_n_Wrap(hadv_mean_annual(17:), 0) - dim_avg_n_Wrap(hadv_mean_annual(0:16), 0)


print("Q_net_diff="+Q_net_diff)
print("ocean_diff="+ocean_diff)
print("dT_wgt_diff="+dT_wgt_diff)
print("ocean_vadv_diff="+ocean_vadv_diff)
print("ocean_hadv_diff="+ocean_hadv_diff)





;;;long-term mean

  Q_net_tr_mean = dim_avg_n_Wrap(Q_net_tr_annual, 0)
  
  total_ocean_mean = dim_avg_n_Wrap(total_ocean_annual, 0)
  dT_wgt_mean = dim_avg_n_Wrap(dT_wgt, 0)

  print("Q_net_tr_mean="+Q_net_tr_mean)
  print("total_ocean_mean="+total_ocean_mean)
  print("dT_wgt_mean="+dT_wgt_mean)

  hadv_clm = dim_avg_n_Wrap(hadv_mean_annual, 0) 
  vadv_clm = dim_avg_n_Wrap(vadv_mean_annual, 0)

  print("hadv_clm="+hadv_clm)
  print("vadv_clm="+vadv_clm)




;;;



     
    

       

   


print("*********开始画图:时间序列********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "all-heat-flux-contribute2BOXSST-wgt-all-component")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 2013
  resxy@trXMinF = 1979
  resxy@tmXBValues = ispan(1979,2013,5)
  resxy@tmXBLabels = ispan(1979,2013,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(1979,2013,1)
  
  ;  resxy@tmYLMode  =  "Explicit"
  ; resxy@tmYLValues = ispan(-2, 10, 1)
  ; resxy@tmYLLabels = ispan(-2,10,1)
  ; resxy@trYMinF = -2
  ; resxy@trYMaxF = 10

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6


  resxy@xyLineThicknesses =  (/ 5,3,2,2,2,2,2/)          ; make second line thicker
  
  resxy@xyLineColors       =  (/"black","gray","blue","red","green","pink"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0,0,0,0,0,0/)
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines","lines"/)
  


  data_all = new((/6,34/),"float")
  data_all(0,:)=dT_wgt
  data_all(1,:)=total_annual
  data_all(2,:)=Q_net_tr_annual
  data_all(3,:)=total_ocean_annual
  data_all(4,:)=hadv_mean_annual
  data_all(5,:)=vadv_mean_annual
  

  ; data_all(6,:)=total_trend
  ; data_all(7,:)=atm_trend
  ; data_all(8,:)=ocean_trend


  plot = gsn_csm_xy(wks,time,data_all, resxy)


  ;;添加直线

  res= True
  res@gsLineThicknessF = 5.0            ; line thickness   
  res@gsLineColor      = "yellow"         ; color of lines
  res@gsLineDashPattern=15

  plot2 = gsn_add_polyline(wks, plot,time(0:16), Q_net_tr_annual_1, res)
  plot3 = gsn_add_polyline(wks, plot,time(0:16), total_ocean_annual_1, res)
 
  plot5 = gsn_add_polyline(wks, plot,time(17:), Q_net_tr_annual_2, res)
  plot6 = gsn_add_polyline(wks, plot,time(17:), total_ocean_annual_2, res)
  
  res@gsLineColor      = "red"         ; color of lines
  
   plot4 = gsn_add_polyline(wks, plot,time(0:16), dT_wgt_mean_1, res)
   plot7 = gsn_add_polyline(wks, plot,time(17:), dT_wgt_mean_2, res)

  





  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       = resxy@xyLineColors
  lineres@lgLineThicknesses = (/ 5,5,5,5,5,5/)
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = resxy@xyDashPatterns


  textres = True
  textres@lgLabels = (/"dT","total","Q_net","Q_ocean","Horizontal advection","Vertical advection"/)     ;加线的标注

  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)
  


  end

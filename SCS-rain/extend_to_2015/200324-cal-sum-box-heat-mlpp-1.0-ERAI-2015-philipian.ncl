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

    dim  = dimsizes(total_month)

    total_annual_test = new(dim(0)/12,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,(dim(0)/12)-2,1
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
    endmon = 201512

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
   lon2 = 170

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

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
    dsst_annual_mean = dim_avg_n_Wrap(dsst_annual(:,{0:15},{125:160}), (/1,2/))
    ;;;
    delete(date)
    delete(timeind)

    sst_mean = dim_avg_n_Wrap(sst(:,{0:15},{125:160}), (/1,2/))



;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    ;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")

    f6=addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_mld_density_criterion_0.8_2015.nc", "r")  

    f7=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_taux_mn_1980-2015_chazhi.nc", "r")   
    f8=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_tauy_mn_1980-2015_chazhi.nc", "r")   

    


    indyear = (endyear -startyear +1)*12-1
    u_r = f2->u(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    v_r = f3->v(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    wt_r = f4->wt(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    T_r = f5->temp(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    mlt = f6->mld(0:indyear,:,:)


    taux_r = f7->taux(0:indyear,{lat1:lat2},{lon1:lon2})
    tauy_r = f8->tauy(0:indyear,{lat1:lat2},{lon1:lon2})
    
    u = tofloat(u_r)
    v = tofloat(v_r)
    wt = tofloat(wt_r)
    T = tofloat(T_r)
    

    taux = tofloat(taux_r)
    tauy = tofloat(tauy_r)
    
    copy_VarCoords(taux_r, taux) 
    copy_VarCoords(tauy_r, tauy)

    copy_VarCoords(u_r, u) 
    copy_VarCoords(v_r, v)
    copy_VarCoords(wt_r, wt)
    copy_VarCoords(T_r, T)

    delete(u_r)
    delete(v_r)
    delete(wt_r)
    delete(T_r)


    delete(taux_r)
    delete(tauy_r)

    depth = T&depth
    printVarSummary(wt)
;;;

;;;计算EKman速度
 
    lat = taux&lat
    pi = get_pi("float")
    omega = 2*pi/86400
    f = 2*omega*sin(pi*lat/180)
    density = 1025
    f_3d = conform_dims(dimsizes(taux), f, 1)
    
    ue = tauy
    ve = taux  
    ue  = (/tauy/(density*f_3d)/)   ;;;EKMAN
    ve = (/-1*taux/(density*f_3d)/)   ;;;Ekman
    vr = uv2vr_cfd((/taux/(density*f_3d)/) ,(/tauy/(density*f_3d)/) ,taux&lat,taux&lon,2)


;;;calculate lateral  and bottom boundaries of the domain Heat transport, bottom is 45m , lateral is 
  ;;; four horizatal boundaries

    ;;parameter 
     pi = get_pi("float")
     dx = u&lat
     dx = 6378388.*cos(0.0184533*u&lat)
     dx!0 = "lat"
     copy_VarCoords(u(1,1,:,1), dx)

     dy = (pi/180)*1*6378388.
 

      dz_3d_test = new(dimsizes(u(:,:,{-5:20},{100:170})),float)
  

      printVarSummary(dz_3d_test)
      printVarSummary(mlt(:,{-5:20},:))
   
      do i = 0,dimsizes(mlt&time)-1
        dz_3d_test(i,:,:,:) = dpres_plevel(depth, mlt(i,{-5:20},:), 0, 0)
      end do 
      copy_VarCoords(u(:,:,{-5:20},{100:170}), dz_3d_test)
      printVarSummary(dz_3d_test)
      ; print(mlt(227,{10},{120}))
    
  
 

    
     dx_3d_test = conform_dims(dimsizes(T(:,:,{0:15},{125:160})) ,dx({0:15}) , 2)


     dx_3d_test  = where(ismissing(T(:,:,{0:15},{125:160})), 0.0, dx_3d_test)
     dz_3d_test  = where(ismissing(T(:,:,{-5:20},{100:170})), 0.0, dz_3d_test)
  
   
     volume = dim_sum_n(dx_3d_test*dz_3d_test(:,:,{0:15},{125:160})*dy, (/1,2,3/))
    

     Tbox =dim_sum_n_Wrap(T(:,:,{0:15},{125:160})*dy*dx_3d_test*dz_3d_test(:,:,{0:15},{125:160}), (/1,2,3/))/volume
     copy_VarCoords(T(:,1,1,1), Tbox)

     ; volume = sum(45*dy*45.58*dx({0:15})) 

     dT = center_finite_diff_n(Tbox, 1, False, 0, 0)
     copy_VarCoords(Tbox, dT)
     dT_annual  = month_to_annual(dT, 1) 
    
     Tbox_annual = month_to_annual(Tbox, 1)


  
     dT_test = center_finite_diff_n(Tbox_annual, 1, False, 0, 0)
     dT_wgt = dT_test
     dT_wgt(1:35) = Tbox_annual(1:35) - Tbox_annual(0:34)   ;;;2015 - 1980
   
     

    ;; lateral heat transport north and south v*(T -Tbox)dxdz,south: 0N ,110-155,north: 15N ,110-155
       x1 = 0
       x2 = 15
       y1  = 125.5
       y2 = 160
       z1  = 45.6
    
       dim =dimsizes(v(:,:,{x1},{y1:y2}))

       dz_3d1 = dz_3d_test(:,:,{x1},{y1:y2})
       dz_3d2 = dz_3d_test(:,:,{x2},{y1:y2})
       Tbox_3d = conform_dims(dim ,Tbox , 0)
      
       
       
       v_lateral_south = dim_sum_n_Wrap(v(:,:,{x1},{y1:y2})*(T(:,:,{x1},{y1:y2}) - Tbox_3d)*dx({0})*dz_3d1, (/1,2/))/volume
       v_lateral_north = dim_sum_n_Wrap(v(:,:,{x2},{y1:y2})*(T(:,:,{x2},{y1:y2}) - Tbox_3d)*dx({15})*dz_3d2, (/1,2/))/volume

    ;; lateral heat transport west and east u*(T -Tbox)dydz,west: 0N ,110-155,north: 15N ,110-155
       dim := dimsizes(v(:,:,{x1:x2},{y1}))
       ; dx_3d =conform_dims(dim ,dx({0:15}) , 2)
       dz_3d1 := dz_3d_test(:,:,{x1:x2},{y1})
       dz_3d2 := dz_3d_test(:,:,{x1:x2},{y2})
       Tbox_3d :=conform_dims(dim, Tbox, 0)

       u_lateral_west = dim_sum_n_Wrap(u(:,:,{x1:x2},{y1})*(T(:,:,{x1:x2},{y1}) - Tbox_3d)*dy*dz_3d1, (/1,2/))/volume
       u_lateral_east = dim_sum_n_Wrap(u(:,:,{x1:x2},{y2})*(T(:,:,{x1:x2},{y2}) - Tbox_3d)*dy*dz_3d2, (/1,2/))/volume
    ;;;
    

    ;;; lateral heat transport vertical w*(T -Tbox)dxdy 
       
       dim:= dimsizes(wt(:,4,{x1:x2},{y1:y2}))
       dx_3d :=conform_dims(dim ,dx({0:15}) , 1)

       Tbox_3d :=conform_dims(dim, Tbox, 0)   

       wtT = wt(:,4,{x1:x2},{y1:y2})
       mlt_scs = mlt(:,{x1:x2},{y1:y2})
       printVarSummary(mlt_scs)
      
       mlt_scs = where(ismissing(mlt_scs), avg(mlt_scs), mlt_scs)
       mlt_scs = where(mlt_scs.lt.5, 5, mlt_scs)
    ;;;;;  



    ;;;;;
        ; print("start triple loops")
        ; printMinMax(mlt_scs, 1)
        ; ;;;modifyed by shuhenglin in 20200823
        ; x1 = 0
        ; x2 = 20
        ; y1  = 110.5
        ; y2 = 170
        ; ;;;;;
        ; wt_mlt = wt(:,4,{x1:x2},{y1:y2})
        ; T_mlt = T(:,4,{x1:x2},{y1:y2})
        ; T_mlt_diff = T(:,4,{x1:x2},{y1:y2})
        ; wt_mlt_dh  = wt(:,4,{x1:x2},{y1:y2})

       ; exit()
       ; ;;;将wt进行垂直方向差值，为了得到更加精确的结果
       ;     print("开始差值")
       ;      xi  = wt&depth 
       ;      xo = ispan(5,330,5)

       ;      wt_chazhi = linint1_n_Wrap(xi,wt(:,:,{x1:x2},{y1:y2}),False,xo,0,1)
       ;      T_chazhi = linint1_n_Wrap(xi,T(:,:,{x1:x2},{y1:y2}),False,xo,0,1)

       ;      dhdt = center_finite_diff_n(mlt, 30*86400, False, 0, 0)
       ;      dhdt!0 = "time"
       ;      dhdt!1 = "lat"
       ;      dhdt!2 = "lon"
       ;      copy_VarCoords(mlt, dhdt)
       ;      printVarSummary(dhdt)
      
       ;      delete(wt)
       ;      delete(T)
       ;    printVarSummary(wt_chazhi)
        
       ;   ; ;    ;;mlt变化引起温度的变化，可以不考虑 

       ;   do i = 0,dim(0)-1
       ;      do j = 0,dim(1)-1
       ;         do k = 0,dim(2)-1
       ;            ; wtT(i,:,:) = (wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2})+dhdt(i,{x1:x2},{y1:y2}))*(T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) - Tbox(i))
       ;             ; wtT(i,:,:) = (wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}))*(T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) - Tbox(i))
       ;             ; print(mlt_scs(i,j,k))
       ;             wt_mlt(i,:,:) = wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2})

       ;             T_mlt(i,:,:) = T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) 
       ;             T_mlt_diff(i,:,:) = T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) - Tbox(i)
       ;             wt_mlt_dh(i,:,:) = wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) + dhdt(i,{x1:x2},{y1:y2})
       ;         end do 
       ;        end do
       ;   end do
       ;    ; wt_mlt_dh = where(wt_mlt_dh.gt.0, wt_mlt_dh, wt_mlt_dh@_FillValue)
       ;    ;;;;写nc文件,简单的方法  
       ;     system("/bin/rm -f SODA3.4.2_wt_T_in_MLD_density_defined_0.8_2015.nc")   ; remove any pre-existing file
            
       ;    ncdf = addfile("SODA3.4.2_wt_T_in_MLD_density_defined_0.8_2015.nc" ,"c")  ; open output netCDF file

       ;    ;===================================================================
       ;    ; create global attributes of the file (optional)
       ;    ;===================================================================
       ;       fAtt               = True            ; assign file attributes
       ;       fAtt@title         = "SODA3.4.2 wt and T in mld,mld is defined by 1.0 density criterion"
       ;       fAtt@Conventions   = "None"
       ;       fAtt@creation_date = systemfunc ("date")
       ;       fileattdef( ncdf, fAtt )            ; copy file attributes

       ;    ;===================================================================
       ;    ; make time an UNLIMITED dimension; recommended  for most applications
       ;    ;===================================================================
       ;       filedimdef(ncdf,"time",-1,True) 

       ;    ;===================================================================
       ;    ; output variables directly; NCL will call appropriate functions
       ;    ; to write the meta data associated with each variable
       ;    ;===================================================================
       ;       ncdf->wt_mlt  = wt_mlt                          ; 4D      
       ;       ncdf->wt_mlt_dh  = wt_mlt_dh                          ; 4D            
       ;       ncdf->T_mlt = T_mlt  
       ;       ncdf->T_mlt_diff = T_mlt_diff                    ; 3D 

       ;    ;;;;
       ;    exit()

;;;;;;
       f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain//extend_to_2015/SODA3.4.2_wt_T_in_MLD_density_defined_0.8_2015.nc", "r")
       ;;;
         x1 = 0
         x2 = 15
         y1  = 125.5
         y2 = 160
       ;;;
       
       wt_mlt = f1->wt_mlt(:,{x1:x2},{y1:y2})
       printVarSummary(wt_mlt)
       ; wt_mlt = where(wt_mlt.le.0, 0, wt_mlt)

       T_mlt_diff = f1->T_mlt_diff(:,{x1:x2},{y1:y2})


       wtT = mlt_scs
       wtT = wt_mlt*T_mlt_diff


;  ;;;;;;
  
       w_bottom = dim_sum_n_Wrap(wtT*dx_3d*dy, (/1,2/))/volume
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

       ; v_lateral_south_west = 30*86400*v_lateral_south_west
       ; v_lateral_south_east = 30*86400*v_lateral_south_east

       ; u_lateral_west_north = 30*86400*u_lateral_west_north
       ; u_lateral_west_south = 30*86400*u_lateral_west_south
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

        ; v_lateral_south_west_annual = wgt_annual(v_lateral_south_west)
        ; v_lateral_south_east_annual = wgt_annual(v_lateral_south_east)

        ; u_lateral_west_south_annual = wgt_annual(u_lateral_west_south)
        ; u_lateral_west_north_annual = wgt_annual(u_lateral_west_north)
       

    ;;;

    ;;;
        total_ocean_annual = hadv_mean_annual
        total_ocean_annual = hadv_mean_annual + vadv_mean_annual
    ;;;

;;;read NCEP-NCAR heat flux data
  
   
;;;read NCEP-NCAR heat flux data
  
    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/annual/ERAI-heatflux-1degree.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.198001)
  
    ;;;region
   lat1 = -10
   lat2 = 20
   lon1 = 100
   lon2 = 170
   ;;;


    lhtfl = f3->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f3->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f3->nlwrs(timeind,{lat1:lat2},{lon1:lon2})


    lhtfl = (/lhtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    shtfl = (/shtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nswrs = (/nswrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nlwrs = (/nlwrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2



    ;;;nswrs_pen 
     R = 0.58
     L1 = 0.35
     L2 = 23
     nswrs_pen = nswrs
     mo = (R*exp(-1*mlt/L1)+ (1-R)*exp(-1*mlt/L2)) 
     ; mo = (1-R)*exp(-1*mlt/20)   ;;;;follow by paper
     nswrs_pen = mo*nswrs
     ; print(mo)
     ; exit()
     ; nswrs_pen =  nswrs*0.05


     ; nswrs_pen = (R*exp(-1*45.6/L1)+ (1-R)*exp(-1*45.6/L2))*nswrs
     ; nswrs_pen = 0.05*nswrs
     


  
    ;;;

    ;;;convert to units degree ℃/yeara
      density = 1026  ;;kg/m3
      cp  = 3993         ;;j/kg℃
      coef  = density*cp
      ; coef = 4.088*10^6

      lhtfl = 30*86400*lhtfl/coef
      shtfl = 30*86400*shtfl/coef
      nswrs = 30*86400*nswrs/coef
      nlwrs = 30*86400*nlwrs/coef
      nswrs_pen = 30*86400*nswrs_pen/coef
    ;;;

    Q_net = lhtfl
    Q_net = lhtfl + shtfl + nswrs + nlwrs - nswrs_pen
    copy_VarCoords(lhtfl, Q_net)
    ;;;

    ;;; Heat transport 

    printVarSummary(Q_net)
    Q_net&lat@units = "degrees_north"
    Q_net&lon@units = "degrees_east"


    
    Q_net(:,{x1:x2},{y1:y2}) = where(ismissing(sst(:,{0:15},{125:160})),Q_net@missing_value, Q_net(:,{x1:x2},{y1:y2}))  ;;;去掉陆地上的辐射 
   
    Q_net_tr = dim_sum_n_Wrap(Q_net(:,{x1:x2},{y1:y2})*dx_3d*dy, (/1,2/))/volume
    copy_VarCoords(Tbox, Q_net_tr)

    ;;;
    
    ; Q_net_tr_annual = month_to_annual(Q_net_tr, 0)   ;;;这里是累加
; 
     Q_net_tr_annual = wgt_annual(Q_net_tr)



  
;;;total 

       total_month = Q_net_tr
       total_month = Q_net_tr + hadv_mean + vadv_mean
       
       total_annual  = Q_net_tr_annual + hadv_mean_annual + vadv_mean_annual

       resi_annual = dT_wgt - Q_net_tr_annual - hadv_mean_annual - vadv_mean_annual
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

  

   ; rc_v_lateral_south_west = regCoef_n(ispan(1, yearnum, 1),v_lateral_south_west_annual, 0, 0)
   ; rc_v_lateral_south_east = regCoef_n(ispan(1, yearnum, 1),v_lateral_south_east_annual, 0, 0)


   ; rc_u_lateral_west_south = regCoef_n(ispan(1, yearnum, 1),u_lateral_west_south_annual, 0, 0)
   ; rc_u_lateral_west_north = regCoef_n(ispan(1, yearnum, 1),u_lateral_west_north_annual, 0, 0)


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

    ; rc_v_lateral_south_west = rc_v_lateral_south_west*yearnum
    ; rc_v_lateral_south_east = rc_v_lateral_south_east*yearnum

    ; rc_u_lateral_west_south = rc_u_lateral_west_south*yearnum
    ; rc_u_lateral_west_north = rc_u_lateral_west_north*yearnum
   

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

   ; print("rc_v_lateral_south_west="+rc_v_lateral_south_west)
   ; print("rc_v_lateral_south_east="+rc_v_lateral_south_east)

   ; print("rc_u_lateral_west_south="+rc_u_lateral_west_south)
   ; print("rc_u_lateral_west_north="+rc_u_lateral_west_north)
   
   rc_Tbox = regCoef_n(ispan(1, yearnum, 1),Tbox_annual, 0, 0)

   print("rc_Tbox="+rc_Tbox)



;;;

;;;;两段平均 80-96 97-13
time = ispan(1980, 2015, 1)*1.0
m = new(36,"float")

Q_net_tr_annual_1 = time(0:17)
Q_net_tr_annual_1= dim_avg_n_Wrap(Q_net_tr_annual(0:17), 0)
Q_net_tr_annual_2  =time(18:)
Q_net_tr_annual_2= dim_avg_n_Wrap(Q_net_tr_annual(18:), 0)

total_ocean_annual_1 = time(0:17)
total_ocean_annual_1= dim_avg_n_Wrap(total_ocean_annual(0:17), 0)
total_ocean_annual_2 = time(18:)
total_ocean_annual_2= dim_avg_n_Wrap(total_ocean_annual(18:), 0)

total_annual_1 = time(0:17)
total_annual_1=  dim_avg_n_Wrap(total_annual(0:17), 0)
total_annual_2 = time(18:)
total_annual_2=  dim_avg_n_Wrap(total_annual(18:), 0)

dT_wgt_mean_1 = time(0:17)
dT_wgt_mean_1  =dim_avg_n_Wrap(dT_wgt(0:17), 0)
dT_wgt_mean_2 = time(18:)
dT_wgt_mean_2  =dim_avg_n_Wrap(dT_wgt(18:), 0)

print("dT_wgt_mean_1="+dim_avg_n_Wrap(dT_wgt(0:17), 0))
print("dT_wgt_mean_2="+dim_avg_n_Wrap(dT_wgt(18:), 0))

Q_net_diff  = dim_avg_n_Wrap(Q_net_tr_annual(18:), 0) - dim_avg_n_Wrap(Q_net_tr_annual(0:17), 0)
ocean_diff = dim_avg_n_Wrap(total_ocean_annual(18:), 0) - dim_avg_n_Wrap(total_ocean_annual(0:17), 0)
dT_wgt_diff = dim_avg_n_Wrap(dT_wgt(18:), 0) - dim_avg_n_Wrap(dT_wgt(0:17), 0)
ocean_vadv_diff = dim_avg_n_Wrap(vadv_mean_annual(18:), 0) - dim_avg_n_Wrap(vadv_mean_annual(0:17), 0)
; ocean_vadv_diff = dim_avg_n_Wrap(vadv_mean(12*18:), 0) - dim_avg_n_Wrap(vadv_mean(0:12*18-1), 0)
; ocean_hadv_diff = dim_avg_n_Wrap(hadv_mean(12*18:), 0) - dim_avg_n_Wrap(hadv_mean(0:12*18-1), 0)
ocean_hadv_diff = dim_avg_n_Wrap(hadv_mean_annual(18:), 0) - dim_avg_n_Wrap(hadv_mean_annual(0:17), 0)


print("Q_net_diff="+Q_net_diff)
print("ocean_diff="+ocean_diff)
print("dT_wgt_diff="+dT_wgt_diff)
print("ocean_vadv_diff="+ocean_vadv_diff)
print("ocean_hadv_diff="+ocean_hadv_diff)


exit()


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
  wks   = gsn_open_wks("pdf", "all-heat-flux-contribute2BOXSST-wgt-all-component-mlpp-0.8-ERAI-2015-phipian")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 2015
  resxy@trXMinF = 1979
  resxy@tmXBValues = ispan(1980,2015,5)
  resxy@tmXBLabels = ispan(1980,2015,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(1979,2015,1)
  
  resxy@tmYLMode  =  "Explicit"
  resxy@tmYLValues = ispan(-3, 3, 1)
  resxy@tmYLLabels = ispan(-3,3,1)
  resxy@trYMinF = -3
  resxy@trYMaxF = 3

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6


  ; resxy@xyLineThicknesses =  (/ 5,5,5,5,3,3,3/)          ; make second line thicker
  ; resxy@xyLineColors       =  (/"black","blue","red","gray","green","pink"/)          ; change line color
  ; resxy@xyMonoDashPattern = False
  ; resxy@xyDashPatterns =(/0,0,0,0,0,0,0/)
  ; resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines","lines"/)


  resxy@xyLineThicknesses =  (/ 5,5,5,5,3,3/)          ; make second line thicker
  resxy@xyLineColors       =  (/"black","blue","red","green","pink"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0,0,0,0,0/)
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines"/)
  


  data_all = new((/5,36/),"float")
  data_all(0,:)=dT_wgt
  data_all(1,:)=Q_net_tr_annual
  data_all(2,:)=total_ocean_annual
  data_all(3,:)=hadv_mean_annual
  data_all(4,:)=vadv_mean_annual
  

  ; data_all(6,:)=total_trend
  ; data_all(7,:)=atm_trend
  ; data_all(8,:)=ocean_trend


  plot = gsn_csm_xy(wks,time,data_all, resxy)


  ;;添加直线

  res= True
  res@gsLineThicknessF = 3.0            ; line thickness   
  res@gsLineColor      = "blue"         ; color of lines
  res@gsLineDashPattern=15

  plot2 = gsn_add_polyline(wks, plot,time(0:17), Q_net_tr_annual_1, res)
  res@gsLineColor      = "red"         ; color of lines
  plot3 = gsn_add_polyline(wks, plot,time(0:17), total_ocean_annual_1, res)
  res@gsLineColor      = "blue"         ; color of lines
  plot5 = gsn_add_polyline(wks, plot,time(18:), Q_net_tr_annual_2, res)
  res@gsLineColor      = "red"         ; color of lines
  plot6 = gsn_add_polyline(wks, plot,time(18:), total_ocean_annual_2, res)
  
  ; res@gsLineColor      = "black"         ; color of lines
  
  ;  plot4 = gsn_add_polyline(wks, plot,time(0:17), dT_wgt_mean_1, res)
  ;  plot7 = gsn_add_polyline(wks, plot,time(18:), dT_wgt_mean_2, res)

  

  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       = resxy@xyLineColors
  ; lineres@lgLineThicknesses = (/ 3,3,3,3,3,3/)
   lineres@lgLineThicknesses = (/ 3,3,3,3,3/)
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = resxy@xyDashPatterns


  textres = True
  ; textres@lgLabels = (/"dT","Qnet","Ocean advection","Residual","Horizontal advection","Vertical advection"/)     ;加线的标注
  textres@lgLabels = (/"dT","Qnet","Ocean advection","Horizontal advection","Vertical advection"/)     ;加线的标注
  textres@lgLineLabelFontHeights = 0.005


  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)
  


  end

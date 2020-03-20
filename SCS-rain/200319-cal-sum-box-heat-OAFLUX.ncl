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
    dim = dimsizes(total_month)
    total_annual_test = new(dim(0)/12,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,dim(0)/12-2,1
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
    
    oayear = ispan(1984,2007,1)
    fils3 = "/home/yangsong3/data-observation/OAflux/monthly/netheat/qnet_"+oayear+".nc"
    print(fils3)
    f3 := addfiles(fils3, "r")

    date:=cd_calendar(f3[:]->time, 1)

    timeind:=ind(date.le.200712.and.date.ge.198401)
  
   

    qnet = short2flt(f3[:]->qnet(timeind,{lat1:lat2},{lon1:lon2}))
  
    qnet@_FillValue = 1e+20
    qnet@missing_value = 1e+20
    printVarSummary(qnet)
    qnet = where(qnet.eq.3276.6, qnet@_FillValue, qnet)

fils4 = "/home/yangsong3/data-observation/OAflux/monthly/radiation/sw_isccp_"+oayear+".nc"
f4 := addfiles(fils4, "r")

nswrs  = short2flt(f4[:]->nswrs(:,{lat1:lat2},{lon1:lon2}))
nswrs@_FillValue = 1e+20
nswrs@missing_value = 1e+20
printVarSummary(nswrs)
nswrs = where(nswrs.eq.3276.6, nswrs@_FillValue, nswrs)

fils5 = "/home/yangsong3/data-observation/OAflux/monthly/radiation/lw_isccp_"+oayear+".nc"
f5 := addfiles(fils5, "r")

nlwrs  = short2flt(f5[:]->nlwrs(:,{lat1:lat2},{lon1:lon2}))
nlwrs@_FillValue = 1e+20
nlwrs@missing_value = 1e+20
printVarSummary(nlwrs)
nlwrs = where(nlwrs.eq.3276.6, nlwrs@_FillValue, nlwrs)


fils6 = "/home/yangsong3/data-observation/OAflux/monthly/turbulence/lh_oaflux_"+oayear+".nc"
f6 := addfiles(fils6, "r")

lhtfl  = short2flt(f6[:]->lhtfl(:,{lat1:lat2},{lon1:lon2}))
lhtfl@_FillValue = 1e+20
lhtfl@missing_value = 1e+20
printVarSummary(lhtfl)
lhtfl = where(lhtfl.eq.3276.6, lhtfl@_FillValue, lhtfl)


fils7 = "/home/yangsong3/data-observation/OAflux/monthly/turbulence/sh_oaflux_"+oayear+".nc"
f7:= addfiles(fils7, "r")

shtfl  = short2flt(f7[:]->shtfl(:,{lat1:lat2},{lon1:lon2}))
shtfl@_FillValue = 1e+20
shtfl@missing_value = 1e+20
printVarSummary(shtfl)
shtfl = where(shtfl.eq.3276.6, shtfl@_FillValue, shtfl)



    ;;;convert to units degree ℃/year
      ; density = 1025  ;;kg/m3
      ; cp  = 3850         ;;j/kg℃
      ; coef  = density*cp
      coef = 4.088*10^6
         ;;;nswrs_pen 
     R = 0.58
     L1 = 0.35
     L2 = 23
     nswrs_pen = nswrs
     ; nswrs_pen =  nswrs*(R*exp(-1*mlt/L1)+ (1-R)*exp(-1*mlt/L2))
     mo = (R*exp(-1*45.6/L1)+ (1-R)*exp(-1*45.6/L2))
     nswrs_pen = mo*nswrs
      

    
      qnet = 30*86400*qnet/coef
      
      lhtfl = 30*86400*lhtfl/coef
      shtfl = 30*86400*shtfl/coef
      nswrs = 30*86400*nswrs/coef
      nlwrs = 30*86400*nlwrs/coef
      nswrs_pen = 30*86400*nswrs_pen/coef
    ;;;
    Q_net = lhtfl
    Q_net = nswrs - nlwrs - lhtfl -shtfl - nswrs_pen
    Q_net = qnet - nswrs_pen
     printVarSummary(Q_net)
    
    ;;;


    ;;; Heat transport 

    printVarSummary(Q_net)
    Q_net&lat@units = "degrees_north"
    Q_net&lon@units = "degrees_east"


    
    Q_net(:,{x1:x2},{y1:y2}) = where(ismissing(sst(12*4:12*28-1,{0:15},{110:155})),Q_net@missing_value, Q_net(:,{x1:x2},{y1:y2}))  ;;;去掉陆地上的辐射 
    dim:= dimsizes(Q_net(:,{x1:x2},{y1:y2}))
    dx_3d :=conform_dims(dim ,dx({0:15}) , 1)
    
    Q_net_tr = dim_sum_n_Wrap(Q_net(:,{x1:x2},{y1:y2})*dx_3d*dy, (/1,2/))/volume(0)
    copy_VarCoords(Tbox, Q_net_tr)

    ;;;
    
    ; Q_net_tr_annual = (/month_to_annual(Q_net_tr, 0)/)   ;;;这里是累加
; 
     Q_net_tr_annual = wgt_annual(Q_net_tr)
     printVarSummary(Q_net_tr_annual)
     
     ;;;


  Q_net_tr_annual_new  = new(34,float)
  Q_net_tr_annual_new(4:27) = Q_net_tr_annual
 
;;
 total_annual = total_ocean_annual + Q_net_tr_annual_new


     
    

    


print("*********开始画图:时间序列********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "OAFLUX-all-heat-flux-contribute2BOXSST-wgt-all-component")
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
  data_all(2,:)=Q_net_tr_annual_new
  data_all(3,:)=total_ocean_annual
  data_all(4,:)=hadv_mean_annual
  data_all(5,:)=vadv_mean_annual
  

  ; data_all(6,:)=total_trend
  ; data_all(7,:)=atm_trend
  ; data_all(8,:)=ocean_trend
time = ispan(1980, 2013, 1)

  plot = gsn_csm_xy(wks,time,data_all, resxy)








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

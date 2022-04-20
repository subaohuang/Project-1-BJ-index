;;********************************************************
;;计算SST倾向以及SST值与净长波辐射，短波辐射，感热，潜热以及平流和对流,
;;其中热通量的方向都是向上为正，向下为负。。
;;;的回归场，热通量数据是NCEP-NCAR数据，海洋数据是SODA2.2.4,SST是哈哈得来SST数据
;;由于数据之间的精度不同，需要进行插值预处理,SODA精度是0.5×0.5,HADISST是1×1
;;NCEP-NCAR的精度为1.87×1.9
;;s
;;                                         lnshuheng 2020.02.18
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

    
    ;;;sst 倾向值
     dsst = center_finite_diff_n(sst, 1, False, 0, 0) 
     copy_VarCoords(sst, dsst)
     dsst_month = clmMonTLL(dsst(:,{lat1:lat2},{lon1:lon2}))
     dsst_month_scs = dim_avg_n_Wrap(dsst_month(:,{0:15},{110:155}), (/1,2/))
     printVarSummary(dsst_month_scs)
        
     
    ;;;
    delete(date)
    delete(timeind)
;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_mlt_mn_1980-2015_chazhi.nc", "r")   


    indyear = (2013-1980+1)*12-1
    u = f2->u(0:indyear,:,{lat1:lat2},{lon1:lon2})
    v = f3->v(0:indyear,:,{lat1:lat2},{lon1:lon2})
    wt = f4->wt(0:indyear,:,{lat1:lat2},{lon1:lon2})
    T = f5->temp(0:indyear,:,{lat1:lat2},{lon1:lon2})
    mlt = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})

 
    ;;;calculate Horizontal and Vertical advction

    gradLatLon  = grad_latlon_cfd (T, T&lat, T&lon, True, False)
                             
    T_grad_lat  = gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list' 
    T_grad_lon  = gradLatLon[1]
    delete(gradLatLon)
    
    hadv = u
    hadv = u*T_grad_lon + v*T_grad_lat
    hadv = (/-1*hadv/)

    
    T_grad_z = center_finite_diff_n(T, T&depth,False,0,1)
    
    vadv = wt
    vadv = wt*T_grad_z

    ;;;convert to units degree ℃/month
      hadv =31*86400*hadv
      vadv =31*86400*vadv
    ;;;
    
    ;;;;垂直积分平均，在混合层内进行
    mlt_clm = dim_avg_n_Wrap(mlt, 0)

    z =  vadv&depth
    z&depth = z
    
    dz = dpres_plevel(z, dim_avg_n_Wrap(mlt, 0), 0, 0)
   
    dz_4d = conform_dims(dimsizes(u), dz, (/1,2,3/))

    hadv_mean = dim_sum_n_Wrap(hadv*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    vadv_mean  = dim_sum_n_Wrap(vadv*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
   
    copy_VarCoords(u(:,1,:,:), hadv_mean)
    copy_VarCoords(u(:,1,:,:), vadv_mean)

    

    ;;;季节气候态
    
       hadv_mean_month = clmMonTLL(hadv_mean)
       vadv_mean_month = clmMonTLL(vadv_mean)

       hadv_month_scs = dim_avg_n_Wrap(hadv_mean_month(:,{0:15},{110:155}), (/1,2/))
       vadv_month_scs = dim_avg_n_Wrap(vadv_mean_month(:,{0:15},{110:155}), (/1,2/))

       total_ocean_month_scs = hadv_month_scs
       total_ocean_month_scs = hadv_month_scs + vadv_month_scs

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
     nswrs_pen =  nswrs*(R*exp(-1*mlt/L1)+ (1-R)*exp(-1*mlt/L2))
     ; nswrs_pen = 0.05*nswrs
    ;;;

    ;;;convert to units degree ℃/year
      density = 1025  ;;kg/m3
      cp  = 3850         ;;j/kg℃
      coef  = mlt*density*cp

      lhtfl = 31*86400*lhtfl/coef
      shtfl = 31*86400*shtfl/coef
      nswrs = 31*86400*nswrs/coef
      nlwrs = 31*86400*nlwrs/coef
      nswrs_pen = 31*86400*nswrs_pen/coef
    ;;;

    
    ;;;转化成方向向下

    lhtfl = (/-1*lhtfl/)
    shtfl = (/-1*shtfl/)
    nswrs = (/-1*nswrs/)
    nlwrs = (/-1*nlwrs/)
    nswrs_pen = (/-1*nswrs_pen/)


    lhtfl_month =clmMonTLL(lhtfl)
    shtfl_month = clmMonTLL(shtfl)
    nswrs_month = clmMonTLL(nswrs)
    nlwrs_month = clmMonTLL(nlwrs)
    nswrs_pen_month = clmMonTLL(nswrs_pen)
    
    Q_net_month = lhtfl_month
    Q_net_month  = lhtfl_month + shtfl_month + nswrs_month + nlwrs_month -nswrs_pen_month
     
    Q_net_month_scs = dim_avg_n_Wrap(Q_net_month(:,{0:15},{110:155}), (/1,2/)) 
;;;

;;;total 
total_month_scs = Q_net_month_scs
total_month_scs = Q_net_month_scs + total_ocean_month_scs
ocean_eddy = total_month_scs
ocean_eddy = dsst_month_scs -  total_month_scs 

;;;



print("*********开始画图:时间序列********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "all-heat-flux-season-cycle")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 11
  resxy@trXMinF = 0
  resxy@tmXBValues = ispan(0,11,1)
  resxy@tmXBLabels = ispan(0,11,1)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(0,11,1)
  

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.3
  resxy@vpWidthF = 0.7


  resxy@xyLineThicknesses =  (/ 5,3,2,2,2,2/)          ; make second line thicker
  
  resxy@xyLineColors       =  (/"black","gray","blue","red","orange","green"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0,0,0,0,0/)
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines","lines"/)
  


  data_all = new((/6,12/),"float")
  data_all(0,:) = dsst_month_scs
  data_all(1,:) = total_month_scs
  data_all(2,:) = Q_net_month_scs
  data_all(3,:) = hadv_month_scs
  data_all(4,:) = vadv_month_scs
  data_all(5,:) = ocean_eddy
  

  time = ispan(0, 11, 1)  
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
  textres@lgLabels = (/"dT","total","Q_net","Horizontal advection","Vertical advection","ocean_eddy"/)     ;加线的标注

  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)
  


  end





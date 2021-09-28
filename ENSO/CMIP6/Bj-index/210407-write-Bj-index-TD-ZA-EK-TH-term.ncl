;;;;;
;;****************************************************************
;; calculate BJ index terms 
;; CMIP6 variables used are uo,vo,wo,tauu,hdfs,thetao 
;;;;*************************************************************

begin

;读入数据

   ;  model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-CSM2-MR",\
   ;            "CAMS-CSM1-0","CanESM5", "CAS-ESM2-0" ,"CESM2","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","CNRM-CM6-1"\
   ; , "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
   ;   "GISS-E2-1-G","HadGEM3-GC31-LL","IPSL-CM6A-LR","MIROC-ES2L","MIROC6", \
   ;  "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3",  "NorESM2-LM"/)
    model_name = (/"CanESM5"/)
    dim_model = dimsizes(model_name)

    plot = new(dim_model,graphic)
    iplot= plot
    plot1 = plot

    yearstart = 190001
    yearend = 201412

    term_BJ = new((/dim_model,5/),float)
    

    startyear_r = str_split_by_length(yearstart, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(yearend, 4 )

    yearnum = stringtoint(endyear_r(0)) - stringtoint(startyear_r(0)) + 1
;;;
    lat1 = -10
    lat2 = 10
    lon1 = 120
    lon2 = 290

do k = 0,dim_model-1
;;;;

   model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils1 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon_levint/thetao_Omon_"+model+"_historical_"+"*")

    f1 := addfiles(fils1, "r")

    time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)   

    tos = lonFlip(f1[:]->thetao(timeind,{0:50},{lat1:lat2},:))
    
    

    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos_ar = rmMonAnnCycLLLT(tos(lev|:,lat|:,lon|:,time|:))
    tos_a = tos_ar(time|:,lev|:,lat|:,lon|:)
    tos_clm = dim_avg_n_Wrap(tos, 0)  ;;气候态的海温
;;;

;;parameter 
     pi = get_pi("float")
     dx = tos&lat
     dx = 6378388.*cos(0.0174533*tos&lat)
     dx!0 = "lat"
     copy_VarCoords(tos(1,1,:,1), dx)
     dy = tos&lon
     dy = (pi/180)*1*6378388.
     dy!0 = "lon"
     copy_VarCoords(tos(1,1,1,:), dy)
     depth = tofloat(tos&lev)
     dz = depth
     dz(1:)= depth(1:) - depth(0:dimsizes(depth)-2)
     dz(0)= depth(0) 
;;;


;;;;计算T_clm的纬向梯度
  
     dtdxdy = grad_latlon_cfd(tos_clm, tos_clm&lat, tos_clm&lon, False, False)
     dtdx = dtdxdy[1]
     copy_VarCoords(tos_clm, dtdx)
;;;;

print("*********************ZA feedback********************************")

;;;计算东box区域dtdx平均，混合层深度取50 ，185.2E - 80w -5-5N
     dtdx_vol = wgt_volave(dtdx(:,{-5:5},{185:280}), dz, 1, 1, 0)
;;;
   
;;;计算东box区域异常海温平均，混合层深度取50 ，185.2E - 80w -5-5N
     Tbox_a = wgt_volave(tos_a(:,:,{-5:5},{185:280}), dz, 1, 1, 0)
;;;


;;;;读入洋流数据
    fils2 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon_levint/uo_Omon_"+model+"_historical_"+"*")
    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    uo = lonFlip(f2[:]->uo(timeind,{0:50},{lat1:lat2},:))

    uo = (/dtrend_msg_n(ispan(1,dimsizes(uo&time),1),uo,False,True,0)/)
    uo_ar = rmMonAnnCycLLLT(uo(lev|:,lat|:,lon|:,time|:))
    uo_a = uo_ar(time|:,lev|:,lat|:,lon|:)
    uo_clm = dim_avg_n_Wrap(uo, 0)
;;;;

;;;计算东box区域异常纬向洋流平均，混合层深度取50 ，185.2E - 80w -5-5N
     ubox_a = wgt_volave(uo_a(:,:,{-5:5},{185:280}), dz, 1, 1, 0)
;;;


;;;;读入纬向风应力数据
    fils3 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/tauu_Amon_"+model+"_historical_"+"*")
    f3 := addfiles(fils3, "r")

    time:=f3[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    tauu = lonFlip(f3[:]->tauu(timeind,{lat1:lat2},:))

    tauu = (/dtrend_msg_n(ispan(1,dimsizes(tauu&time),1),tauu,False,True,0)/)
    tauu_a = rmMonAnnCycTLL(tauu) 
    tauu_clm = dim_avg_n_Wrap(tauu, 0)
;;;;

 
;;;;miua  taux = miua* Tbox_a
    tauu_pac = dim_avg_n_Wrap(tauu_a(:,{-5:5},{120:280}),(/1,2/))
    miua  =  regCoef_n(Tbox_a, tauu_pac, 0, 0) ;; 风对海温异常的响应系数
    print("miua="+miua)
;;;;

;;;betau ubox_a = betau*tauu_pac

    betau  =  regCoef_n(tauu_pac, ubox_a, 0, 0) ;; 区域洋流对热带纬向风的响应
    print("betau="+betau)
;;;

;;;****ZA feedback miua*betau* -1*dtdx_vol
   ZA = miua*betau*-1*dtdx_vol
   print("****ZA value="+ZA)
;;;


print("*********************EK feedback********************************")


;;;;计算T_clm的垂直梯度
  
   dtdz = center_finite_diff_n(tos_clm, tos_clm&lev, False, 1,0) 
   copy_VarCoords(tos_clm, dtdz)
   
;;;计算东box区域dtdx平均，混合层深度取50 ，185.2E - 80w -5-5N
    dtdz_vol = wgt_volave(dtdz(:,{-5:5},{185:280}), dz, 1, 1, 0)
;;;
    


;;;;读入垂直洋流数据
    fils4 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon_levint/wo_Omon_"+model+"_historical_"+"*")
    f4 := addfiles(fils4, "r")

    time:=f4[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    wo = lonFlip(f4[:]->wo(timeind,{0:50},{lat1:lat2},:))

    wo = (/dtrend_msg_n(ispan(1,dimsizes(wo&time),1),wo,False,True,0)/)
    wo_ar = rmMonAnnCycLLLT(wo(lev|:,lat|:,lon|:,time|:))
    wo_a = wo_ar(time|:,lev|:,lat|:,lon|:)
    wo_clm = dim_avg_n_Wrap(wo(:,{50},:,:), 0)  ;;;注意,BJ 公式中W的气候态要取混合层底的速度
    wo_clm_4d =  conform_dims(dimsizes(wo_a),wo_clm, (/2,3/))
                              
; ;;;;
;      depth_w = tofloat(wo&lev)
;      dz_w = depth_w
;      if depth_w(0).eq.0.0 then 
;         depth_w(0) = (depth_w(1) + depth_w(0))/2.0
;      end if 
  
;      dz_w(1:)= depth_w(1:) - depth_w(0:dimsizes(depth_w)-2)
;      dz_w(0)= depth_w(0) 
; ;;;

;;;betaw : wbox_a = -beatw*tauu_pac
    wo_a = where(wo_clm_4d.ge.0, wo_a, 0.0)

   ;;;计算东box区域异常垂直洋流平均，混合层深度取50 ，185.2E - 80w -5-5N
    wbox_a = wgt_volave(wo_a(:,:,{-5:5},{185:280}), dz, 1, 1, 0)
   ;;;
    betaw  =  regCoef_n(-1*tauu_pac, wbox_a, 0, 0) ;; 区域ekman流对热带纬向风的响应
    print("betaw="+betaw)
;;;
   

;;;***EK feedback miua*betaw* -1*dtdz_vol
   EK = miua*betaw*-1*dtdz_vol
   print("****EK value="+EK)
;;;

print("*********************TH feedback********************************")


;;;get ocean heat content anomaly
depth_total = tofloat(f1[0]->lev({0:300}))
 
dz_total = depth_total
dz_total(1:)= depth_total(1:) - depth_total(0:dimsizes(depth_total)-2)
dz_total(0)= depth_total(0) 

hc = dim_avg_wgt_n_Wrap(lonFlip(f1[:]->thetao(timeind,{0:300},{lat1:lat2},:)),dz_total,1,1)
    
hc = (/dtrend_msg_n(ispan(1,dimsizes(hc&time),1),hc,False,True,0)/)
hc_a = rmMonAnnCycTLL(hc)
hc_clm = dim_avg_n_Wrap(hc, 0)  ;;气候态
;;;


;;;计算东西box区域异常hc平均，混合层深度取50 ，185.2E - 80w -5-5N
    hca_e = dim_avg_n_Wrap(hc_a(:,{-5:5},{185:280}),(/1,2/))
    hca_w = dim_avg_n_Wrap(hc_a(:,{-5:5},{120:185}),(/1,2/))

    hac_diff  = hca_e
    hac_diff = hca_e - hca_w
;;;


;;;betah : hcae - hcaw = beath*tauu_pac
    betah  =  regCoef_n(tauu_pac, hac_diff, 0, 0) ;; 区域ekman流对热带纬向风的响应
    print("betah="+betah)
;;;


    ;;;get subsurface sea warter anomaly 50-100m
    depth_sub = tofloat(f1[0]->lev({50:100}))

    dz_sub = depth_sub
    dz_sub(1:)= depth_sub(1:) - depth_sub(0:dimsizes(depth_sub)-2)
    dz_sub(0)= depth_sub(0) 

    tsub = dim_avg_wgt_n_Wrap(lonFlip(f1[:]->thetao(timeind,{50:100},{lat1:lat2},:)),dz_sub,1,1)

    tsub = (/dtrend_msg_n(ispan(1,dimsizes(hc&time),1),hc,False,True,0)/)
    tsub_a = rmMonAnnCycTLL(tsub)
    tsub_clm = dim_avg_n_Wrap(tsub, 0)  ;;气候态
    ;;;

;;;
;;;ah : tsub_a = ah*hca_e
  wo_clm_3d =  conform_dims(dimsizes(tsub_a),wo_clm, (/1,2/))
  printVarSummary(wo_clm_3d)
  printVarSummary(tsub_a)

    ; wo_clm_4d_tsub = conform_dims(dimsizes(tsub_a), wo_clm, ndim)
    tsub_a = where(wo_clm_3d.gt.0, tsub_a, 0.0)

   ;;;计算东box区域异常垂直洋流平均，混合层深度取50 ，185.2E - 80w -5-5N
     tsub_a_box = dim_avg_n_Wrap(tsub_a(:,{-5:5},{185:280}),(/1,2/))
   ;;;
    ah  =  regCoef_n(hca_e, tsub_a_box, 0, 0) ;; 区域ekman流对热带纬向风的响应
    print("ah="+ah)
;;;


 ;;;计算东box区域异常垂直洋流平均，混合层深度取50 ，185.2E - 80w -5-5N
   wbox_clm = avg(wo_clm({-5:5},{185:280}))
   wbox_clm = wbox_clm/50.0
 ;;;


;;***TH term  thermalcline feedback term TH= miua*betah*ah*wbox_clm

  TH = miua*betah*ah*wbox_clm
  print("TH="+TH)
 print("****TH value="+TH)

print("*********************Thermal Damping********************************")
  
;;;;读入Q数据
    fils5 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon/hfds_Omon_"+model+"_historical_"+"*")
    f5 := addfiles(fils5, "r")

    time:=f5[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    hfds = lonFlip(f5[:]->hfds(timeind,{lat1:lat2},:))


    Cp = 3940;1996
    rho0 = 1022.4;1026
    Hm=50.0

  ; density = 1025  ;;kg/m3
  ; cp  = 3850         ;;j/kg℃
  ; coef  = density*cp
    hfds = (/hfds/(Cp*rho0*Hm)/) ;; k/s

    hfds = (/dtrend_msg_n(ispan(1,dimsizes(hfds&time),1),hfds,False,True,0)/)
    hfds_a = rmMonAnnCycTLL(hfds)
    hfds_clm = dim_avg_n_Wrap(hfds, 0)


;;;;

;;;hfds_a_box 
    hfds_a_box = dim_avg_n_Wrap(hfds_a(:,{-5:5},{185:280}), (/1,2/))

;;;alfas  hfds_a_box = -alfas*Tbox_a
    alfas  =  regCoef_n(-1*Tbox_a, hfds_a_box, 0, 0) ;; 区域洋流对热带纬向风的响应
    print("alfas="+alfas)
;;;

;;**TD = -alfas
    TD = -alfas
    print("****TD value="+TD)
;;;;

print("***********start MA term ************************")

;;;;边界平均的海温异常 西：185,东：280E 南：-5 北 5

    Ta_eb = wgt_areaave_Wrap(tos_a(:,:,{-5:5},{280}), dz, 1, 0) 
    printVarSummary(Ta_eb)
    a1  =  regCoef_n(Tbox_a, Ta_eb, 0, 0)
    print("a1="+a1)

    ; Ta_wb = dim_sum_n_Wrap(tos_a(:,:,{-5:5},{185})*dy*dx_3d(:,:,:,1)*dz_3d(:,:,:,1), (/1,2/))/volume_2d
    Ta_wb = wgt_areaave_Wrap(tos_a(:,:,{-5:5},{185}), dz, 1, 0) 
    a2  =  regCoef_n(Tbox_a, Ta_wb, 0, 0)
    print("a2="+a2)

    Ta_nb = wgt_areaave_Wrap(tos_a(:,:,{5},{185:280}), dz, 1, 0) 
    a3  =  regCoef_n(Tbox_a, Ta_nb, 0, 0)
    print("a3="+a3)

    Ta_sb = wgt_areaave_Wrap(tos_a(:,:,{-5},{185:280}), dz, 1, 0) 
    a4  =  regCoef_n(Tbox_a, Ta_sb, 0, 0)
    print("a4="+a4)  
;;;;

    Lx = (280 - 185)*dx({0})
    Ly = 10*dy(0)
;;;

;;;;
    fils6 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon_levint/vo_Omon_"+model+"_historical_"+"*")
    f6 := addfiles(fils6, "r")

    time:=f6[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    vo = lonFlip(f6[:]->vo(timeind,{0:60},{lat1:lat2},:))

    vo = (/dtrend_msg_n(ispan(1,dimsizes(vo&time),1),vo,False,True,0)/)
    vo_ar = rmMonAnnCycLLLT(vo(lev|:,lat|:,lon|:,time|:))
    vo_a = vo_ar(time|:,lev|:,lat|:,lon|:)
    vo_clm = dim_avg_n_Wrap(vo, 0)
;;;;

 
;;;;边界平均的洋流异常 西：185,东：280E 南：-5 北 5

    uclm_eb = avg(dim_avg_wgt_n_Wrap(uo_clm(:,{-5:5},{280}), dz, 1, 0))

    uclm_wb = avg(dim_avg_wgt_n_Wrap(uo_clm(:,{-5:5},{185}), dz, 1, 0))
  
    vclm_nb =  avg(dim_avg_wgt_n_Wrap(vo_clm(:,{5},{185:280}), dz, 1, 0))
    vclm_sb =  avg(dim_avg_wgt_n_Wrap(vo_clm(:,{-5},{185:280}), dz, 1, 0))


    
;;;term MA  -((a1*ua_eb - a2*ua_wb)/Lx + (b1*va_nb - b2*va_sb)/Ly)
    MA =  -1*((a1*uclm_eb  - a2*uclm_wb)/Lx + (a3*vclm_nb  - a4*vclm_sb)/Ly)
    
    print("****MA value="+MA)
;;;
  MA = 86400*365*MA
  TD = 86400*365*TD
  ZA = 86400*365*ZA
  EK = 86400*365*EK
  TH = 86400*365*TH

  term_BJ(k,0) = tofloat(MA)
  term_BJ(k,1) = tofloat(TD)
  term_BJ(k,2) = tofloat(ZA)
  term_BJ(k,3) = tofloat(EK)
  term_BJ(k,4) = tofloat(TH)

  print(MA)
  print(TD)
  print(ZA)
  print(EK)
  print(TH)

 print("********************************finish calculate"+" "+model+"************************************")
end do

print("************"+"写入nc文件"+"*************")
term_BJ!0  = "model_name"
term_BJ!1  = "terms"

system("/bin/rm -f term_BJ_cmip6.nc")
out = addfile("term_BJ_cmip6.nc", "c")
out->term_BJ = term_BJ

   
end 



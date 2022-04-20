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

    term_MA = new(dim_model,float)
    
    startyear_r = str_split_by_length(yearstart, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(yearend, 4 )

    yearnum = stringtoint(endyear_r(0)) - stringtoint(startyear_r(0)) + 1
    print(yearnum)
;;;


 

    lat1 = -10
    lat2 = 10
    lon1 = 120
    lon2 = 290


do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils1 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon/thetao_Omon_"+model+"_historical_"+"*")

    f1 := addfiles(fils1, "r")

    time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    tos = lonFlip(f1[:]->thetao(timeind,{0:60},{lat1:lat2},:))
    printVarSummary(tos)
    exit()
    
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos_ar = rmMonAnnCycLLLT(tos(lev|:,lat|:,lon|:,time|:))
    tos_a = tos_ar(time|:,lev|:,lat|:,lon|:)
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
     ; dx_3d = conform_dims(dimsizes(tos(:,:,{-5:5},{185:280})) ,dx({-5:5}) , 2)
     ; dz_3d = conform_dims(dimsizes(tos(:,:,{-5:5},{185:280})), dz , 1)

     ; copy_VarCoords(tos(:,:,{-5:5},{185:280}), dx_3d)
     ; copy_VarCoords(tos(:,:,{-5:5},{185:280}), dz_3d)

     ; dx_3d = where(ismissing(tos(:,:,{-5:5},{185:280})), 0.0, dx_3d)
     ; dz_3d  = where(ismissing(tos(:,:,{-5:5},{185:280})), 0.0, dz_3d)

     ; volume = dim_sum_n(dx_3d*dz_3d*dy(0), (/1,2,3/))

;;;计算东box区域海温异常平均，混合层深度取50 ，185.2E - 80w -5-5N
     ; Tbox_a = wgt_volave(tos_a(:,:,{-5:5},{185:280}), dz, dx({-5:5}), dy({185:280}), 0)
     Tbox_a = wgt_volave(tos_a(:,:,{-5:5},{185:280}), dz, 1, 1, 0)
     copy_VarCoords(tos_a(:,1,1,1), Tbox_a)
;;;
    
;;;;边界平均的海温异常 西：185,东：280E 南：-5 北 5

    Ta_eb = wgt_areaave_Wrap(tos_a(:,:,{-5:5},{280}), dz, 1, 0) 
    printVarSummary(Ta_eb)
    a1  =  regCoef_n(Tbox_a, Ta_eb, 0, 0)
    print(a1)

    ; Ta_wb = dim_sum_n_Wrap(tos_a(:,:,{-5:5},{185})*dy*dx_3d(:,:,:,1)*dz_3d(:,:,:,1), (/1,2/))/volume_2d
    Ta_wb = wgt_areaave_Wrap(tos_a(:,:,{-5:5},{185}), dz, 1, 0) 
    a2  =  regCoef_n(Tbox_a, Ta_wb, 0, 0)
    print(a2)

    Ta_nb = wgt_areaave_Wrap(tos_a(:,:,{5},{185:280}), dz, 1, 0) 
    a3  =  regCoef_n(Tbox_a, Ta_nb, 0, 0)
    print(a3)

    Ta_sb = wgt_areaave_Wrap(tos_a(:,:,{-5},{185:280}), dz, 1, 0) 
    a4  =  regCoef_n(Tbox_a, Ta_sb, 0, 0)
    print(a4)  
;;;;

    Lx = (280 - 185)*dx({0})
    Ly = 10*dy
;;;

;;;;读入洋流数据
    fils2 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon/uo_Omon_"+model+"_historical_"+"*")
    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    uo = lonFlip(f2[:]->uo(timeind,{0:60},{lat1:lat2},:))

    uo = (/dtrend_msg_n(ispan(1,dimsizes(uo&time),1),uo,False,True,0)/)
    uo_ar = rmMonAnnCycLLLT(uo(lev|:,lat|:,lon|:,time|:))
    uo_a = uo_ar(time|:,lev|:,lat|:,lon|:)
    uo_clm = dim_avg_n_Wrap(uo, 0)
;;;;

;;;;
    fils3 :=  systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/Omon/vo_Omon_"+model+"_historical_"+"*")
    f3 := addfiles(fils3, "r")

    time:=f3[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    vo = lonFlip(f3[:]->vo(timeind,{0:60},{lat1:lat2},:))

    vo = (/dtrend_msg_n(ispan(1,dimsizes(vo&time),1),vo,False,True,0)/)
    vo_ar = rmMonAnnCycLLLT(vo(lev|:,lat|:,lon|:,time|:))
    vo_a = vo_ar(time|:,lev|:,lat|:,lon|:)
    vo_clm = dim_avg_n_Wrap(vo, 0)
;;;;



 
;;;;边界平均的洋流异常 西：185,东：280E 南：-5 北 5

    uclm_eb = avg(dim_avg_wgt_n_Wrap(uo_clm(:,{-5:5},{280}), dz, 1, 0))

    uclm_wb = avg(dim_avg_wgt_n_Wrap(uo_clm(:,{-5:5},{185}), dz, 1, 0))
  
    va_nb =  avg(dim_avg_wgt_n_Wrap(vo_a(:,{5},{185:280}), dz, 1, 0))
    va_sb =  avg(dim_avg_wgt_n_Wrap(vo_a(:,{-5},{185:280}), dz, 1, 0))

    

;;;term MA  -((a1*ua_eb - a2*ua_wb)/Lx + (b1*va_nb - b2*va_sb)/Ly)
   
  MA =  -1*((a1*uclm_eb  - a2*uclm_wb)/Lx + (a3*vclm_nb  - a4*vclm_sb)/Ly)

  term_MA(k) = MA

 print("********************************finish calculate"+" "+model+"************************************")
end do
   
end 



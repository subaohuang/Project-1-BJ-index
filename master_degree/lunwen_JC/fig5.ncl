function  one_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end


;;;********************************************************************
;;;主程序
;;;******************************************************************
begin



diri="/mnt/f/share/CMIP5/msf_history/"
xi_model=(/"CSIRO-Mk3-6-0","INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
; xi_model=(/"CSIRO-Mk3-6-0"/)
; xi_model=(/"INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
dong_model=(/"MPI-ESM-MR","MIROC5"/)
;;;全球
latmin=-30
latmax=30
lonmin=0
lonmax=360
lat=fspan(30, -30, 25)
lon=fspan(0, 359.5, 144)
lat_pr=fspan(30, -30, 25)
lon_pr=fspan(0, 359.5, 144)
;;;

; psl_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/),float)
; psl_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/),float)


; pr_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat_pr),dimsizes(lon_pr)/) ,float)
; pr_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat_pr),dimsizes(lon_pr)/),float)

; uwnd_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/) ,float)
; uwnd_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/), float)

; vwnd_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/), float)
; vwnd_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/),float)

psl_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/),float)
psl_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat),dimsizes(lon)/),float)


pr_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat_pr),dimsizes(lon_pr)/) ,float)
pr_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat_pr),dimsizes(lon_pr)/),float)

uwnd_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/) ,float)
uwnd_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat),dimsizes(lon)/), float)

vwnd_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/), float)
vwnd_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat),dimsizes(lon)/),float)

;;;;AMIP
psl_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/),float)
psl_amip_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat),dimsizes(lon)/),float)


pr_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat_pr),dimsizes(lon_pr)/) ,float)
pr_amip_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat_pr),dimsizes(lon_pr)/),float)

uwnd_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/) ,float)
uwnd_amip_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat),dimsizes(lon)/), float)

vwnd_amip_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat),dimsizes(lon)/), float)
vwnd_amip_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat),dimsizes(lon)/),float)
;;;


do i=0,1
    if i.eq.0 then 
     dim=dimsizes(xi_model)   
     models=xi_model
 else
    dim=dimsizes(dong_model)
    models=dong_model
   end if

 do j=0,dim-1
    
    model=models(j)
    yearstart=197901
    yearend=200812
    
    year1=1979
    year2=2008
    
   ;;;读入降水数据*************************************************************************************
    print("******************"+model+"*********************")
    fils1="/mnt/i/CMIP5/history/pr_chazhi/pr_"+model+"_historical_r1i1p1_197901-200512.nc"
    fils2= "/mnt/i/CMIP5/rcp45/pr_chazhi/pr_"+model+"_rcp45_r1i1p1_200601-201212.nc"

    printVarSummary(fils1)
    printVarSummary(fils2)
    
    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils1)
    delete(fils2)
     
    time=f1[:]->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    pr=f1[:]->pr(timeind,{30:-30},{0:360})
    printVarSummary(pr)
    ;;;换算成mm/day
    pr=(/pr*3600*24.0/)

    delete(f1)
    delete(pr&time)
    pr&time=date(timeind)

    delete(time)
    delete(date)
    printVarSummary(pr)
    delete(timeind)
    
    pr_annual=month_to_annual(pr, 1)
    pr_annual!0="time"
    pr_annual&time=ispan(year1,year2,1)
    delete(pr)
    
    ;;;AMIP的数据
    fils3= "/mnt/f/share/CMIP5/AMIP/pr_chazhi/pr_"+model+"_AMIP_r1i1p1_197901-200812.nc"
    f3=addfile(fils3, "r")
    time=f3->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    pr_amip=f3->pr(timeind,{30:-30},{0:360})
    delete(f3)
    delete(fils3)
    ;;;换算成mm/day
    pr_amip=(/pr_amip*3600*24.0/)
    pr_amip_annual=month_to_annual(pr_amip, 1)
    pr_amip_annual!0="time"
    pr_amip_annual&time=ispan(year1,year2,1)
    delete(pr_amip)
    delete(time)
    delete(date)



    ;***********************************************************************************************
    ;;
    ;;;**********************************************************************************************
    print("******************读入风场*********************")
    fils1=systemfunc("ls /mnt/i/CMIP5/history/uv/uv_"+model+"_history_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/rcp45/uv/uv_"+model+"_rcp45_r1i1p1_"+"*"+".nc")
    print(fils1)
    print(fils1)

    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils1)
    delete(fils2)
     
    time=f1[:]->time    
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd=f1[:]->uwnd(timeind,{30:-30},{0:360})  ;;;
    vwnd=f1[:]->vwnd(timeind,{30:-30},{0:360})  ;;;
    delete(f1)
   
;;;;
    delete(uwnd&time)
    delete(vwnd&time)
    uwnd&time=date(timeind)
    vwnd&time=date(timeind)
     
    uwnd_annual=month_to_annual(uwnd, 1)
    uwnd_annual!0="time"
    uwnd_annual&time=ispan(year1,year2,1)
    vwnd_annual=month_to_annual(vwnd, 1)
    vwnd_annual!0="time"
    vwnd_annual&time=ispan(year1,year2,1)
    delete(uwnd)
    delete(vwnd)
    delete(timeind)
    delete(time)
    delete(date)
    ;;;;;;;

    ;;;AMIP的风场数据
    fils3=systemfunc("ls /mnt/f/share/CMIP5/AMIP/uv/uv_"+model+"_amip_r1i1p1_"+"*"+".nc") 
    f3=addfile(fils3, "r")
    time=f3->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    uwnd_amip=f3->uwnd(timeind,{30:-30},{0:360})
    vwnd_amip=f3->vwnd(timeind,{30:-30},{0:360})
    delete(f3)
    delete(fils3)
    ;;;换算成mm/day
    uwnd_amip_annual=month_to_annual(uwnd_amip, 1)
    uwnd_amip_annual!0="time"
    uwnd_amip_annual&time=ispan(year1,year2,1)
    delete(uwnd_amip)
    vwnd_amip_annual=month_to_annual(vwnd_amip, 1)
    vwnd_amip_annual!0="time"
    vwnd_amip_annual&time=ispan(year1,year2,1)
    delete(vwnd_amip)
    delete(timeind)
    delete(time)
    delete(date)
    ;;;;

;;;**************************************************************************************
;;;;**********************************************************************************
;;;;***************************************************************************


   print("******************读入海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/i/CMIP5/history/slp_chazhi/psl_"+model+"_historical_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/RCP45/slp_chazhi/psl_"+model+"_rcp45_r1i1p1_"+"*"+".nc")
    print(fils1)
    print(fils2)


    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils1)
    delete(fils2)
   
    time=f1[:]->time
    date=cd_calendar(time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    printVarSummary(timeind)
    psl=lonPivot(f1[:]->psl(0:359,:,:),0)   
    delete(f1)
    psl=(/psl/)/100.0
   
     
    psl_annual=month_to_annual(psl(:,{30:-30},{0:360}), 1)
    psl_annual!0="time"
    psl_annual&time=ispan(year1,year2,1)
    delete(timeind)
    delete(psl)
    delete(time)
    delete(date)

      ;;;AMIP的数据
    fils3= "/mnt/f/share/CMIP5/AMIP/slp_chazhi/psl_"+model+"_amip_r1i1p1_197901-200812.nc"
    f3=addfile(fils3, "r")
    time=f3->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    psl_amip=lonPivot(f3->psl(0:359,:,:),0) 
    psl_amip=(/psl_amip/100.0/)
    delete(f3)
    delete(fils3)
    psl_amip_annual=month_to_annual(psl_amip(:,{30:-30},{0:360}), 1)
    psl_amip_annual!0="time"
    psl_amip_annual&time=ispan(year1,year2,1)
    delete(psl_amip)
    delete(time)
    delete(timeind)
    delete(date)
    ;;;
    ;;;;

;;;***********************************************************************************************
;;;
;;;;'**********************************************************************************************'
 printVarSummary(psl_annual)


    if i.eq.0 then 
      psl_ximodel(j,:,:,:)=psl_annual
      pr_ximodel(j,:,:,:)=pr_annual
      uwnd_ximodel(j,:,:,:)=uwnd_annual
      vwnd_ximodel(j,:,:,:)=vwnd_annual

      psl_amip_ximodel(j,:,:,:)=psl_amip_annual
      pr_amip_ximodel(j,:,:,:)=pr_amip_annual
      uwnd_amip_ximodel(j,:,:,:)=uwnd_amip_annual
      vwnd_amip_ximodel(j,:,:,:)=vwnd_amip_annual
   else
    psl_dongmodel(j,:,:,:)=psl_annual
    pr_dongmodel(j,:,:,:)=pr_annual
    uwnd_dongmodel(j,:,:,:)=uwnd_annual
    vwnd_dongmodel(j,:,:,:)=vwnd_annual
    psl_amip_dongmodel(j,:,:,:)=psl_amip_annual
    pr_amip_dongmodel(j,:,:,:)=pr_amip_annual
    uwnd_amip_dongmodel(j,:,:,:)=uwnd_amip_annual
    vwnd_amip_dongmodel(j,:,:,:)=vwnd_amip_annual
    end if

    delete([/psl_annual,pr_annual,uwnd_annual,vwnd_annual/])
    delete([/psl_amip_annual,pr_amip_annual,uwnd_amip_annual,vwnd_amip_annual/])
end do
delete(dim)
delete(models)
end do


psl_xi_mean=dim_avg_n(psl_ximodel,0)
pr_xi_mean=dim_avg_n(pr_ximodel,0)
uwnd_xi_mean=dim_avg_n(uwnd_ximodel,0)
vwnd_xi_mean=dim_avg_n(vwnd_ximodel,0)

psl_dong_mean=dim_avg_n(psl_dongmodel,0)
pr_dong_mean=dim_avg_n(pr_dongmodel,0)
uwnd_dong_mean=dim_avg_n(uwnd_dongmodel,0)
vwnd_dong_mean=dim_avg_n(vwnd_dongmodel,0)

psl_xi_mean!0="time"
psl_xi_mean!1="lat"
psl_xi_mean!2="lon"

psl_xi_mean&time=ispan(year1, year2, 1)
psl_xi_mean&lat=lat
psl_xi_mean&lon=lon
psl_xi_mean&lat@units="degrees_north"
psl_xi_mean&lon@units="degrees_east"
printVarSummary(psl_xi_mean)

copy_VarCoords(psl_xi_mean, psl_dong_mean)
copy_VarCoords(psl_xi_mean, pr_xi_mean)
copy_VarCoords(psl_xi_mean, pr_dong_mean)
copy_VarCoords(psl_xi_mean, uwnd_xi_mean)
copy_VarCoords(psl_xi_mean, uwnd_dong_mean)
copy_VarCoords(psl_xi_mean, vwnd_xi_mean)
copy_VarCoords(psl_xi_mean, vwnd_dong_mean)



psl_xi_mean_amip=dim_avg_n(psl_amip_ximodel,0)
pr_xi_mean_amip=dim_avg_n(pr_amip_ximodel,0)
uwnd_xi_mean_amip=dim_avg_n(uwnd_amip_ximodel,0)
vwnd_xi_mean_amip=dim_avg_n(vwnd_amip_ximodel,0)

psl_dong_mean_amip=dim_avg_n(psl_amip_dongmodel,0)
pr_dong_mean_amip=dim_avg_n(pr_amip_dongmodel,0)
uwnd_dong_mean_amip=dim_avg_n(uwnd_amip_dongmodel,0)
vwnd_dong_mean_amip=dim_avg_n(vwnd_amip_dongmodel,0)

psl_xi_mean_amip!0="time"
psl_xi_mean_amip!1="lat"
psl_xi_mean_amip!2="lon"

psl_xi_mean_amip&time=ispan(year1, year2, 1)
psl_xi_mean_amip&lat=lat
psl_xi_mean_amip&lon=lon
psl_xi_mean_amip&lat@units="degrees_north"
psl_xi_mean_amip&lon@units="degrees_east"
printVarSummary(psl_xi_mean_amip)

copy_VarCoords(psl_xi_mean_amip, psl_dong_mean_amip)
copy_VarCoords(psl_xi_mean_amip, pr_xi_mean_amip)
copy_VarCoords(psl_xi_mean_amip, pr_dong_mean_amip)
copy_VarCoords(psl_xi_mean_amip, uwnd_xi_mean_amip)
copy_VarCoords(psl_xi_mean_amip, uwnd_dong_mean_amip)
copy_VarCoords(psl_xi_mean_amip, vwnd_xi_mean_amip)
copy_VarCoords(psl_xi_mean_amip, vwnd_dong_mean_amip)



    print("*****************计算趋势************************")  

    pr_rc1_xi_mean=regCoef_n(ispan(1,dimsizes(pr_xi_mean&time),1),pr_xi_mean(:,:,:),0,0)  
    copy_VarCoords(pr_xi_mean(0,:,:), pr_rc1_xi_mean)
    delete(pr_xi_mean)
    printVarSummary(pr_rc1_xi_mean)

    pr_rc1_dong_mean=regCoef_n(ispan(1,dimsizes(pr_dong_mean&time),1),pr_dong_mean(:,:,:),0,0)  
    copy_VarCoords(pr_dong_mean(0,:,:), pr_rc1_dong_mean)
    delete(pr_dong_mean)
    printVarSummary(pr_rc1_dong_mean)
 
    uwnd_rc1_xi_mean=regCoef_n(ispan(1,dimsizes(uwnd_xi_mean&time),1),uwnd_xi_mean(:,:,:),0,0)  
    copy_VarCoords(uwnd_xi_mean(0,:,:), uwnd_rc1_xi_mean)
    delete(uwnd_xi_mean)
    printVarSummary(uwnd_rc1_xi_mean)

     
    uwnd_rc1_dong_mean=regCoef_n(ispan(1,dimsizes(uwnd_dong_mean&time),1),uwnd_dong_mean(:,:,:),0,0)  
    copy_VarCoords(uwnd_dong_mean(0,:,:), uwnd_rc1_dong_mean)
    delete(uwnd_dong_mean)
    printVarSummary(uwnd_rc1_dong_mean)
    
    vwnd_rc1_xi_mean=regCoef_n(ispan(1,dimsizes(vwnd_xi_mean&time),1),vwnd_xi_mean(:,:,:),0,0)  
    copy_VarCoords(vwnd_xi_mean(0,:,:), vwnd_rc1_xi_mean)
    delete(vwnd_xi_mean)
    printVarSummary(vwnd_rc1_xi_mean)

     
    vwnd_rc1_dong_mean=regCoef_n(ispan(1,dimsizes(vwnd_dong_mean&time),1),vwnd_dong_mean(:,:,:),0,0)  
    copy_VarCoords(vwnd_dong_mean(0,:,:), vwnd_rc1_dong_mean)
    delete(vwnd_dong_mean)
    printVarSummary(vwnd_rc1_dong_mean)


    psl_rc1_xi_mean=regCoef_n(ispan(1,dimsizes(psl_xi_mean&time),1),psl_xi_mean(:,:,:),0,0)  
    copy_VarCoords(psl_xi_mean(0,:,:), psl_rc1_xi_mean)
    printVarSummary(pr_rc1_xi_mean)
    delete(psl_xi_mean)
    psl_rc1_dong_mean=regCoef_n(ispan(1,dimsizes(psl_dong_mean&time),1),psl_dong_mean(:,:,:),0,0)  
    copy_VarCoords(psl_dong_mean(0,:,:), psl_rc1_dong_mean)
    delete(psl_dong_mean)
    printVarSummary(pr_rc1_dong_mean)

    prob1=one_test(psl_rc1_xi_mean)
    prob2=one_test(psl_rc1_dong_mean)
 

    uwnd_xi_prob1=one_test(uwnd_rc1_xi_mean)
    vwnd_xi_prob1=one_test(vwnd_rc1_xi_mean)
    
    uwnd_dong_prob1=one_test(uwnd_rc1_xi_mean)
    vwnd_dong_prob1=one_test(vwnd_rc1_dong_mean)

    uwnd_rc1_xi_mean=where(uwnd_xi_prob1.le.0.1, uwnd_rc1_xi_mean, uwnd_rc1_xi_mean@_FillValue)
    uwnd_rc1_dong_mean=where(uwnd_dong_prob1.le.0.1, uwnd_rc1_dong_mean, uwnd_rc1_dong_mean@_FillValue)
    delete(uwnd_xi_prob1)
    delete(vwnd_xi_prob1)
    delete(uwnd_dong_prob1)
    delete(vwnd_dong_prob1)

    psl_rc1_xi_mean=(/psl_rc1_xi_mean/)*10.0
    psl_rc1_dong_mean=(/psl_rc1_dong_mean/)*10.0

    pr_rc1_xi_mean=(/pr_rc1_xi_mean/)*10.0
    pr_rc1_dong_mean=(/pr_rc1_dong_mean/)*10.0

    uwnd_rc1_xi_mean=(/uwnd_rc1_xi_mean/)*10.0
    uwnd_rc1_dong_mean=(/uwnd_rc1_dong_mean/)*10.0

    vwnd_rc1_xi_mean=(/vwnd_rc1_xi_mean/)*10.0
    vwnd_rc1_dong_mean=(/vwnd_rc1_dong_mean/)*10.0

;;;******************************************************************************************************************************
;;;*****************************************************************************************************************************4
  ;;;AMIP实验
    pr_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(pr_xi_mean_amip&time),1),pr_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(pr_xi_mean_amip(0,:,:), pr_rc1_xi_mean_amip)
    delete(pr_xi_mean_amip)
    printVarSummary(pr_rc1_xi_mean_amip)

    pr_rc1_dong_mean_amip=regCoef_n(ispan(1,dimsizes(pr_dong_mean_amip&time),1),pr_dong_mean_amip(:,:,:),0,0)  
    copy_VarCoords(pr_dong_mean_amip(0,:,:), pr_rc1_dong_mean_amip)
    delete(pr_dong_mean_amip)
    printVarSummary(pr_rc1_dong_mean_amip)
 
    uwnd_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(uwnd_xi_mean_amip&time),1),uwnd_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(uwnd_xi_mean_amip(0,:,:), uwnd_rc1_xi_mean_amip)
    delete(uwnd_xi_mean_amip)
    printVarSummary(uwnd_rc1_xi_mean_amip)

     
    uwnd_rc1_dong_mean_amip=regCoef_n(ispan(1,dimsizes(uwnd_dong_mean_amip&time),1),uwnd_dong_mean_amip(:,:,:),0,0)  
    copy_VarCoords(uwnd_dong_mean_amip(0,:,:), uwnd_rc1_dong_mean_amip)
    delete(uwnd_dong_mean_amip)
    printVarSummary(uwnd_rc1_dong_mean_amip)
    
    vwnd_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(vwnd_xi_mean_amip&time),1),vwnd_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(vwnd_xi_mean_amip(0,:,:), vwnd_rc1_xi_mean_amip)
    delete(vwnd_xi_mean_amip)
    printVarSummary(vwnd_rc1_xi_mean_amip)

     
    vwnd_rc1_dong_mean_amip=regCoef_n(ispan(1,dimsizes(vwnd_dong_mean_amip&time),1),vwnd_dong_mean_amip(:,:,:),0,0)  
    copy_VarCoords(vwnd_dong_mean_amip(0,:,:), vwnd_rc1_dong_mean_amip)
    delete(vwnd_dong_mean_amip)
    printVarSummary(vwnd_rc1_dong_mean_amip)


    psl_rc1_xi_mean_amip=regCoef_n(ispan(1,dimsizes(psl_xi_mean_amip&time),1),psl_xi_mean_amip(:,:,:),0,0)  
    copy_VarCoords(psl_xi_mean_amip(0,:,:), psl_rc1_xi_mean_amip)
    printVarSummary(pr_rc1_xi_mean_amip)
    delete(psl_xi_mean_amip)
    psl_rc1_dong_mean_amip=regCoef_n(ispan(1,dimsizes(psl_dong_mean_amip&time),1),psl_dong_mean_amip(:,:,:),0,0)  
    copy_VarCoords(psl_dong_mean_amip(0,:,:), psl_rc1_dong_mean_amip)
    delete(psl_dong_mean_amip)
    printVarSummary(pr_rc1_dong_mean_amip)

    prob1_amip=one_test(psl_rc1_xi_mean_amip)
    prob2_amip=one_test(psl_rc1_dong_mean_amip)
 

    uwnd_xi_prob1=one_test(uwnd_rc1_xi_mean_amip)
    vwnd_xi_prob1=one_test(vwnd_rc1_xi_mean_amip)
    
    uwnd_dong_prob1=one_test(uwnd_rc1_xi_mean_amip)
    vwnd_dong_prob1=one_test(vwnd_rc1_dong_mean_amip)

    uwnd_rc1_xi_mean_amip=where(uwnd_xi_prob1.le.0.1, uwnd_rc1_xi_mean_amip, uwnd_rc1_xi_mean_amip@_FillValue)
    uwnd_rc1_dong_mean_amip=where(uwnd_dong_prob1.le.0.1, uwnd_rc1_dong_mean_amip, uwnd_rc1_dong_mean_amip@_FillValue)

    psl_rc1_xi_mean_amip=(/psl_rc1_xi_mean_amip/)*10.0
    psl_rc1_dong_mean_amip=(/psl_rc1_dong_mean_amip/)*10.0

    pr_rc1_xi_mean_amip=(/pr_rc1_xi_mean_amip/)*10.0
    pr_rc1_dong_mean_amip=(/pr_rc1_dong_mean_amip/)*10.0

    uwnd_rc1_xi_mean_amip=(/uwnd_rc1_xi_mean_amip/)*10.0
    uwnd_rc1_dong_mean_amip=(/uwnd_rc1_dong_mean_amip/)*10.0

    vwnd_rc1_xi_mean_amip=(/vwnd_rc1_xi_mean_amip/)*10.0
    vwnd_rc1_dong_mean_amip=(/vwnd_rc1_dong_mean_amip/)*10.0
;;;*******************************************************************************************************************
;;;**********************************************************************************************************************

;;*********************************************************************************************
;;读入观测的海温和风场的情况
;;;*********************************************************************************************
    yearstart=197901
    yearend=200812

;;;
    print("*****************读入MRE的风场和海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_u_v_chazhi_1000.nc")

    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd_mre_r=f1->uwnd(timeind,{30:-30},{0:360})
    vwnd_mre_r=f1->vwnd(timeind,{30:-30},{0:360})

    delete(uwnd_mre_r&time)
    delete(vwnd_mre_r&time)
    uwnd_mre_r&time=date(timeind)
    vwnd_mre_r&time=date(timeind)
    
    uwnd_mre=month_to_annual(uwnd_mre_r, 1)
    vwnd_mre=month_to_annual(vwnd_mre_r, 1)

    uwnd_mre!0="time"
    vwnd_mre!0="time"
    uwnd_mre&time=ispan(year1, year2,1)
    vwnd_mre&time=ispan(year1, year2,1)


    uwnd_mre&lat@units="degree_north"
    uwnd_mre&lon@units="degree_east"
    vwnd_mre&lat@units="degree_north"
    vwnd_mre&lon@units="degree_east"

    delete(uwnd_mre_r)
    delete(vwnd_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)

  
    ;;;气压场
    ; fils1=systemfunc("ls /mnt/f/share/MRE/mre_slp_chazhi.nc")
    fils1=systemfunc("ls /mnt/f/share/Hadley/slp.mnmean.real.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    ; slp_r1=f1->slp(timeind,{30:-30},{0:360})
    slp_r=f1->slp(timeind,{30:-30},{0:360})
    ; slp_r=tofloat(slp_r1/100.0)
    ; copy_VarMeta(slp_r1, slp_r)
    delete(slp_r&time)
    slp_r&time=date(timeind)
   
    slp_mre=month_to_annual(slp_r, 1)

    slp_mre!0="time"  
    slp_mre&time=ispan(year1, year2,1)

    slp_mre&lat@units="degree_north"
    slp_mre&lon@units="degree_east"
    delete(slp_r)
    delete(date)
    delete(f1)
    delete(fils1)
;;;
  ;;;降水场
    fils1=systemfunc("ls /mnt/f/share/GPCP/precip.mon.mean.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    pr_mre_r=f1->precip(timeind,{30:-30},{0:360})

    delete(pr_mre_r&time)
    printVarSummary(pr_mre_r)
    pr_mre_r&time =  date(timeind)
    pr_mre=month_to_annual(pr_mre_r, 1)
    
    pr_mre!0="time"
    pr_mre&time=ispan(year1, year2,1)
    pr_mre&lat@units="degree_north"
    pr_mre&lon@units="degree_east"

    delete(pr_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)

    ;;;;;*******************************************************************
    ;;;;;;;读入AMIP总的海表气压场，降水场与风场
    ;;;;;********************************************************************
    fils1=systemfunc("ls /mnt/f/share/CMIP5/AMIP/uv/uv_amip_amip_r1i1p1_197901-200812.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd_amipall_r=f1->uwnd(timeind,{30:-30},{0:360})
    vwnd_amipall_r=f1->vwnd(timeind,{30:-30},{0:360})

    delete(uwnd_amipall_r&time)
    delete(vwnd_amipall_r&time)
    uwnd_amipall_r&time=date(timeind)
    vwnd_amipall_r&time=date(timeind)
    
    uwnd_amipall=month_to_annual(uwnd_amipall_r, 1)
    vwnd_amipall=month_to_annual(vwnd_amipall_r, 1)

    uwnd_amipall!0="time"
    vwnd_amipall!0="time"
    uwnd_amipall&time=ispan(year1, year2,1)
    vwnd_amipall&time=ispan(year1, year2,1)


    uwnd_amipall&lat@units="degree_north"
    uwnd_amipall&lon@units="degree_east"
    vwnd_amipall&lat@units="degree_north"
    vwnd_amipall&lon@units="degree_east"

    delete(uwnd_amipall_r)
    delete(vwnd_amipall_r)
    delete(date)
    delete(f1)
    delete(fils1)
    
    ;;;;;;;;;
    ;;;气压场
    fils1=systemfunc("ls /mnt/f/share/CMIP5/AMIP/slp_chazhi/psl_amip_amip_r1i1p1_197901-200812.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    slp_r=lonFlip(f1->psl(timeind,{30:-30},:))
    slp_r=(/slp_r/100.0/)
    delete(slp_r&time)
    slp_r&time=date(timeind)

    slp_amipall=month_to_annual(slp_r, 1)

    slp_amipall!0="time"  
    slp_amipall&time=ispan(year1, year2,1)

    slp_amipall&lat@units="degree_north"
    slp_amipall&lon@units="degree_east"
    delete(slp_r)
    delete(date)
    delete(f1)
    delete(fils1)

    ;;;;;;;;;;
      ;;;降水场
    fils1=systemfunc("ls /mnt/f/share/CMIP5/AMIP/pr_chazhi/pr_amip_amip_r1i1p1_197901-200812.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    pr_amipall_r=f1->pr(timeind,{30:-30},{0:360})
    pr_amipall_r=(/pr_amipall_r*3600.0*24.0/)
    delete(pr_amipall_r&time)
    printVarSummary(pr_amipall_r)
    pr_amipall_r&time =  date(timeind)
    pr_amipall=month_to_annual(pr_amipall_r, 1)
    
    pr_amipall!0="time"
    pr_amipall&time=ispan(year1, year2,1)
    pr_amipall&lat@units="degree_north"
    pr_amipall&lon@units="degree_east"

    delete(pr_amipall_r)
    delete(date)
    delete(f1)
    delete(fils1)














    ;;;
    print("*****************计算趋势************************")  
    pr_rc1_mre=regCoef_n(ispan(1,dimsizes(pr_mre&time),1),pr_mre(:,:,:),0,0)  
    copy_VarCoords(pr_mre(0,:,:), pr_rc1_mre)
    printVarSummary(pr_rc1_mre)
    

    uwnd_rc1_mre=regCoef_n(ispan(1,dimsizes(uwnd_mre&time),1),uwnd_mre(:,:,:),0,0)  
    copy_VarCoords(uwnd_mre(0,:,:), uwnd_rc1_mre)
    printVarSummary(uwnd_rc1_mre)
  
  
    vwnd_rc1_mre=regCoef_n(ispan(1,dimsizes(vwnd_mre&time),1),vwnd_mre(:,:,:),0,0)  
    copy_VarCoords(vwnd_mre(0,:,:), vwnd_rc1_mre)
    printVarSummary(vwnd_rc1_mre)
  
    ;;;***********************************
    slp_rc1_mre=regCoef_n(ispan(1,dimsizes(slp_mre&time),1),slp_mre(:,:,:),0,0)  
    copy_VarCoords(slp_mre(0,:,:), slp_rc1_mre)
    printVarSummary(slp_rc1_mre)
    

    prob_mre=one_test(slp_rc1_mre)
 
    uwnd_mre_prob1=one_test(uwnd_rc1_mre)
    vwnd_mre_prob1=one_test(vwnd_rc1_mre)

    uwnd_rc1_mre=where(uwnd_mre_prob1.le.0.1, uwnd_rc1_mre, uwnd_rc1_mre@_FillValue)
    
    pr_rc1_mre=(/pr_rc1_mre/)*10
    uwnd_rc1_mre=(/uwnd_rc1_mre/)*10
    vwnd_rc1_mre=(/vwnd_rc1_mre/)*10
    slp_rc1_mre=(/slp_rc1_mre/)*10



    ;;;
    print("*****************计算所有AMIP平均的趋势************************")  
    pr_rc1_amipall=regCoef_n(ispan(1,dimsizes(pr_amipall&time),1),pr_amipall(:,:,:),0,0)  
    copy_VarCoords(pr_amipall(0,:,:), pr_rc1_amipall)
    printVarSummary(pr_rc1_amipall)
    

    uwnd_rc1_amipall=regCoef_n(ispan(1,dimsizes(uwnd_amipall&time),1),uwnd_amipall(:,:,:),0,0)  
    copy_VarCoords(uwnd_amipall(0,:,:), uwnd_rc1_amipall)
    printVarSummary(uwnd_rc1_amipall)
  
  
    vwnd_rc1_amipall=regCoef_n(ispan(1,dimsizes(vwnd_amipall&time),1),vwnd_amipall(:,:,:),0,0)  
    copy_VarCoords(vwnd_amipall(0,:,:), vwnd_rc1_amipall)
    printVarSummary(vwnd_rc1_amipall)
  
    ;;;***********************************
    slp_rc1_amipall=regCoef_n(ispan(1,dimsizes(slp_amipall&time),1),slp_amipall(:,:,:),0,0)  
    copy_VarCoords(slp_amipall(0,:,:), slp_rc1_amipall)
    printVarSummary(slp_rc1_amipall)
    

    prob_amipall=one_test(slp_rc1_amipall)
 
    uwnd_amipall_prob1=one_test(uwnd_rc1_amipall)
    vwnd_amipall_prob1=one_test(vwnd_rc1_amipall)

    uwnd_rc1_amipall=where(uwnd_amipall_prob1.le.0.1, uwnd_rc1_amipall, uwnd_rc1_amipall@_FillValue)
    
    pr_rc1_amipall=(/pr_rc1_amipall/)*10
    uwnd_rc1_amipall=(/uwnd_rc1_amipall/)*10
    vwnd_rc1_amipall=(/vwnd_rc1_amipall/)*10
    slp_rc1_amipall=(/slp_rc1_amipall/)*10




 print("**********************画图*****************************")

wks1=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/tos_group/变强变弱模式的集合以及再分析平均结果全球")
 

    ; colormap=read_colormap_file("GMT_jet")
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)

    ; cmap(0:6,:)=colormap(2:92:15,:)
    ; ; cmap(0:6,:)=colormap(2:50:8,:)
    ; cmap(7:17,:)=colormap(105:255:15,:)
   
      cmap=read_colormap_file("MPL_PuOr")
      cmap=(/cmap(::-1,0:3)/)

    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
  
    res@cnFillPalette="GMT_polar"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnFillDrawOrder="PreDraw"
    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4
    res@cnLevelSpacingF      =  0.05

    res@tmXBLabelFontHeightF=0.012
    res@tmYLLabelFontHeightF=0.012
    res@tmYLMajorThicknessF=2.0
    res@tmXBMajorThicknessF=2.0
    res@tmBorderThicknessF=3.0

    res@tmXTOn=False
    res@tmYROn=False
    res@mpFillDrawOrder="PostDraw"
    ; res@mpFillDrawOrder="PreDraw"
    ; res@mpLandFillColor="Transparent"
    res@mpFillOn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@vpWidthF=0.65
    res@vpHeightF=0.35
  

    ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.03
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"

    ;;;CurlyVector
    res@vcLineArrowThicknessF=1.5
    res@vcMinFracLengthF=0.3
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@gsnCenterString ="Observation"
    res@vcRefMagnitudeF          = 0.6
    res@vcMinDistanceF           = 0.03           ; thin vectors\
    res@vcRefLengthF             = 0.04     ; define length of vec ref


    plot=new(6,graphic)

    plot(0) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_mre,vwnd_rc1_mre,slp_rc1_mre,res)
    res@vcRefMagnitudeF          = 0.4 
    res@vcRefLengthF             = 0.04     ; define length of vec ref
    res@gsnCenterString ="All AMIP"
    plot(1) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_amipall,vwnd_rc1_amipall,slp_rc1_amipall,res)
    res@gsnCenterString ="CMIP5(+PWCI)"
    plot(2) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_xi_mean,vwnd_rc1_xi_mean,psl_rc1_xi_mean,res)
    res@gsnCenterString ="CMIP5(-PWCI)"
    plot(3) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_dong_mean,vwnd_rc1_dong_mean,psl_rc1_dong_mean,res)
    res@gsnCenterString ="AMIP(+CMIP5)"
    plot(4) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_xi_mean_amip,vwnd_rc1_xi_mean_amip,psl_rc1_xi_mean_amip,res)
    res@gsnCenterString ="AMIP(-CMIP5)"
    plot(5) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_dong_mean_amip,vwnd_rc1_dong_mean_amip,psl_rc1_dong_mean_amip,res)


         ;;;打点
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@cnLinesOn=False
        res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

        res_pro@cnLineLabelsOn=False
        res_pro@tmYRMode = "Automatic"
        res_pro@cnInfoLabelOn=False

        iplot=new(6,graphic)
        iplot(0)=gsn_csm_contour(wks1,prob_mre,res_pro)
        iplot(1)=gsn_csm_contour(wks1,prob_amipall,res_pro)
        iplot(2)=gsn_csm_contour(wks1,prob1,res_pro)
        iplot(3)=gsn_csm_contour(wks1,prob2,res_pro)
        iplot(4)=gsn_csm_contour(wks1,prob1_amip,res_pro)
        iplot(5)=gsn_csm_contour(wks1,prob2_amip,res_pro)
        
     
        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 1.0         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0025      ; make dots larger
        do i=0,5
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
        ;;;  
        overlay(plot(i), iplot(i))
        end do

        ;;;画slp的等值线
        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -1.0        ; set min contour level
        res_slp@cnMaxLevelValF      = 1.0      ; set max contour level
        res_slp@cnLevelSpacingF     = 0.2      ; set contour spacin

        res_slp@cnLineLabelsOn     =False
        res_slp@cnLineLabelBackgroundColor="Transparent"

        
        colors=new(11,string)
        ; colors(0:3)="darkgoldenrod"
        ; colors(4)="black"
        ; colors(5:)="royalblue2"
         
        colors(0:4)="blueviolet"
        ; colors(0:4)="springgreen"
        ; colors(0:4)="darkgreen"
        
        ; colors(5)="sandybrown"
        colors(5)="Transparent"
        ; colors(5)="white"
        ; colors(6:)="seagreen4"
        ; colors(6:)="springgreen"
        colors(6:)="green"
        ; colors(6:)="deepskyblue"

        res_slp@cnLineColors=colors
        res_slp@cnMonoLineColor=False
        res_slp@gsnContourNegLineDashPattern=12
        res_slp@cnLineLabelInterval        = 2             ; default = 2
        res_slp@cnLineLabelPlacementMode   = "constant"    ; choose constant label method
        res_slp@cnLineLabelFormat= "0@;*.1f" 
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineThicknessF=2.5
        res_slp@gsnContourZeroLineThicknessF=2.0

        plot_slp=new(6,graphic)

        plot_slp(0)=gsn_csm_contour(wks1,pr_rc1_mre ,res_slp)
        plot_slp(1)=gsn_csm_contour(wks1,pr_rc1_amipall ,res_slp)
        plot_slp(2)=gsn_csm_contour(wks1,pr_rc1_xi_mean,res_slp)
        plot_slp(3)=gsn_csm_contour(wks1,pr_rc1_dong_mean,res_slp)
        plot_slp(4)=gsn_csm_contour(wks1,pr_rc1_xi_mean_amip,res_slp)
        plot_slp(5)=gsn_csm_contour(wks1,pr_rc1_dong_mean_amip,res_slp)
        
        do i=0,5
        overlay(plot(i), plot_slp(i))
    end do




pres=True
pres@gsnPanelLabelBar=True
pres@gsnPanelFigureStrings="("+(/"a","b","c","d","e","f"/)+")"
pres@gsnPanelFigureStringsJust="TopLeft"
pres@pmLabelBarWidthF=0.8
gsn_panel(wks1, plot, (/3,2/), pres)


end
    
  

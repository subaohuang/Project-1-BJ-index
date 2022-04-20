
function one_test(ssta_rc)
begin
    
    n=dimsizes(ssta_rc)
    var_ssta_rc    =  dim_variance_n(ssta_rc(:,:,:), 0)    ; 13.477
    var_ssta_rc  =  sqrt(var_ssta_rc/n(0)) 
    var_ssta_rc  =where(var_ssta_rc.eq.0, var_ssta_rc@_FillValue, var_ssta_rc)
    
    ssta_rc_mean=dim_avg_n(ssta_rc, 0)
    t_ssta_rc = ssta_rc_mean/var_ssta_rc
    confi_sst_rc=student_t(t_ssta_rc , n(0)-1)
 return confi_sst_rc

 end
function  sst_cold_annual(sst,timeind,year1,year2)
begin 

    ;;;先得到秋季和冬季的SONDJF海温
    sst_son=month_to_season(sst(time|:,lat|:,lon|:), "SON")
    sst_djf=month_to_season(sst(time|:,lat|:,lon|:), "DJF")
    printVarSummary(sst_son)

    sst_son&time=ispan(year1,year2,1)
    sst_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195lonmin1,195lonmin2的平均
    ;;注意点
    sst_djf(0:dimsizes(sst_djf&time)-2,:,:)=(/sst_djf(1:,:,:)/)
    sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/sst(dimsizes(timeind)-1,:,:)/)
    
    ;;;
    sst_mean_annual=new((/dimsizes(timeind)/12,dimsizes(sst&lat),dimsizes(sst&lon)/),typeof(sst))
    printVarSummary(sst_mean_annual)
    ; sst_mean_annual=month_to_annual(sst_mean, 1)
     ; sst_mean_annual=(/sst_son(:,:,:)+sst_djf(:,:,:)/)*0.5
     sst_mean_annual=(/sst_djf(:,:,:)/)*1.0   ;取冬季
    ; sst_mean_annual=(/month_to_annual(sst(time|:,lat|:,lon|:), 1)/)
    delete(sst_son)
    delete(sst_djf)

    sst_mean_annual!0="time"
    sst_mean_annual!1="lat"
    sst_mean_annual!2="lon"
    sst_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    sst_mean_annual&lon=sst&lon
    sst_mean_annual&lat=sst&lat
    printVarSummary(sst_mean_annual)
return sst_mean_annual
end

;;;********************************************************************
;;;主程序
;;;******************************************************************
begin

xi_model=(/"GFDL-ESM2M","INMCM4","FIO-ESM","FGOALS-s2","GFDL-CM3","IPSL-CM5A-MR"/)
dong_model=(/"BCC-CSM1-1","BCC-CSM1-1-M","MPI-ESM-MR","CMCC-CMS","HadGEM2-CC","MIROC5"/)

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
        "CSIRO-Mk3-6-0","FGOALS-g2","FGOALS-s2","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
         "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
;;;强度
 ;;str model all:

 diri="/mnt/f/share/CMIP5/msf_history/"
 dim=dimsizes(models)
 
 width_trend=asciiread(diri+"cmip_msf_inten_trend.ascii" ,dim, "float")
 
;  ind1=ind(width_trend.ge.0)
;  ind2=ind(width_trend.le.0)

; xi_model=models(ind1)

; dong_model=models(ind2)
delete(models)
print("xi_model="+xi_model)
print("**********")
print("dong_model="+dong_model)



; latmin=-50
; latmax=50
; lonmin=90
; lonmax=300
; lat=fspan(50, -50, 41)
; lon=fspan(90, 300, 85)
; lat_sst=fspan(49.5, -49.5, 100)
; lon_sst=fspan(90.5, 299.5, 210)

; latmin=-30
; latmax=30
; lonmin=90
; lonmax=300
; lat=fspan(30, -30, 25)
; lon=fspan(90, 300, 85)
; lat_sst=fspan(29.5, -29.5, 60)
; lon_sst=fspan(90.5, 299.5, 210)

; ;;;有印度洋
; latmin=-30
; latmax=30
; lonmin=30
; lonmax=300
; lat=fspan(30, -30, 25)
; lon=fspan(30, 300, 109)
; lat_sst=fspan(29.5, -29.5, 60)
; lon_sst=fspan(30.5, 299.5, 270)
; ;;;

;;;全球
latmin=-30
latmax=30
lonmin=0
lonmax=360
lat=fspan(30, -30, 25)
lon=fspan(0, 359.5, 144)
lat_sst=fspan(29.5, -29.5, 60)
lon_sst=fspan(0.5, 359.5, 360)
;;;



psl_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/),float)
psl_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/),float)


sst_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat_sst),dimsizes(lon_sst)/) ,float)
sst_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat_sst),dimsizes(lon_sst)/),float)

uwnd_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/) ,float)
uwnd_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/), float)

vwnd_rc_ximodel=new((/dimsizes(xi_model),dimsizes(lat),dimsizes(lon)/), float)
vwnd_rc_dongmodel=new((/dimsizes(dong_model),dimsizes(lat),dimsizes(lon)/),float)


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
    yearend=200512
    
    year1=1979
    year2=2005
    
     if model.eq."FGOALS-s2".or.model.eq."HadGEM2-CC".or.model.eq."HadGEM2-ES" then

    yearstart=197901
    yearend=200412

    year1=1979
    year2=2004
    end if

    fils1=systemfunc("ls /mnt/i/CMIP5/history/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_historical_r1i1p1_"+"*"+".nc")
    printVarSummary(fils1)
    print("******************"+model+"*********************")
    f1=addfile(fils1, "r")
    delete(fils1)
    vNames= getfilevarnames(f1) ; get variable names of all groups on file
    print(vNames)
    
    if (any(vNames.eq."LONN179_180")) then 


    time=f1->TIME
 
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    if model.eq."FGOALS-s2" then
    sst_r=f1->TOS(timeind,:,:)
    else
    sst_r=lonFlip(f1->TOS(timeind,:,:))
    end if    
    sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r)
    dNames = getvardims(f1->TOS)
      delete(f1)
    print(dNames)
    sst_r!0="time"
    sst_r!1="lat"
    sst_r!2="lon"
    end if

  
    delete(vNames)
    delete(dNames)

    delete(sst_r&time)
    sst_r&time=date(timeind)
      if model.eq."FGOALS-s2".or.model.eq."HadGEM2-CC".or.model.eq."HadGEM2-ES" then
    sst=sst_cold_annual(sst_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2004)
    else 
    sst=sst_cold_annual(sst_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    end if
    delete(sst_r)
    delete(time)
    delete(date)
    printVarSummary(sst)
    delete(timeind)
 
    print("******************读入风场*********************")
    fils1=systemfunc("ls /mnt/i/CMIP5/history/uv/uv_"+model+"_history_r1i1p1_"+"*"+".nc")

    print(fils1)
    f1=addfile(fils1, "r")

    time=f1->time    
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    
    uwnd_r=f1->uwnd(timeind,:,:)  ;;;
    vwnd_r=f1->vwnd(timeind,:,:)  ;;;
    delete(f1)

;;;;
    delete(uwnd_r&time)
    delete(vwnd_r&time)
    uwnd_r&time=date(timeind)
    vwnd_r&time=date(timeind)
     if model.eq."FGOALS-s2".or.model.eq."HadGEM2-CC".or.model.eq."HadGEM2-ES" then
    uwnd=sst_cold_annual(uwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2004)
    vwnd=sst_cold_annual(vwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2004)
    else 
    uwnd=sst_cold_annual(uwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    vwnd=sst_cold_annual(vwnd_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    end if
    delete(timeind)
    delete(uwnd_r)
    delete(vwnd_r)
    delete(time)
    delete(date)

   print("******************读入海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/i/CMIP5/history/slp_chazhi/psl_"+model+"_historical_r1i1p1_"+"*"+".nc")
    print(fils1)
    f1=addfile(fils1, "r")
   
    time=f1->time

    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    

      ; psl_r=f1->PSL(timeind,:,:)
      psl_r=lonPivot(f1->psl(timeind,:,:),0)
 
     print(psl_r&lon)
      psl_r=(/psl_r/)/100.0
    printVarSummary(psl_r)
    delete(psl_r&time)
    psl_r&time=date(timeind)
    
    lonind=ind(psl_r&lon.le.lonmax.and.psl_r&lon.ge.lonmin)

     if model.eq."FGOALS-s2".or.model.eq."HadGEM2-CC".or.model.eq."HadGEM2-ES" then
    psl=sst_cold_annual(psl_r(:,{latmax:latmin},lonind),timeind,1979,2004)
    else
    psl=sst_cold_annual(psl_r(:,{latmax:latmin},lonind),timeind,1979,2005)
end if
    delete(timeind)
    delete(psl_r)
    delete(time)
    delete(date)
    delete(lonind)

    print("*****************计算趋势************************")  

    sst_rc1=regCoef_n(ispan(1,dimsizes(sst&time),1),sst(:,:,:),0,0)  
    copy_VarCoords(sst(0,:,:), sst_rc1)
    delete(sst)
    printVarSummary(sst_rc1)
 
    

    uwnd_rc1=regCoef_n(ispan(1,dimsizes(uwnd&time),1),uwnd(:,:,:),0,0)  
    copy_VarCoords(uwnd(0,:,:), uwnd_rc1)
    delete(uwnd)
    printVarSummary(uwnd_rc1)

    
   
    vwnd_rc1=regCoef_n(ispan(1,dimsizes(vwnd&time),1),vwnd(:,:,:),0,0)  
    copy_VarCoords(vwnd(0,:,:), vwnd_rc1)
    delete(vwnd)
    printVarSummary(vwnd_rc1)

    
    psl_rc1=regCoef_n(ispan(1,dimsizes(psl&time),1),psl(:,:,:),0,0)  
    copy_VarCoords(psl(0,:,:), psl_rc1)
    delete(psl)
    printVarSummary(psl_rc1)   

    
    if i.eq.0 then 
      psl_rc_ximodel(j,:,:)=psl_rc1
      sst_rc_ximodel(j,:,:)=sst_rc1
      uwnd_rc_ximodel(j,:,:)=uwnd_rc1
      vwnd_rc_ximodel(j,:,:)=vwnd_rc1
   else
      psl_rc_dongmodel(j,:,:)=psl_rc1
      sst_rc_dongmodel(j,:,:)=sst_rc1
      uwnd_rc_dongmodel(j,:,:)=uwnd_rc1
      vwnd_rc_dongmodel(j,:,:)=vwnd_rc1
    end if

    delete([/psl_rc1,sst_rc1,uwnd_rc1,vwnd_rc1/])
end do
delete(dim)
delete(models)
end do

psl_rc1_xi_mean=dim_avg_n(psl_rc_ximodel,0)
sst_rc1_xi_mean=dim_avg_n(sst_rc_ximodel,0)
uwnd_rc1_xi_mean=dim_avg_n(uwnd_rc_ximodel,0)
vwnd_rc1_xi_mean=dim_avg_n(vwnd_rc_ximodel,0)

psl_rc1_dong_mean=dim_avg_n(psl_rc_dongmodel,0)
sst_rc1_dong_mean=dim_avg_n(sst_rc_dongmodel,0)
uwnd_rc1_dong_mean=dim_avg_n(uwnd_rc_dongmodel,0)
vwnd_rc1_dong_mean=dim_avg_n(vwnd_rc_dongmodel,0)


psl_rc1_xi_mean!0="lat"
psl_rc1_xi_mean!1="lon"

psl_rc1_xi_mean&lat=lat
psl_rc1_xi_mean&lon=lon
psl_rc1_xi_mean&lat@units="degrees_north"
psl_rc1_xi_mean&lon@units="degrees_east"
printVarSummary(psl_rc1_xi_mean)

sst_rc1_xi_mean!0="lat"
sst_rc1_xi_mean!1="lon"

sst_rc1_xi_mean&lat=lat_sst
sst_rc1_xi_mean&lon=lon_sst
sst_rc1_xi_mean&lat@units="degrees_north"
sst_rc1_xi_mean&lon@units="degrees_east"

printVarSummary(psl_rc1_xi_mean)


copy_VarCoords(psl_rc1_xi_mean, psl_rc1_dong_mean)
copy_VarCoords(sst_rc1_xi_mean, sst_rc1_dong_mean)
copy_VarCoords(psl_rc1_xi_mean, uwnd_rc1_xi_mean)
copy_VarCoords(psl_rc1_xi_mean, uwnd_rc1_dong_mean)
copy_VarCoords(psl_rc1_xi_mean, vwnd_rc1_xi_mean)
copy_VarCoords(psl_rc1_xi_mean, vwnd_rc1_dong_mean)

; prob1=one_test(sst_rc_ximodel)
; copy_VarCoords(sst_rc1_xi_mean, prob1)
; printVarSummary(prob1)

; prob2=one_test(sst_rc_dongmodel)
; copy_VarCoords(sst_rc1_dong_mean, prob2)

; uwnd_prob1=one_test(uwnd_rc_ximodel)
; vwnd_prob1=one_test(vwnd_rc_ximodel)
; uwnd_rc1_xi_mean=where(uwnd_prob1.le.0.1.or.vwnd_prob1.le.0.1, uwnd_rc1_xi_mean, uwnd_rc1_xi_mean@_FillValue)

; uwnd_prob2=one_test(uwnd_rc_ximodel)
; vwnd_prob2=one_test(vwnd_rc_dongmodel)
; uwnd_rc1_dong_mean=where(uwnd_prob2.le.0.2.or.vwnd_prob2.le.0.2, uwnd_rc1_dong_mean, uwnd_rc1_dong_mean@_FillValue)

psl_rc1_xi_mean=(/psl_rc1_xi_mean/)*10.0
psl_rc1_dong_mean=(/psl_rc1_dong_mean/)*10.0

sst_rc1_xi_mean=(/sst_rc1_xi_mean/)*10.0
sst_rc1_dong_mean=(/sst_rc1_dong_mean/)*10.0

uwnd_rc1_xi_mean=(/uwnd_rc1_xi_mean/)*10.0
uwnd_rc1_dong_mean=(/uwnd_rc1_dong_mean/)*10.0

vwnd_rc1_xi_mean=(/vwnd_rc1_xi_mean/)*10.0
vwnd_rc1_dong_mean=(/vwnd_rc1_dong_mean/)*10.0
;;*********************************************************************************************
;;读入观测的海温和风场的情况
;;;*********************************************************************************************
   yearstart=197901
    yearend=200512
;;;
   print("*****************读入MRE的风场和海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_u_v_chazhi_1000.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    uwnd_mre_r=f1->uwnd(timeind,:,:)
    vwnd_mre_r=f1->vwnd(timeind,:,:)

    delete(uwnd_mre_r&time)
    delete(vwnd_mre_r&time)
    uwnd_mre_r&time=date(timeind)
    vwnd_mre_r&time=date(timeind)
    uwnd_mre=sst_cold_annual(uwnd_mre_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    vwnd_mre=sst_cold_annual(vwnd_mre_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
     
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
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_slp_chazhi.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    slp_r=f1->slp(timeind,:,:)

    delete(slp_r&time)
    slp_r&time=date(timeind)
    slp_mre=sst_cold_annual(slp_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
      
    slp_mre&lat@units="degree_north"
    slp_mre&lon@units="degree_east"
    delete(slp_r)
    delete(date)
    delete(f1)
    delete(fils1)
;;;
  ;;;海温场
    fils1=systemfunc("ls /mnt/f/share/MRE/haderscobe_sst_chazhi.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_mre_r=f1->sst(timeind,:,:)

    delete(sst_mre_r&time)
    printVarSummary(sst_mre_r)
    sst_mre_r&time =  date(timeind)
    sst_mre=sst_cold_annual(sst_mre_r(:,{latmax:latmin},{lonmin:lonmax}),timeind,1979,2005)
    sst_mre=(/month_to_annual(sst_mre_r(:,{latmax:latmin},{lonmin:lonmax}), 1)/)

    sst_mre&lat@units="degree_north"
    sst_mre&lon@units="degree_east"
    delete(sst_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)
    sst_mre_cha= dim_avg_n_Wrap(sst_mre(:,{-5:5},{80:160}),(/1,2/)) - dim_avg_n_Wrap(sst_mre(:,{-5:5},{200:280}),(/1,2/))
    
    sst_cha_trend=regCoef_n(ispan(1,dimsizes(sst_mre&time),1),sst_mre_cha,0,0)  
    sst_cha_trend=(/sst_cha_trend*10/)
    print(sst_cha_trend)
    exit()

    ;;;
    print("*****************计算趋势************************")  

    sst_rc1_mre=regCoef_n(ispan(1,dimsizes(sst_mre&time),1),sst_mre(:,:,:),0,0)  
    copy_VarCoords(sst_mre(0,:,:), sst_rc1_mre)
    printVarSummary(sst_rc1_mre)
    
    sst_rc1_mre=(/sst_rc1_mre/)*10

    uwnd_rc1_mre=regCoef_n(ispan(1,dimsizes(uwnd_mre&time),1),uwnd_mre(:,:,:),0,0)  
    copy_VarCoords(uwnd_mre(0,:,:), uwnd_rc1_mre)
    printVarSummary(uwnd_rc1_mre)
  
    uwnd_rc1_mre=(/uwnd_rc1_mre/)*10
      
    vwnd_rc1_mre=regCoef_n(ispan(1,dimsizes(vwnd_mre&time),1),vwnd_mre(:,:,:),0,0)  
    copy_VarCoords(vwnd_mre(0,:,:), vwnd_rc1_mre)
    printVarSummary(vwnd_rc1_mre)
    vwnd_rc1_mre=(/vwnd_rc1_mre/)*10
    ;;;***********************************
    slp_rc1_mre=regCoef_n(ispan(1,dimsizes(slp_mre&time),1),slp_mre(:,:,:),0,0)  
    copy_VarCoords(slp_mre(0,:,:), slp_rc1_mre)
    printVarSummary(slp_rc1_mre)
    slp_rc1_mre=(/slp_rc1_mre/)*10


 print("**********************画图*****************************")

wks1=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/tos_group/变强变弱模式的集合以及再分析平均结果全球")
    ; colormap=read_colormap_file("BlWhRe")
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)
    ; cmap(0:4,:)=colormap(0:40:10,:)
    ; cmap(5:17,:)=colormap(52:100:4,:)

    colormap=read_colormap_file("GMT_jet")
    printVarSummary(colormap)
    cmap=new((/18,4/),typeof(colormap))
    printVarSummary(colormap)

    cmap(0:6,:)=colormap(2:92:15,:)
    ; cmap(0:6,:)=colormap(2:50:8,:)
    cmap(7:17,:)=colormap(105:255:15,:)

     
  ; colormap=read_colormap_file("GMT_jet")
  ;   printVarSummary(colormap)
  ;   cmap=new((/20,4/),typeof(colormap))
  ;   printVarSummary(colormap)
  ;   cmap(0:6,:)=colormap(2:92:15,:)
  ;   cmap(7:19,:)=colormap(106:238:11,:)


    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterString ="CMIP5(+PWC trend)"
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
  
    res@cnFillPalette=cmap
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnFillDrawOrder="PreDraw"
   
    
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False
    
    res@mpFillDrawOrder="PostDraw"
    
    res@mpLandFillColor="gray40"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@vpWidthF=0.7
    res@vpHeightF=0.35
  

    ;;;LineArrow
    res@gsnScalarContour        = True
    res@vcRefMagnitudeF          = 0.4 
    res@vcRefLengthF             = 0.03     ; define length of vec ref
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.05
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01


    ;;;CurlyVector

    res@vcLineArrowThicknessF=1.5
    res@vcMinFracLengthF=0.3
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
        


    plot=new(3,graphic)
    plot(0) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_xi_mean,vwnd_rc1_xi_mean,sst_rc1_xi_mean,res)
     res@gsnCenterString ="CMIP5(-PWC trend)"
    plot(1) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_dong_mean,vwnd_rc1_dong_mean,sst_rc1_dong_mean,res)
     res@gsnCenterString ="Observation "
     res@vcRefMagnitudeF          = 0.6
    plot(2) = gsn_csm_vector_scalar_map(wks1,uwnd_rc1_mre,vwnd_rc1_mre,sst_rc1_mre,res)

        ;  ;;;打点
        ; res_pro=True
        ; res_pro@gsnDraw=False
        ; res_pro@gsnFrame=False
        ; res_pro@cnLinesOn=False
        ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        ; res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

        ; res_pro@cnLineLabelsOn=False
        ; res_pro@tmYRMode = "Automatic"
        ; res_pro@cnInfoLabelOn=False

        ; iplot=new(2,graphic)

        ; iplot(0)=gsn_csm_contour(wks1,prob1,res_pro)
        ; iplot(1)=gsn_csm_contour(wks1,prob2,res_pro)
     
        ; opt                      = True
        ; opt@gsnShadeFillType     = "pattern"
        ; opt@gsnShadeMid     = 17         ; stipple
        ; ; opt@gsnShadeHigh   = "green"
        ; opt@gsnShadeFillScaleF   = 0.2         ; add extra density
        ; opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
        ; do i=0,1
        ; iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
        ; ;;;  
        ; overlay(plot(i), iplot(i))
        ; end do

        ;;;画slp的等值线
        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -0.8        ; set min contour level
        res_slp@cnMaxLevelValF      = 0.8        ; set max contour level
        res_slp@cnLevelSpacingF     = 0.05       ; set contour spacing
        res_slp@cnLineLabelsOn     =True
        res_slp@cnLineLabelBackgroundColor="Transparent"
        res_slp@gsnContourNegLineDashPattern=2
        res_slp@cnLineLabelInterval        = 2             ; default = 2
        res_slp@cnLineLabelPlacementMode   = "constant"    ; choose constant label method
        res_slp@cnLineLabelFormat= "0@;*.1f" 
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineColor="black"
        res_slp@cnLineThicknessF=1.5

        plot_slp=new(3,graphic)

        
        plot_slp(0)=gsn_csm_contour(wks1,psl_rc1_xi_mean,res_slp)
        plot_slp(1)=gsn_csm_contour(wks1,psl_rc1_dong_mean,res_slp)
        res_slp@cnLevelSpacingF     = 0.1      ; set contour spacing
        plot_slp(2)=gsn_csm_contour(wks1,slp_rc1_mre ,res_slp)

        do i=0,2
        overlay(plot(i), plot_slp(i))
    end do




pres=True
pres@gsnPanelLabelBar=True
pres@gsnPanelFigureStrings="("+(/"a","b","c"/)+")"
pres@gsnPanelFigureStringsJust="TopLeft"
gsn_panel(wks1, plot, (/3,1/), pres)


end
    
  

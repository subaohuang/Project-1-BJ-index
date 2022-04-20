begin
;读入数据
; model_name = (/"CESM2",  "EC-Earth3"  , "GFDL-ESM4"  ,"MPI-ESM1-2-HR" ,\
;  "CESM2-WACCM" , "EC-Earth3-Veg",  "GISS-E2-1-G", "IPSL-CM6A-LR" ,"MRI-ESM2-0","SAM0-UNICON",\
; "BCC-ESM1","CNRM-CM6-1","FGOALS-f3-L","GISS-E2-1-G-CC","MCM-UA-1-0","NESM3"/)
; model_name = (/ "ACCESS1-0","ACCESS1-3", "CanESM2","CMCC-CESM","CMCC-CM","CMCC-CMS",\
;             "CSIRO-Mk3-6-0", "FGOALS-s2","GFDL-CM2p1","GFDL-CM3" ,"GFDL-ESM2G","GFDL-ESM2M",\
;             "GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadCM3","HadGEM2-CC","HadGEM2-ES",\
;             "inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR",\
;             "MPI-ESM-P","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
model_name = (/ "MPI-ESM-MR",\
            "MPI-ESM-P","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)


dim_model = dimsizes(model_name)

a=6371000
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180

diri_output="/home/yangsong3/data-observation/linshh/CMIP5-derived/msf/1950-2004/"
do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart = 195001
    yearend = 200412
    fils1 := "/home/yangsong3/data-observation/linshh/CMIP5-derived/vd/1950-2004/ud_vd_"+model+"_history_r1i1p1_195001-200412.nc"
    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    print(date)
 
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    vd := f1->vd(timeind,:,:,:)
    
    printVarSummary(vd)
    
    ;;;为保证垂直方向上的质量守恒，我们先去掉垂直方向质量平均的值
    ps = 1013.0
    plev = vd&level
    dp = dpres_plevel(plev,ps,100.0,0)
    printVarSummary(dp)
    print(plev)
    print(dp)

    vd_avg = vd
    dp_4d = conform_dims(dimsizes(vd), dp, 1)
    vd_avg = (/vd*dp_4d/sum(dp)/)
    delete(dp_4d)
    delete(dp)
    printMinMax(vd_avg, 1)
    vd = (/vd - vd_avg/)
    printMinMax(vd, 1)
    delete(vd_avg)
    ;;;


    ;;;对四维的vd进行垂直积分
    ; print("开始积分")
    ; pbot=dble2flt(vd&level)
    ; pres_new=new((/dimsizes(vd&time),dimsizes(vd&lat),dimsizes(vd&lon)/),float)
    ; pres_new=1013.0
    ; vint=new((/dimsizes(vd&level),dimsizes(vd&time),dimsizes(vd&lat),dimsizes(vd&lon)/),float)

  
    ; do i=0,(dimsizes(vd&level)-1)
    ; vint(i,:,:,:)=dble2flt(vibeta(vd&level, vd(time|:,lat|:,lon|:,level|:), 1, pres_new,pbot(i),0))
    ; end do
    ; printVarSummary(vint)
    ; exit()

    ;;;;对四维的vd进行垂直积分



    ;;;先得到层厚和每层厚度的风场
    

    level2 = plev
    level2(0:dimsizes(plev)-2) = plev(1:)
    level2(dimsizes(plev)-1) = 0
    dp = plev - level2
    
    delete(plev)
    delete(level2)

    vd_2 = vd
    vd_in = vd

    vd_2(:,0:dimsizes(vd&level)-2,:,:) = vd(:,1:,:,:)
    vd_2(:,dimsizes(vd&level)-1,:,:) = 0.0
    vd_in = (/(vd + vd_2)*0.5/)
    printVarSummary(vd_in)
    delete(vd_2)


   ;;; 积分
    vint = vd_in
    do i=0,dimsizes(vd_in&level)-1
        print("积分开始")
        dp_new = conform_dims(dimsizes(vd_in(:,i:,:,:)), dp(i:), 1)
        vint(:,i,:,:)= dim_sum_n(vd_in(:,i:,:,:)*dp_new,1)
        delete(dp_new)
    end do 
    delete(dp)
    delete(vd_in)
    printVarSummary(vint)

    cosfi = cos((pi/180.0)*vd&lat)

    cosfi_new = conform_dims(dimsizes(vint), cosfi, 2)

    delete(cosfi)
    
    msf = lonFlip(vint)
    delete(vint)
    msf =(/tofloat((/((2.0*pi*a)/g)*cosfi_new*msf*100.0/))/)
    delete(cosfi_new)
    
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,{240:270}), (/3/))
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,:), (/3/))
    ; msf_mean_clm = dim_avg_n_Wrap(msf_mean, (/0/))
    ; msf_mean_clm = (/msf_mean_clm/10^10/)
; ;;;;;draw 
;     res = True
;     res@gsnDraw=False
;     res@gsnFrame=False
;     res@cnFillPalette = "BlRe"

;     res@cnLinesOn=True
 
;     res@cnFillOn =True
   
;     res@cnLinesOn=True
;     res@cnLevelSelectionMode ="ManualLevels"
;     res@cnMinLevelValF = -12
;     res@cnMaxLevelValF = 12
;     res@cnLevelSpacingF = 2

;     res@lbLabelBarOn =True
;     res@cnLineLabelsOn = False
;     res@cnInfoLabelOn = False
  
    
;     res@tiYAxisString=""
;     res@tmXTOn=False
;     res@tmYROn=False
;     res@tmXBMode="Explicit"
 
;     res@tmXBMinorOn=True
;     res@tmXBMinorValues=ispan(-90,90,10)

  
;     res@tmYRMode="Automatic"
     
;     res@vpHeightF=0.3
;     res@vpWidthF=0.7

;     ;;;;
;     wks2 = gsn_open_wks("pdf", "CMIP6中哈德莱环流气候态")
    
;     res@gsnLeftString = "Spring"
;     plot = gsn_csm_pres_hgt(wks2,msf_mean_clm,res)
; draw(plot)
; frame(wks2)
; exit()

    ntime = dimsizes(msf&time)
    nlev = dimsizes(msf&level)
    nlat = dimsizes(msf&lat)
    nlon = dimsizes(msf&lon)


    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="msf_hc_"+model+"_history_r1i1p1_195001-200412.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "CMIP5 historical run Mass stream Function for Hadley"    
    fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time",typeof(msf&time), (/"time"/))
    filevardef(fout, "level",typeof(msf&level), (/"level"/))
    filevardef(fout, "lat", typeof(msf&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(msf&lon),(/"lon"/))

   
    filevardef(fout, "msf", typeof(msf),(/"time","level","lat","lon"/))

   
    filevarattdef(fout, "level", msf&level)
    filevarattdef(fout, "lat", msf&lat)
    filevarattdef(fout, "lon", msf&lon)
    filevarattdef(fout, "time", msf&time)
    filevarattdef(fout, "msf", msf)
 
    
    setfileoption(fout,"DefineMode",False)


    fout->time=(/msf&time/)
    fout->level=(/msf&level/)
    fout->lat=(/msf&lat/)
    fout->lon=(/msf&lon/)
    fout->msf = (/msf/)
    delete(msf)
 print("********************************finish calculate"+" "+model+"************************************")
 end do 
end 




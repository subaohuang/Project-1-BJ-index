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
model_name = (/"GFDL-CM3" ,"GFDL-ESM2G","GFDL-ESM2M",\
            "GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadCM3","HadGEM2-CC","HadGEM2-ES",\
            "inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR",\
            "MPI-ESM-P","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)


dim_model = dimsizes(model_name)

a=6371000.0
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180.0

diri_output="/home/yangsong3/data-observation/linshh/CMIP5-derived/msf_wc/1850-1949/"

do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart = 185001
    yearend = 194912
    fils1 := "/home/yangsong3/data-observation/linshh/CMIP5-derived/vd/1850-1949/ud_vd_"+model+"_history_r1i1p1_185001-194912.nc"
    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
  
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    ud := f1->ud(timeind,:,{-30:30},:)
    plev = ud&level
    
    printVarSummary(ud)

     
    ; ; ;;;;;;方法1
    ;     ud_new = ud(time|:,lat|:,lon|:,level|:)  
    ;     delete(ud)
    ;     pbot = ud_new&level

    ;     pres_new = new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    ;     pres_new = 1013.0
    ;     vint = new((/dimsizes(ud_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    ;     do i = 0,(dimsizes(ud_new&level)-1)
    ;         print("积分开始")
    ;       vint(i,:,:,:) = vibeta(ud_new&level, ud_new, 1, pres_new,pbot(i),0)    ;;;intergrate from specific level to 0 ,the order of ud_new's level must 
    ;                                                                     ;;; Bottom to top
    ;     end do
    ; ;;;;;
    
 ;;;先得到层厚和每层厚度的风场，方法2
    level2 = plev
    level2(0:dimsizes(plev)-2) = plev(1:)
    level2(dimsizes(plev)-1) = 0
    dp = plev - level2

    delete(level2)

    ud_2 = ud
    ud_in = ud

    ud_2(:,0:dimsizes(ud&level)-2,:,:) = ud(:,1:,:,:)
    ud_2(:,dimsizes(ud&level)-1,:,:) = 0.0
    ud_in = (/(ud + ud_2)*0.5/)
    printVarSummary(ud_in)
    delete(ud_2)
    ;;;;

    ;;; 积分
    vint = ud_in
    do i=0,dimsizes(ud_in&level)-1
        print("积分开始")
        dp_new = conform_dims(dimsizes(ud_in(:,i:,:,:)), dp(i:), 1)
        vint(:,i,:,:)= dim_sum_n(ud_in(:,i:,:,:)*dp_new,1)
        delete(dp_new)
    end do 
   delete(ud_in)
    fi_time = ((a*dfi)/g)*vint*100.0
    
    copy_VarCoords(vint, fi_time)
    delete(vint)
    printVarSummary(fi_time)


    ; fi_time!0 = "level"
    ; fi_time!1 = "time"
    ; fi_time!2 = "lat"
    ; fi_time!3 = "lon"

    ; fi_time&level = ud_new&level
    ; fi_time&level@units = "hpa"
    ; fi_time&lon = ud_new&lon
    ; fi_time&time = ud_new&time
    ; fi_time&lat = ud_new&lat
    msf = lonFlip(fi_time(time|:,level|:,lat|:,lon|:))
    printMinMax(msf,1)
    delete(fi_time)

    ; msf = lonFlip(fi_time)
    printVarSummary(msf)
 
    ntime=dimsizes(msf&time)
    nlev=dimsizes(msf&level)
    nlat=dimsizes(msf&lat)
    nlon=dimsizes(msf&lon)


;     msf_mean_clm = dim_avg_n_Wrap(msf(time|:,{lat|-5:5},level|:,lon|:), (/0,1/))

;     printMinMax(msf, 1)

;     msf_mean_clm = (/msf_mean_clm/10^9/)
;     ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,:), (/3/))
;     ; msf_mean_clm = dim_avg_n_Wrap(msf_mean, (/0/))
;     ; msf_mean_clm = (/msf_mean_clm/10^10/)
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
 
;     ; res@tmXBMinorOn=True
;     ; res@tmXBMinorValues=ispan(-90,90,10)

  
;     res@tmYRMode="Automatic"
     
;     res@vpHeightF=0.3
;     res@vpWidthF=0.7

;     ;;;;
;     wks2 = gsn_open_wks("pdf", "CMIP5中沃克环流气候态")
    
;     res@gsnLeftString = "Spring"
;     plot = gsn_csm_pres_hgt(wks2,msf_mean_clm,res)
; draw(plot)
; frame(wks2)
; exit()




    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="msf_wc_"+model+"_history_r1i1p1_185001-194912.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "CMIP5 historical run Mass stream Function for walker circulation"    
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




 ;;;************************************************************
  ;; 计算质量流函数
  ;; 作为funcion调用计算MSF
  ;; using function 
  ;; input values must be in ascending latitude order
  ;; input array must be on a global grid
  ;; input data must  have the dimension name of lat ,lon ,time ,lev,can use cdo setname 
  ;; to change ,the data variable must be named U,V,OMEGA
  ;; 
  ;;;********************************************************、
function get_msf(path1,path2,timestart,timeend,outpath,outputname)
local a, g, pi ,dfi
  begin
   ;;;流函数公式需要的参数
    a = 6371000
    g = 9.8
    pi = get_pi("double")
    dfi = (pi*10)/180


    f1 = addfile(path1, "r")
    f2 = addfile(path2, "r")
    time = f1->time
    ; time = (/time - 1/)

    lat = f1->lat
    lon = f1->lon
    lev= f1->level   ;level反着取

    lev@units = "hpa"
    lev@long_name = "pressure level"

    date = cd_calendar(time, 1)    ;;;output calendar date 

    timeind = ind(date.ge.timestart.and.date.le.timeend)
    time_wr = time(timeind)

    levwant= (/1000.0, 925.0, 850.0, 700.0, 600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)   ;level反着取


    levind1 = get1Dindex(lev,levwant)

    level  = f1->level(levind1)

    delete(date)     

   
    vwnd = f2->vwnd(timeind,levind1,:,:)   ;;;variable must have ascending lat and lon order coordinates
    ; uwnd = f1->uwnd(timeind,levind1,:,:)   ;;;variable must have ascending lat and lon order coordinates
 
    ; printVarSummary(uwnd)
    ; ;;;
    ; print("开始计算辐散风")
    ; vrdvx = uv2vrdvF(uwnd,vwnd)
    ; uvd = dv2uvF(vrdvx(1,:,:,:,:))
    ; delete(vrdvx)
    ; vd = uvd(1,:,:,:,:)
    ; print("结束计算辐散风")

    ; copy_VarCoords(vwnd, vd)
    ; copy_VarAtts(vwnd, vd)
    vd = vwnd
    vd!1 = "lev"

    ;;;先得到层厚和每层厚度的风场
    ps = 1013.0
    plev = vwnd&level

    dp = dpres_plevel(plev,ps,min(plev),0)
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
    ;

    level2 = plev
    level2(0:dimsizes(plev)-2) = plev(1:)
    level2(dimsizes(plev)-1) = 0
    dp = plev - level2
    
    delete(plev)
    delete(level2)

    vd_2 = vd
    vd_in = vd

    vd_2(:,0:dimsizes(vd&lev)-2,:,:) = vd(:,1:,:,:)
    vd_2(:,dimsizes(vd&lev)-1,:,:) = 0.0
    vd_in = (/(vd + vd_2)*0.5/)
    printVarSummary(vd_in)
    delete(vd_2)
    
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
  
   ;;; 积分
    vint = vd_in
    do i=0,dimsizes(vd_in&lev)-1
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
    nlev = dimsizes(msf&lev)
    nlat = dimsizes(msf&lat)
    nlon = dimsizes(msf&lon)

    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro = outpath
    filo1 = outputname
    system("/bin/rm -f "+diro + filo1+"")
    ; setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
    ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "HC MSF of CESM "    
    fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copy file attributes    
 
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lev","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time",typeof(msf&time), (/"time"/))
    filevardef(fout, "lev",typeof(level), (/"lev"/))
    filevardef(fout, "lat", typeof(msf&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(msf&lon),(/"lon"/))

   
    filevardef(fout, "msf", typeof(msf),(/"time","lev","lat","lon"/))
    

   
    filevarattdef(fout, "lev", level)
    filevarattdef(fout, "lat", msf&lat)
    filevarattdef(fout, "lon", msf&lon)
    filevarattdef(fout, "time", msf&time)
    filevarattdef(fout, "msf", msf)
    
 
    
    setfileoption(fout,"DefineMode",False)


    fout->time=(/msf&time/)
    fout->lev=(/msf&lev/)
    fout->lat=(/msf&lat/)
    fout->lon=(/msf&lon/)
    fout->msf = (/msf/)
    
    delete(msf)
   echo = "msf creat successfully"
return(echo)
end 

begin

;function(path,filename,timestart,timeend,outpath,outputname)
  
       inpath1="/home/yangsong3/data-observation/NCEP_NCAR/monthly/uwnd.mon.mean.nc"
       inpath2="/home/yangsong3/data-observation/NCEP_NCAR/monthly/vwnd.mon.mean.nc"
       
       outputpath="/home/yangsong3/data-observation/linshh/obs-derived/msf/"
       outputname="msf_HC_NCEP1_197901-202012-test.nc"
  m = get_msf(inpath1,inpath2,197901,202012,outputpath,outputname)
  print(m)

end 




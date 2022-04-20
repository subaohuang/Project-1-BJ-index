
;;;********************************************************************
;;;主程序
;;;******************************************************************
begin

; models=(/"CSIRO-Mk3-6-0","INMCM4","GFDL-CM3","Ipr-CM5A-MR"/)
; dong_model=(/"MIROC5","MPI-ESM-MR"/)

    models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)


    ;;;全球
    lat_pr=fspan(89.5, -89.5, 71)
    lon_pr=fspan(0.5, 359.5, 144)
    ;;;

    pr_ximodel=new((/dimsizes(models),360,dimsizes(lat_pr),dimsizes(lon_pr)/),float)
    dim=dimsizes(models)

    diri_output="/mnt/f/share/CMIP5/AMIP/pr_chazhi/"

 do j=0,dim-1
    
    model=models(j)

    yearstart=197901
    yearend=200812
    
    year1=1979
    year2=2008
    
    print("******************"+model+"*********************")
    fils1=systemfunc("ls /mnt/f/share/CMIP5/AMIP/pr_chazhi/pr_"+model+"_amip_r1i1p1_197901-200812.nc")
    print(fils1)
    f1=addfile(fils1, "r")
    time=f1->time
    date=cd_calendar(time, 1)
    delete(fils1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    delete(date)
    pr_r=f1->pr(timeind,:,:)
    printVarSummary(pr_r)
    pr_ximodel(j,:,:,:)=pr_r
    delete(pr_r)
    end do
    
    pr_xi_mean=dim_avg_n(pr_ximodel,0)
    copy_VarMeta(pr_ximodel(0,:,:,:), pr_xi_mean)
    printVarSummary(pr_xi_mean)

    ntime=dimsizes(pr_xi_mean&time)
    nlat=dimsizes(pr_xi_mean&lat)
    nlon=dimsizes(pr_xi_mean&lon)
    printVarSummary(pr_xi_mean)





 print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="pr_amip_amip_r1i1p1_197901-200812.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "pr插值后的结果，取1000hpa"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(pr_xi_mean&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(pr_xi_mean&lon),(/"lon"/))
    filevardef(fout, "time", typeof(pr_xi_mean&time),(/"time"/))

    filevardef(fout, "pr", typeof(pr_xi_mean),(/"time","lat","lon"/))



    filevarattdef(fout, "lat", pr_xi_mean&lat)
    filevarattdef(fout, "lon", pr_xi_mean&lon)
    filevarattdef(fout, "time", pr_xi_mean&time)
    
    filevarattdef(fout, "pr", pr_xi_mean)
    setfileoption(fout,"DefineMode",False)

    fout->lat=(/pr_xi_mean&lat/)
    fout->lon=(/pr_xi_mean&lon/)
    fout->time=(/pr_xi_mean&time/)
    fout->pr=(/pr_xi_mean/)

    delete(pr_xi_mean)
    delete(ntime)
    delete(nlon)
    delete(nlat)
    delete(time)
 end 




    
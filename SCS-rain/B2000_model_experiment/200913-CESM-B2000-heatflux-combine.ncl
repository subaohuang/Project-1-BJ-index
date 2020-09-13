
undef("combine_B2000_heatflux")
  ;;;************************************************************
  ;; 差值CESM混合坐标系，将其输出到标准气压坐标系上
  ;; vinth2p(datai, hbcofa, hbcofb, plevo, psfc, intyp, p0, ii, kxtrp)
  ;; CESM默认的P0 = 1000
  ;; 具体参数看ncl官网这个函数的介绍
  ;;functions cesm_hybird2pre(infilepath,outfilepath)
  ;; infilepath : absolute path of cesm data interpolated from hybrid to pressure level
  ;; OUTPUT variable : U,V,OMEGA,PSL,PRECL
  ;; linshuheng 2019.9.25
  ;;;*********************************************************
  function combine_B2000_heatflux(modelname,infilepath,outfilepath)
  begin
 ;;数据

    fils1 = systemfunc("ls " +infilepath+modelname+".cam.h0"+"*"+".nc")



    f1 = addfiles(fils1, "r")
 
    FLNS = f1[:]->FLNS                   ; read to memory [optional]
    FSNS = f1[:]->FSNS
    LHFLX = f1[:]->LHFLX
    SHFLX = f1[:]->SHFLX
  
    lat = f1[:]->lat
    lon = f1[:]->lon  
    time = f1[:]->time
    
    printVarSummary(FLNS)
    printVarSummary(FSNS)
    printVarSummary(LHFLX)
    printVarSummary(SHFLX)

    FLNS@missing_value = 9.96921e+36
    FSNS@missing_value = 9.96921e+36
    LHFLX@missing_value = 9.96921e+36
    SHFLX@missing_value = 9.96921e+36

    ; time = (/time - 1.0/) ;;;注意cesm得到的时候，time_bnd的第一位，所以读出来有问题
  
    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f "+outfilepath+"" )   ; remove any pre-existing file
    ncdf = addfile(outfilepath ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications
   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================
    fAtt               = True            ; assign file attributes
    fAtt@title         = "将通量数据合并"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes==
    ;=================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
    ncdf->FLNS  =  FLNS
    ncdf->FSNS  = FSNS
    ncdf->LHFLX = LHFLX
    ncdf->SHFLX = SHFLX
 
   printVarSummary(FLNS&lat)
    ncdf->lat   = FLNS&lat
    ncdf->lon = FLNS&lon 
    h = time - 1.0
    printVarSummary(h)  
    ; copy_VarCoords(time, h)
    ; printVarSummary(h)

    ncdf->time = h     ;;;;不要把time的信息赋值过去，不然time会修改，因为time这个变量默认会取coordinate中的情况
    ;;这种简单的写nc文件很容易出现这个问题，如果time这个变量没有坐标变量变量time，那么文件变量会用time这个变量的值作为
    ;;坐标变量
    print("****"+h(0))
    k= "处理完成"
    print("处理完成")
    return(k)
  end 

;;***************************************************************************
;;main program 
;;*************************************************************************

begin 
   ; inpath = "~/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611.nc"
   ; outpath = "~/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611_chazhi.nc"
   ; inpath = "/home/ys17-19/lsh/HEAT_MC_SON.cam.h0.25001-28912.nc"
   ; outpath = "/home/ys17-19/lsh/HEAT_MC_SON.cam.h0.25001-28912_chazhi.nc"

   m = combine_B2000_heatflux(modelname,inpath,outpath)
end 



    





    ;;
 
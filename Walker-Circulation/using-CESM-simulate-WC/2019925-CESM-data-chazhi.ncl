
undef("cesm_hybird2pre")
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
  function cesm_hybird2pre(infilepath,outfilepath)
  begin
 ;;数据
    
    fils1 = infilepath
    f1 = addfile(fils1, "r")
 
    hyam = f1->hyam                   ; read to memory [optional]
    hybm = f1->hybm
    U    = f1->U                     ; say "T" is (lat,lev,lon)
    V    = f1->V                     ; say "T" is (lat,lev,lon)
    PSL = f1->PSL
    OMEGA = f1->OMEGA
    PRECL = f1->PRECL

    lat = f1->lat
    lon = f1->lon  
    time = f1->time


    P0mb = 1000.0                        ; reference pressure [mb]
    pnew = (/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)

    lev = pnew
    
    ;;------------------将sigma混合坐标系进行差值---------------------------------
    Unew = vinth2p (U(time|:,lev|:,lat|:,lon|:),hyam,hybm,pnew,f1->PS,1,P0mb,1,True)


    Vnew = vinth2p (V(time|:,lev|:,lat|:,lon|:),hyam,hybm,pnew,f1->PS,1,P0mb,1,True)

    Wnew = vinth2p (OMEGA(time|:,lev|:,lat|:,lon|:),hyam,hybm,pnew,f1->PS,1,P0mb,1,True)

    Unew!1 = "lev"
    Unew&lev = lev 
   
    Vnew!1 = "lev"
    Vnew&lev = lev 
 
    Wnew!1 = "lev"
    Wnew&lev = lev 
   
    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f" +"" + outfilepath +"")   ; remove any pre-existing file
    ncdf = addfile(outfilepath ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications
   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================
    fAtt               = True            ; assign file attributes
    fAtt@title         = "差值CESM混合坐标系，将其输出到标准气压坐标系上"
    fAtt@source_file   =  "f19_f19_FAMIP.cam.h0.197901-200611.nc"
    fAtt@Conventions   = "CHAZHI"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes==
    ;=================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
    ncdf->U  =  Unew
    ncdf->V  = Vnew

    ncdf->OMEGA   =  Wnew
    ncdf->PRECL = PRECL
    ncdf->PSL = PSL

    ncdf->lat   =  lat
    ncdf->lon = lon 
    ncdf->lev = lev
    ncdf->time = time

    k= "插值处理完成"
    print(k)
    return(k)
  end 

;;***************************************************************************
;;main program 
;;*************************************************************************

begin 
   ; inpath = "~/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611.nc"
   ; outpath = "~/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611_chazhi.nc"
   m = cesm_hybird2pre(inpath,outpath)
end 



    





    ;;
 
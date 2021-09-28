
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
    vNames = getfilevarnames (f1) ; get names of all variables on file

    nNames = dimsizes (vNames)   ; number of variables on the file

    print (vNames)               ; print all variable names on file

    hyam = f1->hyam                   ; read to memory [optional]
    hybm = f1->hybm
    lat = f1->lat
    lon = f1->lon 
    time = f1->time
    

    P0mb = 1000.0                        ; reference pressure [mb]
    pnew = (/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)

    lev = pnew
    lev!0 = "lev"
    lev&lev = pnew
    lev&lev@units = "hPa"
    

    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f "+outfilepath+"" )   ; remove any pre-existing file
    ncdf = addfile(outfilepath ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications
   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================
    fAtt               = True            ; assign file attributes
    fAtt@title         = "差值CESM混合坐标系，将其输出到标准气压坐标系上"
    fAtt@Conventions   = "CHAZHI"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes==
    ;=================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
    ncdf->lat   =  lat
    ncdf->lon = lon 
    ncdf->lev = lev
    h = time - 1.0
    ncdf->time = h     ;;;;不要把time的信息赋值过去，不然time会修改，因为time这个变量默认会取coordinate中的情况
    ;;这种简单的写nc文件很容易出现这个问题，如果time这个变量没有坐标变量变量time，那么文件变量会用time这个变量的值作为

   do n=0,nNames-1              ; loop thru each variable

     dnames := getfilevardims(f1,vNames(n))
 
     vardimnum:=dimsizes(dnames)

    if vardimnum.ge.3 then 
        v := f1->$vNames(n)$        ; read the varible to memory
        if vardimnum.ge.4 then 
    ;;------------------将sigma混合坐标系进行差值---------------------------------
        vnew = vinth2p(v,hyam,hybm,pnew,f1->PS,1,P0mb,1,True)

        vnew!1 = "lev"
        vnew&lev = lev 
        copy_VarAtts(v, vnew)
        ncdf->$vNames(n)$   =  vnew
        delete(vnew)
      else  
        ncdf->$vNames(n)$   =  v
      end if
     delete(v)
    end if
   end do
    ;;坐标变量
    k= "插值处理完成"
    print("插值处理完成")
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
   m = cesm_hybird2pre(inpath,outpath)
end 



    





    ;;
 
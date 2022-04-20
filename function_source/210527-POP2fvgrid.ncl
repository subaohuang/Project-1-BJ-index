load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 

undef("cesm_gx1v6tofvgrid")
  ;;;************************************************************
  ;; 差值CESM POP文件的水平精度gx1v6 to fv grid
  ;; PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
  ;; 具体参数看ncl官网这个函数的介绍PopLatLon
  ;;functions cesm_hybird2pre(infilepath,outfilepath)
  ;; infilepath : absolute path of cesm data interpolated from hybrid to pressure level
  ;; OUTPUT variable : U,V,OMEGA,PSL,PRECL
  ;; linshuheng 2021.5.27
  ;;;*********************************************************
  function cesm_gx1v6tofvgrid(infilepath,outfilepath)
  begin
 ;;数据
    fils1 = infilepath
    f1 = addfile(fils1, "r")
    vNames = getfilevarnames (f1) ; get names of all variables on file

    nNames = dimsizes (vNames)   ; number of variables on the file

    print (vNames)               ; print all variable names on file
  
    time = f1->time
    
    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f "+outfilepath+"" )   ; remove any pre-existing file
    ncdf = addfile(outfilepath ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications
   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================
    fAtt               = True            ; assign file attributes
    fAtt@title         = "POP模式插值"
    fAtt@Conventions   = "CHAZHI"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes==
    ;=================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
    ncdf->time = time     ;;;;不要把time的信息赋值过去，不然time会修改，因为time这个变量默认会取coordinate中的情况
    if any(vNames.eq."z_w_top") then 
      z_w_top = f1->z_w_top
      ; z_w_top = (/z_w_top/100.0/)
      ncdf->z_w_top =z_w_top
    end if
    if any(vNames.eq."z_t") then 
      z_t = f1->z_t  
      ; z_t = (/z_t/100.0/)
      ncdf->z_t =z_t
    end if
    ;;这种简单的写nc文件很容易出现这个问题，如果time这个变量没有坐标变量变量time，那么文件变量会用time这个变量的值作为

   do n=0,nNames-1              ; loop thru each variable

     dnames := getfilevardims(f1,vNames(n))
     vardimnum:=dimsizes(dnames)

    if vardimnum.ge.3 then 
      v := f1->$vNames(n)$        ; read the varible to memory
      v_chazhi = PopLatLon(v,"gx1v6","fv1.9x2.5","bilin","da","090206")
      printVarSummary(v_chazhi)
      ncdf->$vNames(n)$   =  v_chazhi
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
   m = cesm_gx1v6tofvgrid(infilepath,outfilepath)
end 



    





    ;;
 
;;********************************************************
;;检索SODA数据的MLD位置
;;                                         lnshuheng 2020.04.16
;;
;**********************************************************

;;**************************************************************
;;;main program
;;***************************************************************
begin


;;;region
   lat1 = -10
   lat2 = 20
   lon1 = 100
   lon2 = 170
;;;




    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_prho_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_potensial_density_mld_0.8_2015.nc", "r") 

    indyear = (2015-1980+1)*12-1
  
    pd_r = f5->prho(0:indyear,0:20,{lat1:lat2},{lon1:lon2})


    pd = tofloat(pd_r)
    copy_VarCoords(pd_r, pd)
    delete(pd_r)
    printVarSummary(pd)
    

    pd_mld_r = dble2flt(f6->phro_mld)
    
    pd_mld = pd_mld_r(time|:,lat|:,lon|:)
    pd_mld = where(pd_mld.eq."-nan", 10^9, pd_mld)
    pd_mld@_FillValue = 10^9



    pd_mld_scs = pd_mld(:,{lat1:lat2},{lon1:lon2})
   
   

;;;差值

 ;;;将temp进行垂直方向差值，为了得到更加精确的结果
    print("开始差值")
    xi  = pd&depth 
    xo = ispan(5,330,2)
  



    pd_chazhi = linint1_n_Wrap(xi,pd(:,:,:,:),False,xo,0,1)
    pd_chazhi!1 = "depth"
    pd_chazhi&depth = xo
    printVarSummary(pd_chazhi)
    delete(pd)

    
    dim = dimsizes(pd_chazhi)
    printVarSummary(pd_chazhi)
    pd_mld_4d = conform_dims(dim,pd_mld_scs, (/0,2,3/))
    copy_VarMeta(pd_chazhi, pd_mld_4d)
    printVarSummary(pd_mld_4d)
    pd_diff = pd_chazhi
    pd_diff =  pd_mld_4d- pd_chazhi
    printVarSummary(pd_diff)
    depth  = pd_diff&depth
    depth_4d = pd_diff
    depth_4d = conform_dims(dim, depth, 1)
    mlt_4d = depth_4d
    mlt_4d = where(pd_diff.ge.0, mlt_4d, mlt_4d@_FillValue)

    mlt_t = mlt_4d(:,1,:,:)
    mlt_t = dim_max_n(mlt_4d, 1) + 1

    ; f7=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_mlp_mn_1980-2015_chazhi.nc", "r") 
    ; mlp = dble2flt(f7->mlp(0:indyear,{lat1:lat2},{lon1:lon2}))
    ; printVarSummary(mlp)

    ; mld_diff = mlp - mlt_t
    ; print(mld_diff(1,:,:))

 print("写入nc文件")   
 ;;;;写nc文件,简单的方法  
   system("/bin/rm -f SODA3.4.2_mld_density_criterion_0.8_2015.nc")   ; remove any pre-existing file
    
  ncdf = addfile("SODA3.4.2_mld_density_criterion_0.8_2015.nc" ,"c")  ; open output netCDF file

  ;===================================================================
  ; create global attributes of the file (optional)
  ;===================================================================
     fAtt               = True            ; assign file attributes
     fAtt@title         = "SODA3.4.2 mld by using density criterion 1 degree ,198001-201312"
     fAtt@Conventions   = "None"
     fAtt@creation_date = systemfunc ("date")
     fileattdef( ncdf, fAtt )            ; copy file attributes

  ;===================================================================
  ; make time an UNLIMITED dimension; recommended  for most applications
  ;===================================================================
     filedimdef(ncdf,"time",-1,True) 

  ;===================================================================
  ; output variables directly; NCL will call appropriate functions
  ; to write the meta data associated with each variable
  ;===================================================================
     ncdf->mld  = mlt_t                       ; 4D      

  ;;;;
  exit()



  
   end 

   





      


        
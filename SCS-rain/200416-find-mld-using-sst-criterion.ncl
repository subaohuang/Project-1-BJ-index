;;********************************************************
;;检索SODA数据的MLD位置
;;                                         lnshuheng 2020.04.16
;;
;**********************************************************

;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
;;;region
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    sst_scs = sst(:,{lat1:lat2},{lon1:lon2})


;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    ;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
  
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_mlt_mn_1980-2015_chazhi.nc", "r") 

    indyear = (2013-1980+1)*12-1
  
    T_r = f5->temp(0:indyear,0,:,:)

  
    mlt_r = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})

    T = tofloat(T_r)

    mlt = tofloat(mlt_r)

   copy_VarCoords(T_r, T)
   copy_VarCoords(mlt_r, mlt)

       ;;;;写nc文件,简单的方法  
           system("/bin/rm -f SODA3.4.2_temp_5m.nc")   ; remove any pre-existing file
            
          ncdf = addfile("SODA3.4.2_temp_5m.nc" ,"c")  ; open output netCDF file

          ;===================================================================
          ; create global attributes of the file (optional)
          ;===================================================================
             fAtt               = True            ; assign file attributes
             fAtt@title         = "SODA3.4.2_temp at 5m,198001-201312"
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
             ncdf->temp  = T                       ; 4D      
        
          ;;;;
          exit()






   delete(mlt_r)
   delete(T_r)
  
 
;;;

;;;差值

 ;;;将temp进行垂直方向差值，为了得到更加精确的结果
    print("开始差值")
    xi  = T&depth 
    xo = ispan(5,330,2)

    T_chazhi = linint1_n_Wrap(xi,T(:,:,:,:),False,xo,0,1)
    T_chazhi!1 = "depth"
    T_chazhi&depth = xo
    printVarSummary(T_chazhi)
    delete(T)

    mlt = T_chazhi(:,1,:,:)
     
    dim = dimsizes(T_chazhi)
    printVarSummary(sst_scs)
    printVarSummary(T_chazhi)
    sst_4d = conform_dims(dim,sst_scs, (/0,2,3/))

    printVarSummary(sst_4d)
    T_diff = T_chazhi
    T_diff =  sst_4d - T_chazhi
    depth  = T_diff&depth
    depth_4d = T_diff
    depth_4d = conform_dims(dim, depth, 1)
    mlt_4d = depth_4d
    mlt_4d = where(T_diff.gt.0.5, mlt_4d, mlt_4d@_FillValue)
    mlt_t = dim_min_n(mlt_4d, 1)
  

    mlt_diff = mlt - mlt_t
    print(mlt_diff)
    exit()

   end 

   





      


        
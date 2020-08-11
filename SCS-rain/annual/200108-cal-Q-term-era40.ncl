;;********************************************************
;;计算ERA40中的水汽通量算度，并将其输出
;;数据的精度是1.5×1.5°
;; ERA-40:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2020.01.08
;**********************************************************
begin
  ;;;read era data
		path_wind = "/home/yangsong3/data-observation/ERA-40/monthly/"
		var_name = (/"uwnd","vwnd","shum"/)
		filein = path_wind + var_name +".mon.mean.deg1.5.nc"
        
        ; ;;;region 
        ;   lat1 = -10
        ;   lat2 = 50
        ;   lon1 = 60
        ;   lon2 = 160
        ; ;;;
        
        ;;;region 
          lat1 = -90
          lat2 = 90
          lon1 = -180
          lon2 = 180
        ;;;



        ;;;time span
        startmon =195701
        endmon = 200212

        startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
        endyear_r = str_split_by_length(endmon, 4 )

        startyear = stringtoint(startyear_r(0))
        endyear = stringtoint(endyear_r(0))

        yearnum = endyear - startyear +1
        ;;


		f1 = addfile(filein(0),"r")

		date=cd_calendar(f1->time, 1)
		timeind=ind(date.ge.startmon.and.date.le.endmon)

		u = short2flt(f1->u(timeind,:,{lat1:lat2},{lon1:lon2}))
    
    u = lonFlip(u)
    printVarSummary(u)

		f2 = addfile(filein(1),"r")
		v = short2flt(f2->v(timeind,:,{lat1:lat2},{lon1:lon2}))
    v = lonFlip(v)
		
		f3 = addfile(filein(2),"r")
		date :=cd_calendar(f3->time, 1)
		timeind :=ind(date.ge.startmon.and.date.le.endmon)
		q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
    q = lonFlip(q)
   ;;;
   print("=====finish read data==========")
;==============================================================
; Mission 1. calculate Q
; ==============================================================
    era_qu  = q
    era_qv  = q
    era_qu  = u*q
    printVarSummary(era_qu)
    era_qv  = v*q
    era_wgt_qu = era_qu(:,0,:,:)
    era_wgt_qv = era_qv(:,0,:,:)

    era_dQdx_flt  = era_wgt_qu
    era_dQdx   = flt2dble(era_dQdx_flt)
    era_dQdy_flt  = era_wgt_qv
    era_dQdy   = flt2dble(era_dQdy_flt)


   ;;;------Cal DivQ-------
    ;-----------Follow is the vertical integral-----------
    pver=era_qu&level
    pver=pver*100
    pver@units="Pa"

    dp=dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(era_qu),dp,1)          ; confirm dp to caculate sum(var*dp)

    era_wgt_qu=(/dim_sum_n(dp_cfm*era_qu,1)/)               ; sum(var*dp)
    era_wgt_qv=(/dim_sum_n(dp_cfm*era_qv,1)/)               ; sum(var*dp)

    era_wgt_qu = (/era_wgt_qu/9.8/)
    era_wgt_qv = (/era_wgt_qv/9.8/)

    g_lon   =   era_wgt_qu&longitude
    g_lat   =   era_wgt_qu&latitude
    g_nlat  =   dimsizes(g_lat)

    DivQ = uv2dv_cfd(era_wgt_qu,era_wgt_qv,g_lat,g_lon,2)
    ; DivQ = uv2dvF_Wrap(era_wgt_qu(:,::-1,:),era_wgt_qv(:,::-1,:)) 
    DivQ   =   (/-DivQ*86400/) 
    copy_VarCoords(era_wgt_qu, DivQ)
    printVarSummary(DivQ)
    ; print(DivQ(0,:,:))
    print("=====finish DivQ calculate==========")

    dim = dimsizes(DivQ)
    nlat = dim(1)
    nlon = dim(2)


    ;;;写成nc 文件输出

    outpath = "/home/ys17-19/lsh/Project/SCS-rain/era-40-DivQ.nc"
       ;;-----------------------写入nc文件---------------------------
    system("rm " + outpath)   ; remove any pre-existing file
    fout = addfile(outpath ,"c")  ; open output netCDF file
    setfileoption(fout,"DefineMode",True)

   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================

    fAtt               = True            ; assign file attributes
    fAtt@title         = "1957-2002，精度为1.5乘1.5"    
    fAtt@creation_date = "2017.12.21"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/-1,nlat,nlon/)
    dimunlim=(/True,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    filevardef(fout, "lat", typeof(DivQ&latitude),(/"lat"/))
    filevardef(fout, "lon", typeof(DivQ&longitude),(/"lon"/))
    filevardef(fout, "time", typeof(DivQ&time),(/"time"/))
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "DivQ",typeof(DivQ),(/"time","lat","lon"/))

      ;;属性
      filevarattdef(fout, "lon", DivQ&longitude)
      filevarattdef(fout, "lat", DivQ&latitude)
      filevarattdef(fout, "time", DivQ&time)
      filevarattdef(fout, "DivQ", DivQ)

      ;;

      fout->lon=(/DivQ&longitude/)
      fout->lat=(/DivQ&latitude/)
      fout->time=(/DivQ&time/)
      fout->DivQ=(/DivQ/)
       print("=====finish NC write==========")
    ;;;

end


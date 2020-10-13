  begin
  diri="/media/sf_share/ncep/20CR/"
    fils1=systemfunc("ls "+diri+"uwnd.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"vwnd.mon.mean.nc")
    fils4=systemfunc("ls "+diri+"omega.mon.mean.nc")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f4=addfile(fils4, "r")
    a=6371000
    g=9.8
    pi=get_pi("double")
    dfi=(pi*10)/180
    uwnd=f1->uwnd(1295:1703,:,::-1,:)
    vwnd=f2->vwnd(1295:1703,:,::-1,:)
    level=f1->level(:)
    w=f4->omega(1295:1703,:,48:42,50:145)
    w_new=w(time|:,lat|:,level|:,lon|:) 
    w_avg=dim_avg_n_Wrap(w_new, (/0,1/))
   
    vrdvx=uv2vrdvF(uwnd,vwnd)

    uvd=dv2uvF(vrdvx(1,:,:,:,:))
    delete(vrdvx)
    ud=uvd(0,:,:,:,:)
    delete(uvd)
    ud!0="time"
    ud!1="level"
    ud!2="lat"
    ud!3="lon"
    ud_new=ud(time|:,lat|48:42,lon|50:145,level|:)  

    ud_new&lon=w&lon
    ud_new&lat=fspan(-6,6,dimsizes(w&lat))
    ud_new&time=ispan(1,dimsizes(w&time),1)
    ud_new&level=level
    
    ud_avg=dim_avg_n_Wrap(ud_new,(/0,1/))
    ud_avg_new=ud_avg(level|0:dimsizes(w&level)-1,lon|:)  
    printVarSummary(ud_avg_new)
    w_avg_new=w_avg*(-50)
    copy_VarCoords(w_avg, w_avg_new)
    printVarSummary(w_avg_new)
   
    pbot=w_avg_new&level
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    do i=0,(dimsizes(w_avg_new&level)-1)
    vint(i,:,:,:)=vibeta(level, ud_new, 1, pres_new,pbot(i),0)
    end do
    fi=((a*dfi)/g)*vint*100
    delete(vint)

    fi_mean=dim_avg_n_Wrap(fi, (/1,2/))
    delete(fi)

    fi_mean!0="level"
    fi_mean!1="lon"
    fi_mean&level=pbot
    fi_mean&level@units="hpa"
    fi_mean&lon=ud_new&lon
  
    nlev=dimsizes(fi_mean&level)
    nlon=dimsizes(fi_mean&lon)
     printVarSummary(fi_mean) 
  
     ;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro ="/media/sf_share/ncep/20CR/"
    filo="msf_20cr.nc"
    system("/bin/rm -f "+diro + filo+"")
    fout=addfile(diro + filo,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量，取6到-6平均的结果"    
    fAtt@creation_date = "2017.7.31"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon"/)
    dimsize=(/nlev,nlon/)
    dimunlim=(/False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "msf", typeof(fi_mean),getvardims(fi_mean))
    filevardef(fout, "ud", typeof(ud_avg_new),getvardims(ud_avg_new))
    filevardef(fout, "w", typeof(w_avg_new),getvardims(w_avg_new))

    filevarattdef(fout, "msf", fi_mean)
    filevarattdef(fout, "ud", ud_avg_new)
    filevarattdef(fout, "w", w_avg_new)

    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_mean&level/)
    fout->lon=(/fi_mean&lon/)
    fout->msf=(/fi_mean/)
    fout->ud=(/ud_avg_new/)
    fout->w=(/w_avg_new/)

end 




;     wks=gsn_open_wks("png", ""+diri+"mass_stream_20cr")

;     res = True
;     res@gsnDraw=False
;     res@gsnFrame=False
;     res@gsnLeftString=""
;     res@gsnRightString=""

;     res@cnFillPalette  = "BlueDarkRed18"
;     res@cnLineColors="Transparent"
;    ; res@cnLevelSelectionMode ="ManualLevels"
;    ; res@cnMinLevelValF = -1
;    ; res@cnMaxLevelValF = 1
;    ; res@cnLevelSpacingF = 0.1
;     res@cnLineLabelsOn=False
;     res@cnLineLabelFontHeightF=0.01
;     res@cnFillOn = True 
;     res@cnLinesOn = False

;     res@tmXTOn=False
;     res@tmYROn=False

;    res@vcRefMagnitudeF = 2.0                ; define vector ref mag
;    res@vcRefLengthF    = 0.045              ; define length of vec ref
;    ; res@vcGlyphStyle    = "CurlyVector"      ; turn on curly vectors
;    res@vcMinDistanceF  = 0.01               ; thin out vectors
;    res@vcMapDirection  = False


;    plot=gsn_csm_pres_hgt_vector(wks,fi_mean,ud_avg_new,w_avg_new,res)
; draw(plot)
; frame(wks)



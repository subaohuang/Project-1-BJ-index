 function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin
    ; data_new=doubletofloat(dim_standardize_n_Wrap(data(::-1,:,:),1,(/0,1/)))       ;标准化
    ; print(data_new)
    ; copy_VarMeta(data(::-1,:,:), data_new)
    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)
    ;-----------------------------------------------
    rc = regCoef_n(data&time,data,0,2)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data&level
    rc&lon=data&lon
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
     ; print(rc)
return(rc)
end
;-------------------------------------
function normoalized(data,mean)
begin

msf_enino_tr=data/mean
copy_VarCoords(data, msf_enino_tr)
msf_enino_tr!0="lev"
msf_enino_tr&lev@units = "hpa"
return msf_enino_tr
end 
begin 

diri="/mnt/f/share/era/data/"
fils0=systemfunc("ls "+diri+"era_nino3_4.ascii")
nino3_4=asciiread(fils0, 456, "float")
nino3_4_de=dtrend_n(nino3_4(0:407), True,0)
printVarSummary(nino3_4_de)
nino3_4_de!0="time"
nino3_4_de&time=ispan(1,408,1)

nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)

 enino=ind(nino3_4_de.gt.1.0)
 lanila=ind(nino3_4_de.lt.-1.0)
 Neutral=ind(nino3_4_de.lt.1.0.and.nino3_4_de.gt.-1.0)

 ; str_enino=ind(nino3_4_nor.gt.1.0)
 ; str_lanila=ind(nino3_4_nor.lt.-1.0)
 ; weak_enino=ind(nino3_4_nor.gt.0.5.and.nino3_4_nor.lt.1.0)
 ; weak_lanila=ind(nino3_4_nor.gt.-1.0.and.nino3_4_nor.lt.-0.5)

 str_enino=ind(nino3_4_de.gt.1.0)
 str_lanila=ind(nino3_4_de.lt.-1.0)
 weak_enino=ind(nino3_4_de.gt.0.5.and.nino3_4_de.lt.1.0)
 weak_lanila=ind(nino3_4_de.gt.-1.0.and.nino3_4_de.lt.-0.5)



  diri3="/mnt/f/share/era/data/"
  fils3=systemfunc("ls "+diri3+"msf_time_de_re.nc")
  f3=addfile(fils3, "r")
  fils4=systemfunc("ls "+diri3+"msf_time_era.nc")
  f4=addfile(fils4, "r")
  msf_time=f4->msf_time(:,:,0:407)
  msf_time&level@units="hpa"
  msf_clm=dim_avg_n_Wrap(msf_time, 2)
  msf_scale_clm=msf_clm/10^9
  copy_VarMeta(msf_clm, msf_scale_clm)
   msf_time_de_re=f3->msf_time_de_re
   msf_mon_clm=clmMonLLT(msf_time)/10^9
   printVarSummary(msf_mon_clm)
   printVarSummary(msf_time_de_re)
   msf_time_re=new(dimsizes(msf_time), typeof(msf_time))
   m=dimsizes(msf_time&time)-1

    do i=0,11
      do j= i,m,12
    msf_time_re(:,:,j)=(msf_mon_clm(:,:,i)+msf_time_de_re(:,:,j))*10^9
     end do
   end do

   printVarSummary(msf_time_re)
   copy_VarCoords(msf_time_de_re, msf_time_re)

   zo=(msf_time- msf_time_de_re)/10^9
   print("z0="+zo(:,:,2))

;---------------------------------------------------------------------
;写文件，将得到的文件用location插值，用在locatio插值函数。进去确定上升与下沉支位置。
;--------------------------------------------------------------------
;随时间变化的量
    diro=diri
    filo1="re_msf_time_era.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf重构的去趋势的异常场"    
    fAtt@creation_date = "2017.9.5"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/dimsizes(msf_time_re)/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(msf_time_re&level), (/"level"/))
    filevardef(fout, "lon", typeof(msf_time_re&lon),(/"lon"/))
    filevardef(fout, "time", typeof(msf_time_re&time),(/"time"/))
    filevardef(fout, "msf_time_re", typeof(msf_time_re),(/"level","lon","time"/))

    filevarattdef(fout, "msf_time_re", msf_time_re)
  
    setfileoption(fout,"DefineMode",False)

    fout->level=(/msf_time_re&level/)
    fout->lon=(/msf_time_re&lon/)
    fout->time=(/msf_time_re&time/)
    fout->msf_time_re=(/msf_time_re/)
    

  end 
;************************************************************
 ;计算只有时间维的月气候态
;**********************************************************
function clm (x[*]:numeric)
begin

   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check
;
; Compute all 12 monthly averages first. 
;
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,aveMonth)
; copy spatial (eg, lat/lon) coordinate variables
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   return (aveMonth)
end

function  msfzeroline_location(msf_time)
 begin 
 pout=ispan(100,1000,10)
 msf_int=int2p_n_Wrap(msf_time&level, msf_time,pout(::-1), 1, 0)

 xout=fspan(100, 290, 1901)
 msf_int_f=linint1_n_Wrap(msf_int&lon, msf_int, False,xout,0,1)
 delete(msf_int)

 msf_scale=msf_int_f/10^9
 copy_VarMeta(msf_int_f, msf_scale)
 delete(msf_int_f)
 ; levind1=ind(msf_scale&LEVEL.eq.600)
 ; levind2=ind(msf_scale&LEVEL.eq.400)
  levind1=ind(msf_scale&LEVEL.eq.1000)
  levind2=ind(msf_scale&LEVEL.eq.150)

 LON=msf_scale&LON
 locationxi=new((/(levind2-levind1)+1,dimsizes(msf_scale&time)/),float)
 locationdong=new((/(levind2-levind1)+1,dimsizes(msf_scale&time)/),float)

  do t=0,dimsizes(msf_scale&time)-1,1
   do l=levind1,levind2
    do j=0,ind(msf_scale&LON.eq.220)-1,1    

      if  (msf_scale(l,j+1,t).gt.0).and.((msf_scale(l,j,t).lt.0).and.(msf_scale(l,j+400,t).gt.0))  then 
          locationxi(l-levind1,t)=(LON(j+1)+LON(j))/2
      end if    
      if ismissing(locationxi(l-levind1,t)).ne.True  then 
         break
     end if

    end do
   end do
end do

  do t=0,dimsizes(msf_scale&time)-1,1
   do l=levind1,levind2
    do j=ind(msf_scale&LON.eq.180),dimsizes(msf_scale&LON)-2,1   

     if ((msf_scale(l,j+1,t).lt.0).and.(msf_scale(l,j,t).gt.0)) then 
             locationdong(l-levind1,t)=(LON(j+1)+LON(j))/2
     end if
  
      if ismissing(locationdong(l-levind1,t)).ne.True  then 
         break
     end if

     end do 
   end do
end do
 
 locationxi!0="level"
 locationxi!1="time"
 locationxi&level=msf_scale&LEVEL(levind1:levind2)
 locationxi&time=msf_scale&time

 locationdong!0="level"
 locationdong!1="time"
 locationdong&level=msf_scale&LEVEL(levind1:levind2)
 locationdong&time=msf_scale&time


 location=new((/2,dimsizes(locationxi&level),dimsizes(locationxi&time)/),float)
 ; locationxi_a!0="time"
 ; locationdong_a!0="time"
 ; locationxi_a&time=msf_scale&time
 ; locationdong_a&time=msf_scale&time
 location(0,:,:)=locationxi
 location(1,:,:)=locationdong
return (location)
end 

begin 

 diri="/mnt/f/share/era/data/"
 ; fils1=systemfunc("ls "+diri+"msf_time_era_.nc")
 fils1=systemfunc("ls "+diri+"re_msf_time_era.nc")  ;重构的场的位置。
 f1=addfile(fils1, "r")

 ; msf_time=f1->msf_time
 msf_time=f1->msf_time_re
 level=f1->level(:)
 print(level)
 lon=f1->lon(:)
 time=f1->time 
 levind1=ind(msf_time&level.eq.1000)     ;所求0线位置的level
 levind2=ind(msf_time&level.eq.150)      ;所求0线位置的level
 ntime=dimsizes(time)
 msf_time&level@units="hpa"
 
 location=msfzeroline_location(msf_time)
 locationxi=location(0,:,:)
 locationdong=location(1,:,:)
 printVarSummary(locationxi)

 locationxi_a=dim_avg_n_Wrap(locationxi(ind(locationxi&level.eq.600):ind(locationxi&level.eq.400),:),0)
 locationdong_a=dim_avg_n_Wrap(locationdong(ind(locationdong&level.eq.600):ind(locationdong&level.eq.400),:),0)
 
 ; locationxi_clm=clm(locationxi_a(0:455))
  locationxi_clm=clm(locationxi_a(0:407))

 ; h=ind(ismissing(locationxi_a))
 ; print("h="+h)
 ; m=mod(h+1,12)
 ; locationxi_a =where(ismissing(locationxi_a), locationxi_clm(m-1), locationxi_a)


 lolevel=locationxi&level
 nlevel=dimsizes(locationxi&level)
 

;---------------------------------------------------------------------
;写文件
;--------------------------------------------------------------------
    filo="msf_era_level_location1979-2017_reconsitution.nc"
    system("/bin/rm -f "+diri + filo+"")
    fout=addfile(diri + filo,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf的0线位置，插值后的结果,并有时间项"    
    fAtt@creation_date = "2017.9.5"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","time"/)
    dimsize=(/nlevel,ntime/)
    dimunlim=(/False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level", "float",(/"level"/))
    filevardef(fout, "time", typeof(time),getvardims(time))
    filevardef(fout, "locationxi", typeof(locationxi),getvardims(locationxi))
    filevardef(fout, "locationdong", typeof(locationdong),getvardims(locationdong))
    filevardef(fout, "locationxi_a", typeof(locationxi_a),getvardims(locationxi_a))
    filevardef(fout, "locationdong_a", typeof(locationdong_a),getvardims(locationdong_a))

    filevarattdef(fout, "locationxi", locationxi)
    filevarattdef(fout, "locationdong", locationdong)
    filevarattdef(fout, "locationxi_a", locationxi_a)
    filevarattdef(fout, "locationdong_a", locationdong_a)
   

    setfileoption(fout,"DefineMode",False)

    fout->locationxi=(/locationxi/)
    fout->locationdong=(/locationdong/)
    fout->locationxi_a=(/locationxi_a/)
    fout->locationdong_a=(/locationdong_a/)
    fout->level=(/lolevel/)
    fout->time=(/time/)

end 










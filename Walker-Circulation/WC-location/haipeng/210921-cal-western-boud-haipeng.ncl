;;;
undef("level_mean")
function level_mean(msf_time)   ;;;;select 850hpa-200hpa 
begin
; msf_time=(/msf_time/10^9/)
plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
; dp=(/50/)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end
;;;
undef("search_location")
function  search_location(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:290})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+8).gt.0))  then 
         a=(msf_mean_annual_xi(t,i+1)-msf_mean_annual_xi(t,i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi(t)=tofloat(msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(t,i)/a)
        ; locationxi(t)=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if ismissing(locationxi(t)).ne.True  then 
        break
        end if
     end do
end do

;;;judge locationdong
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0))  then 
         a=(msf_mean_annual_dong(t,i+1)-msf_mean_annual_dong(t,i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong(t)=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a)
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
   if ismissing(locationdong(t)).eq.True  then 
   ; co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
   co_in=closest_val(0,msf_mean_annual_dong(t,:))
   locationdong(t)=msf_mean_annual_dong&lon(co_in)
   delete(co_in)
  end if

end do
locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;********************************************************
begin
    diri="/home/ys17-19/lsh/Project/Walker-Circulation/WC-location/haipeng/"

    fils1=systemfunc("ls "+diri+"msf_haipeng_1979-2017.nc")

    ;;;********************************************* 
    f1=addfile(fils1, "r")
    time = f1->time 
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    ; msf_time=f1->msf_time({500},:,timeind)
    msf_time= dble2flt(f1->msf_time)
    levind=ind(msf_time&lev.eq.500)
    msf_time=(/msf_time/10^9/)

    msf_mean = dim_avg_n_Wrap(msf_time(:,{400:600},:), 1)
    printVarSummary(msf_mean)
    delete(msf_time)

    location=search_location(msf_mean)
    locationxi= dble2flt(location[0])
    ; locationdong=dble2flt(location[1])
    ; locationdong_rea(i,:)=locationdon
    ;;;;
    delete(msf_mean)

    locationxi_rea = locationxi
    locationxi_rea!0="time"
    locationxi_rea&time:=time

    printVarSummary(locationxi_rea)

    ntime=dimsizes(locationxi_rea&time)
    print(locationxi_rea)
    

;***************************************************
;写入nc文件
;*************************************************
    diro=diri
    filo1="location_1980-2017.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "WalKer Circulation western Boundary (longitude) which is defined by the location where 400-600 hPa mean Mass stream fuction equal zero"    
    fAtt@creation_date = "2021.09.21"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time"/)
    dimsize=(/ntime/)
    dimunlim=(/False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time", typeof(locationxi_rea&time),(/"time"/))
    filevardef(fout, "locationxi", typeof(locationxi_rea),(/"time"/))
    filevarattdef(fout, "locationxi", locationxi_rea)
    ; filevarattdef(fout, "locationdong", locationdong_rea)
    setfileoption(fout,"DefineMode",False)

    fout->time=(/locationxi_rea&time/)
    fout->locationxi=(/locationxi_rea/)
end
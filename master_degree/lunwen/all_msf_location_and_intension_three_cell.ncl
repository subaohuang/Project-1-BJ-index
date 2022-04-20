;********************************************************
;計算三个环流的强度变化
;*******************************************************
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
function  ploy(locationxi_rea)
begin
dim=dimsizes(locationxi_rea)
dim1=dim(1)
mnmx      = new ((/3,dim(1)/), typeof(locationxi_rea))
mnmx(0,:) = dim_avg_n_Wrap(locationxi_rea, 0)
mnmx(1,0:dim1-1) = dim_min_n_Wrap(locationxi_rea, 0)
mnmx(2,0:dim1-1) = dim_max_n_Wrap(locationxi_rea, 0)
return mnmx
end 
;;;
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
;;; 
undef("search_location_alt") ;;;altlantic
function  search_location_alt(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{270:330})
msf_mean_annual_dong=msf_mean_annual(:,{320:370})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-3
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+2).gt.0))  then 
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
   if ismissing(locationxi(t)).eq.True  then 
   co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_xi(t,:))
   test=msf_mean_annual_xi(t,{290:330})
   co_in=closest_val(0,msf_mean_annual_xi(t,{290:330}))
   locationxi(t)=test&lon(co_in)
   delete(co_in)
  end if
end do
; locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)
location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;;;************************************************************
;;; main code
;;************************************************************
begin
    rea= (/"20cr","era","jra55","MERRA","ncep1","ncep2"/)

    id_inten=new((/6,34/),float)
    pa_inten=new((/6,34/),float)
    al_inten=new((/6,34/),float)
    
    width_pa=new((/6,34/),float)
    width_id=new((/6,34/),float)
    width_alt=new((/6,34/),float)

    do i=0,dimsizes(rea)-1
	diri="/mnt/f/share/rea_msf/"
	LeftString=rea(i)
	f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
	date=cd_calendar(f1->time, 1)
	timeind=ind(date.ge.197901.and.date.le.201212)
	
    year1=1979
    year2=2012
    yearend=201212
	msf_time=f1->msf_time(:,:,timeind)
	
    delete(f1)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    delete(timeind)
    printVarSummary(msf_time)
    msf_time=(/msf_time/10^9/)

    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)

    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")

    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    ; msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    ; msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    delete(msf_mean)
    ;;;
    inten_id= dim_avg_n_Wrap(level_mean(msf_djf(:,:,{60:100})),1)
    inten_pa=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{150:210})),1)
    inten_al=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{310:350})),1)
   

    inten_id=(/dim_rmvmean(inten_id)*-1.0/)
    inten_pa=(/dim_rmvmean(inten_pa)/)
    inten_al=(/dim_rmvmean(inten_al)/)
    ;;;
    
    id_inten(i,:)=tofloat(inten_id)
    pa_inten(i,:)=tofloat(inten_pa)
    al_inten(i,:)=tofloat(inten_al)

    id_inten(i,:)=(/runave(id_inten(i,:), 3, 0)/)
    pa_inten(i,:)=(/runave(pa_inten(i,:), 3, 0)/)
    al_inten(i,:)=(/runave(al_inten(i,:), 3, 0)/)
    delete([/inten_id,inten_pa,inten_al/])

    ;;;;宽度的变化

    msf_cold_mean=msf_djf(:,{500},:)
     delete(msf_djf)
    location=search_location(msf_cold_mean)


    ;;
    pa_w=location[1] - location[0]
    pa_w=(/dim_rmvmean(pa_w)/)
    id_w=(/dim_rmvmean(location[0])/)

    printVarSummary(msf_cold_mean)
    msf_cold_mean_pi=lonPivot(msf_cold_mean, 30)
    delete(msf_cold_mean)
    printVarSummary(msf_cold_mean_pi)

    location_al=search_location_alt(msf_cold_mean_pi)
    delete(msf_cold_mean_pi)

    alt_w=location_al[1] - location_al[0]
    alt_w=(/dim_rmvmean(alt_w)/)


    width_pa(i,:)=tofloat(pa_w)
    width_id(i,:)=tofloat(id_w)
    width_alt(i,:)=tofloat(alt_w)

    width_id(i,:)=(/runave(width_id(i,:), 3, 0)/)
    width_pa(i,:)=(/runave(width_pa(i,:), 3, 0)/)
    width_alt(i,:)=(/runave(width_alt(i,:), 3, 0)/)
    delete([/pa_w,id_w,alt_w/])
 
  end do

    id=ploy(id_inten)
    pa=ploy(pa_inten)    
    al=ploy(al_inten)

    id_width=ploy(width_id)
    pa_width=ploy(width_pa)    
    al_width=ploy(width_alt)
    

    id_inten_avg=id(0,:)
    pa_inten_avg=pa(0,:)
    al_inten_avg=al(0,:)

    id_width_avg=id_width(0,:)
    pa_width_avg=pa_width(0,:)
    al_width_avg=al_width(0,:)

    ;;;计算趋势
    year=ispan(1, 34, 1)
    rc1=regCoef(ispan(1,34, 1),id_inten_avg)
    rc2=regCoef(ispan(1,34, 1),pa_inten_avg)
    rc3=regCoef(ispan(1,34, 1),al_inten_avg)
    
    paindex_rc =rc2*year+rc2@yintercept
    idindex_rc =rc1*year+rc1@yintercept
    alindex_rc =rc3*year+rc3@yintercept
    
    h1=new((/8,34/),float)
    h2=new((/8,34/),float)
    h3=new((/8,34/),float)
    h4=new((/7,34/),float)
    h5=new((/7,34/),float)
    h6=new((/7,34/),float)

    h1(0:5,:)=id_inten
    h1(6,:)=id_inten_avg
    h1(7,:)=idindex_rc
    
    h2(0:5,:)=pa_inten
    h2(6,:)=pa_inten_avg
    h2(7,:)=paindex_rc

    h3(0:5,:)=al_inten
    h3(6,:)=al_inten_avg
    h3(7,:)=alindex_rc
    

    h4(0:5,:)=width_id
    h4(6,:)=id_width_avg
    
    h5(0:5,:)=width_pa
    h5(6,:)=pa_width_avg

    h6(0:5,:)=width_alt
    h6(6,:)=al_width_avg



    time=ispan(1, 34, 1)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/纬向沃克环流的强度与宽度时间序列")
     ; wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/太平洋沃克环流的强度变化_再分析分开")
;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnYRefLine=0.0
    ; res@gsnYRefLineDashPattern=1
    res@gsnLeftString="(a)-IWCI"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@trXMaxF=34
    res@trXMinF=1

    res@trYMaxF=2.0
    res@trYMinF=-2.0

    res@tmXBValues=ispan(1,34,10)
    res@tmXBLabels=ispan(1979,2012,10)
    res@tmXBMinorValues=ispan(1,34,1)


    res@vpHeightF=0.2
    res@vpWidthF=0.4


    res@xyLineColors=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black"/)
    ; res@xyLineColors=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
    res@xyDashPattern=(/0,0,0,0,0,0,0,1/)
    res@xyLineThicknesses=(/2,2,2,2,2,2,2,2/)

    ;;;
    ;;;画图
    plot=new(6,graphic)
    plot(0)=gsn_csm_xy(wks, time, h1(0:6,:), res)


    res@gsnLeftString="(d)IWCW"
    res@trYMaxF=30
    res@trYMinF=-30

    plot(1)=gsn_csm_xy(wks, time, h4(0:6,:) , res)
     
    res@gsnLeftString="(b)PWCI"

    res@trYMaxF=2.0
    res@trYMinF=-2.0

    plot(2)=gsn_csm_xy(wks, time, h2(0:6,:) , res)
    
    res@gsnLeftString="(e)PWCW"
    res@trYMaxF=30
    res@trYMinF=-30

    plot(3)=gsn_csm_xy(wks, time, h5(0:6,:) , res)

    res@gsnLeftString="(c)AWCI"
    res@trYMaxF=1.6
    res@trYMinF=-1.6
    
    plot(4)=gsn_csm_xy(wks, time, h3(0:6,:) , res)

    res@gsnLeftString="(f)AWCW"
    res@trYMaxF=16
    res@trYMinF=-16

    plot(5)=gsn_csm_xy(wks, time, h6(0:6,:) , res)


    text=(/"20CR","ERA","JRA-55","MERRA","NCEP1","NCEP2","ENS"/)

    ; x=(/6,10,15,20,24,28,32/)
    x=(/4,8,12,17,22,28,32/)
    y=(/1.6,1.6,1.6,1.6,1.6,1.6,1.6/)
    y1=(/27,27,27,27,27,27,27/)

    ; color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
       color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black"/)
    m=new(7,graphic)
    textres=True
    ; textres@txFontHeightF =12.0
    textres@txFontHeightF =10.0

    do i=0,6
    textres@txFontColor=color(i)
    m(i)=gsn_add_text(wks,plot(0),text(i),x(i),y(i),textres)
    m(i)=gsn_add_text(wks,plot(1),text(i),x(i),y1(i),textres)
    end do
    pres=True
    gsn_panel(wks, plot, (/3,2/), pres)

    ; draw(plot(1))
    ; frame(wks)

    end
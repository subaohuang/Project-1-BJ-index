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
begin
    rea= (/"20cr","era","jra55","MERRA","ncep1","ncep2"/)

    id_inten=new((/6,34/),float)
    pa_inten=new((/6,34/),float)
    al_inten=new((/6,34/),float)


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
    delete(msf_djf)

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
  end do

    id=ploy(id_inten)
    pa=ploy(pa_inten)    
    al=ploy(al_inten)
    
    id_inten_avg=id(0,:)
    pa_inten_avg=pa(0,:)
    al_inten_avg=al(0,:)
    
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

    h1(0:5,:)=id_inten
    h1(6,:)=id_inten_avg
    h1(7,:)=idindex_rc
    
    h2(0:5,:)=pa_inten
    h2(6,:)=pa_inten_avg
    h2(7,:)=paindex_rc

    h3(0:5,:)=al_inten
    h3(6,:)=al_inten_avg
    h3(7,:)=alindex_rc

    
    time=ispan(1, 34, 1)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/纬向沃克环流的强度变化_再分析分开")
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
    plot=new(3,graphic)
    plot(0)=gsn_csm_xy(wks, time, h1(0:6,:), res)
    res@gsnLeftString="(b)PWCI"
    res@trYMaxF=2.0
    res@trYMinF=-2.0

    plot(1)=gsn_csm_xy(wks, time, h2(0:6,:) , res)
    res@gsnLeftString="(c)AWCI"
    res@trYMaxF=1.6
    res@trYMinF=-1.6
    plot(2)=gsn_csm_xy(wks, time, h3(0:6,:) , res)

    text=(/"20CR","ERA","JRA-55","MERRA","NCEP1","NCEP2","ENS"/)

    ; x=(/6,10,15,20,24,28,32/)
    x=(/4,8,12,17,22,28,32/)
    y=(/1.6,1.6,1.6,1.6,1.6,1.6,1.6/)

    ; color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
       color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black"/)
    m=new(7,graphic)
    textres=True
    ; textres@txFontHeightF =12.0
    textres@txFontHeightF =10.0

    do i=0,6
    textres@txFontColor=color(i)
    m(i)=gsn_add_text(wks,plot(0),text(i),x(i),y(i),textres)
    end do

pres=True
gsn_panel(wks, plot, (/3,1/), pres)

; draw(plot(1))
; frame(wks)

   end
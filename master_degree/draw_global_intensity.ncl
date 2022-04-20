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
begin
	diri="/mnt/f/share/rea_msf/"
	LeftString="mre"
	f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
	date=cd_calendar(f1->time, 1)
	timeind=ind(date.ge.197901.and.date.le.201212)
	level=f1->level
	levind=ind(f1->level.eq.500)
	levind2=ind(level.ge.300.and.level.le.700)
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


  ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_jja=month_to_season(msf_mean(time|:,level|:,lon|:), "JJA")
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ;;;
    
    ;;;计算冬季的强度随时间的的变化，选择3个box
    ; ;;;box1:69-90E,300-700hpa box2:180E-230E  BOX3:40w-20w
    ; inten_id=dim_avg_n_Wrap(msf_djf(:,{300:700},{60:90}),(/1,2/))
    ; inten_pa=dim_avg_n_Wrap(msf_djf(:,{300:700},{150:210}),(/1,2/))
    ; inten_al=dim_avg_n_Wrap(msf_djf(:,{300:700},{320:340}),(/1,2/))

    ;;level_mean

    inten_id= dim_avg_n_Wrap(level_mean(msf_djf(:,:,{60:90})),1)
    inten_pa=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{150:210})),1)
    inten_al=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{320:340})),1)
    ;;;
    
    ;;;计算趋势
    year=ispan(1, 34, 1)
	rc1=regCoef(ispan(1,34, 1),inten_id)
	rc2=regCoef(ispan(1,34, 1),inten_pa)
	rc3=regCoef(ispan(1,34, 1),inten_al)
	

	paindex_rc =rc2*year+rc2@yintercept
	idindex_rc =rc1*year+rc1@yintercept
	alindex_rc =rc3*year+rc3@yintercept


    ;;;读入WASwind
        diri2="/mnt/f/share/WASwind/"
        f1=addfile(diri2+"waswind_v1_0_1.monthly.nc", "r")
        date=cd_calendar(f1->time, 1)
        timeind=ind(date.ge.197901.and.date.le.201112)
        tx_time=f1->tx(timeind,:,:)
     
      
        tx_time&time=date(timeind)
        tx_mean=tx_time
        tx_djf=month_to_season(tx_mean(time|:,lat|:,lon|:), "DJF")
        tx_djf&time=ispan(year1,2011,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
        ;;注意点
        tx_djf(0:dimsizes(tx_djf&time)-2,:,:)=(/tx_djf(1:,:,:)/)
        tx_djf(dimsizes(tx_djf&time)-1,:,:)=(/tx_mean({201112},:,:)/)
        ;;;


        u_inten_id=dim_avg_n_Wrap(tx_djf(:,{-5:5},{70:100}),(/1,2/))
        u_inten_pa=dim_avg_n_Wrap(tx_djf(:,{-5:5},{150:210}),(/1,2/))
        u_inten_al=dim_avg_n_Wrap(tx_djf(:,{-5:5},{330:360}),(/1,2/))
       printVarSummary(u_inten_id)

        ;;;计算趋势
        delete(year)
        year=ispan(1, 33, 1)
        rc1=regCoef(ispan(1,33, 1),u_inten_id)
        rc2=regCoef(ispan(1,33, 1),u_inten_pa)
        rc3=regCoef(ispan(1,33, 1),u_inten_al)
        u_paindex_rc =rc2*year+rc2@yintercept
        u_idindex_rc =rc1*year+rc1@yintercept
        u_alindex_rc =rc3*year+rc3@yintercept
    ;;;
    delete(year)
 year=ispan(1, 34, 1)



    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/纬向沃克环流的强度变化")
     wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/waswind的变化")
	res=True
	res@gsnDraw=False
	res@gsnFrame=False
	res@gsnYRefLine=0.0

	res@gsnLeftString="(a)IWC"
	res@tmXTOn=False
	res@tmYROn=False
	res@tmXBMode="Explicit"
	res@trXMaxF=34
	res@trXMinF=1
	; res@trYMaxF=6
	; res@trYMinF=-6
	res@tmXBValues=ispan(1,34,10)
	res@tmXBLabels=ispan(1979,2012,10)
	res@tmXBMinorValues=ispan(1,34,1)
	res@vpHeightF=0.4
	res@vpWidthF=0.6
    res@xyLineColors=(/"green","blue","red","green","blue","red"/)
    res@xyLineThicknesses=(/2,2,2,1,1,1/)
    res@xyDashPattern=(/0,0,0,1,1,1/)
	
	plot=gsn_csm_xy(wks, year, (/inten_id,inten_pa,inten_al,idindex_rc,paindex_rc,alindex_rc/), res)

	plot1=gsn_csm_xy(wks1, ispan(1,33,1), (/u_inten_id,u_inten_pa,u_inten_al,u_idindex_rc,u_paindex_rc,u_alindex_rc/), res)
   
    ; res@tmXBValues=ispan(1,34,10)
    ; res@tmXBLabels=ispan(1979,2012,10)
    ; res@tmXBMinorValues=ispan(1,34,1)
   draw(plot)
   frame(wks)

   draw(plot1)
   frame(wks1)
   end
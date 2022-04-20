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


  ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    ; msf_jja=month_to_season(msf_mean(time|:,level|:,lon|:), "JJA")
    ; msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    ; printVarSummary(msf_son)
    ; msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    delete(msf_mean)
    ;;;
    
    ;;;计算冬季的强度随时间的的变化，选择3个box
    ; ;;;box1:69-90E,300-700hpa box2:180E-230E  BOX3:40w-20w
    ; inten_id=dim_avg_n_Wrap(msf_djf(:,{300:700},{60:90}),(/1,2/))
    ; inten_pa=dim_avg_n_Wrap(msf_djf(:,{300:700},{150:210}),(/1,2/))
    ; inten_al=dim_avg_n_Wrap(msf_djf(:,{300:700},{320:340}),(/1,2/))

    ;;level_mean

    inten_id= dim_avg_n_Wrap(level_mean(msf_djf(:,:,{70:100})),1)
    inten_pa=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{150:210})),1)
    inten_al=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{320:340})),1)
    delete(msf_djf)

    ;;;
    
    id_inten(i,:)=tofloat(inten_id)
    pa_inten(i,:)=tofloat(inten_pa)
    al_inten(i,:)=tofloat(inten_al)
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



    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/纬向沃克环流的强度变化_allrea2")
    res                    = True              ; plot mods desired
    res@gsnDraw            = False             ; don't draw yet
    res@gsnFrame           = False             ; don't advance fram
    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnYRefLine=0.0

    ; res@gsnLeftString="(a)IWC"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@trXMaxF=34
    res@trXMinF=1
    res@tmXBValues=ispan(1,34,5)
    res@tmXBLabels=ispan(1979,2012,5)
    res@tmXBMinorValues=ispan(1,34,1)
    res@vpHeightF=0.4
    res@vpWidthF=0.6
    res@xyLineColors=(/"green","blue","red","green","blue","red"/)
    res@xyLineThicknesses=(/3,3,3,1,1,1/)
    res@xyDashPattern=(/0,0,0,1,1,1/)
    res@tiYAxisString="Intensity(10~S~9~N~Kg s~S~-1)"

    mid_plot=gsn_csm_xy(wks1, year, (/id(0,:),pa(0,:),al(0,:),idindex_rc,paindex_rc,alindex_rc/), res)

    delete(res@xyLineColors)
    res@gsnXYFillColors = "green"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    res@gsnXYFillOpacities=0.1

    id_plot=gsn_csm_xy(wks1,year,id(1:2,:),res)       ; create line plot
    res@gsnXYFillColors = "blue"
    pa_plot=gsn_csm_xy(wks1,year,pa(1:2,:),res)       ; create line plot
    res@gsnXYFillColors = "red"
    al_plot=gsn_csm_xy(wks1,year,al(1:2,:),res)       ; create line plot

    ;;添加legend 
    gres = True
    gres@YPosPercent = 96.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    lineres@lgLineColors = (/"green","blue","red"/) ; line colors
    lineres@lgLineThicknesses = 5.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = (/"IWC","PWC","AWC"/)
    textres@lgLabelFontHeights = (/0.012,0.012,0.012/)                           ; label font heights



    plot = simple_legend(wks1,mid_plot,gres,lineres,textres)   





    overlay(mid_plot, id_plot)
    overlay(mid_plot, pa_plot)
    overlay(mid_plot, al_plot)
    draw(mid_plot)
    frame(wks1)
   end
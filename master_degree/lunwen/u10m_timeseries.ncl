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
;************************************************
; 画u10m的时间序列
;************************************************

begin
  rea= (/"ncep/20CR/uwnd.10m.mon.mean.nc","era/data/data_series/era_sur_global.nc",\
  "JRA-55/uwnd.10m.mon.mean.nc", "MERRA/MERRA_10m/uwnd.10m.mon.mean.nc",\
  "ncep/ncep1/uwnd.10m.mon.mean.nc" , "ncep/ncep2/uwnd.10m.mon.mean.nc","WASwind/waswind_v1_0_1.monthly.nc"/)
  
   u10_rea_id=new((/dimsizes(rea),34/),float)
   u10_rea_pa=new((/dimsizes(rea),34/),float)
   u10_rea_alt=new((/dimsizes(rea),34/),float)

do j=0,dimsizes(rea)-1,1
  title=rea(j)
  diri="/mnt/f/share/"
   print("***********************"+title+"**************************")
  
    title_uwnd=(/"ncep/20CR/uwnd.10m.mon.mean.nc","ncep/ncep1/uwnd.10m.mon.mean.nc",\
      "ncep/ncep2/uwnd.10m.mon.mean.nc","era/data/data_series/era_sur_global.nc","ICOADS/uwnd.mean.nc",\
      "WASwind/waswind_v1_0_1.monthly.nc"/)
    if any(title.eq.title_uwnd) then
    fils1=diri+title 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    if title.eq."era/data/data_series/era_sur_global.nc" then
    u_r=f1->u10(timeind,:,:)
    else if  title.eq."ncep/ncep2/uwnd.10m.mon.mean.nc" then 
    u_r=f1->uwnd(timeind,0,:,:)      
    else if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc" then
      u_r= f1->u(timeind,:,:)
    else
    u_r=f1->uwnd(timeind,:,:)
    end if
  end if
  end if

    if isatt(u_r, "add_offset") then 
      u=u_r@scale_factor*u_r + u_r@add_offset
      print("**************"+"jius"+"***")
    else 
      u=u_r
    end if

    copy_VarCoords(u_r, u)
    copy_VarAtts(u_r, u)
    u!0="time"
    u!1="lat"
    u!2="lon"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)
    delete(f1)
    delete(fils1)
    delete(date)
    delete(timeind)
    delete(u_r)
    end if
    
    if title.eq."MERRA/MERRA_10m/uwnd.10m.mon.mean.nc"

    ;;;读入MERRA的数据
    m=""+ispan(1979,2012,1)+""
    fils=new((/34,12/),string)
    do i=0,33
    if i.le.13 then
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_10m/MERRA100.prod.assim.tavgM_2d_slv_Nx."+m(i)+"*"+".SUB.nc4")
    end if
    if i.ge.14.and.i.le.21 then 
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_10m/MERRA200.prod.assim.tavgM_2d_slv_Nx."+m(i)+"*"+".SUB.nc4")
    end if
    if i.ge.22 then
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_10m/MERRA300.prod.assim.tavgM_2d_slv_Nx."+m(i)+"*"+".SUB.nc4")
    end if
    fils(i,:)= fils3
    end do
    fils_new=onedtond(fils, 34*12)
    fileh=new(408,string)
    fileh(0:407)=fils_new

    f1=addfiles(fileh, "join")
    printVarSummary(f1)
    date=cd_calendar(f1[:]->TIME, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    u=f1[:]->U10M(timeind,:,:)
  
    u!0="time"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)
    delete(f1)
    delete(fils)
    delete(date)
    end if

    if title.eq."JRA-55/uwnd.10m.mon.mean.nc" then

    fils1=diri+title 
    fils2=diri+"ncep/20CR/uwnd.10m.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)

      u_r=f1->uwnd(252:659,:,:)
    if isatt(u_r, "add_offset") then 
      u=u_r@scale_factor*u_r + u_r@add_offset
    else 
      u=u_r
    end if
    copy_VarCoords(u_r, u)
    copy_VarAtts(u_r, u)
    u!0="time"
    u!1="lat"
    u!2="lon"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)

    delete(f1)
    delete(fils1)
    delete(u_r)
    delete(date)
    end if
    
    if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc" then
    year1=1979
    year2=2011
    yearend=201112
    else
    year1=1979
    year2=2012
    yearend=201212
    end if
    
    if title.eq. "WASwind/waswind_v1_0_1.monthly.nc"  then 
      ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
      opt=(/0,5/)
      ;;先插值lon方向的
      u=(/linmsg_n(u, opt, 2)/)
      ;;插值lat方向的
      opt=(/0,1/)
      u=(/linmsg_n(u, opt, 1)/)
    ;;;mask掉降水和风场的陆地数据
    fils2=systemfunc("ls /mnt/f/share/shixiong/landsea.nc")
    a=addfile(fils2, "r")
    lsdata = a->LSMASK
    lsm = landsea_mask(lsdata,u&lat,u&lon)   ;;得到的结果，1表示land   
    u = mask(u,lsm.eq.1,False)
    u = mask(u,lsm.eq.3,False) 
    end if



    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    u_djf=month_to_season(u, "DJF")
    ; u_djf=(/month_to_annual(u, 1)/)
    
    printVarSummary(u_djf)
    delete(u_djf&time)
    u_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    ; u_djf(0:dimsizes(u_djf&time)-2,:,:)=(/u_djf(1:,:,:)/)
    ; u_djf(dimsizes(u_djf&time)-1,:,:)=(/u({yearend},:,:)/)
    ;;;
    delete(u)
    ;;;区域平均
    if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc" then
    ; u_mean_id=dim_avg_n_Wrap(u_djf(:,{-6:6},{70:100}),(/1,2/))
    ; u_mean_pa=dim_avg_n_Wrap(u_djf(:,{-6:6},{150:210}),(/1,2/))
    ; u_mean_alt=dim_avg_n_Wrap(u_djf(:,{-6:6},{320:340}),(/1,2/))    
     u_mean_id=dim_avg_n_Wrap(u_djf(:,{-4:4},{50:100}),(/1,2/))
    u_mean_pa=dim_avg_n_Wrap(u_djf(:,{-4:4},{150:210}),(/1,2/))
    u_mean_alt=dim_avg_n_Wrap(u_djf(:,{-4:4},{320:340}),(/1,2/))    

    else
    u_mean_id=dim_avg_n_Wrap(u_djf(:,{-5:5},{50:100}),(/1,2/))
    u_mean_pa=dim_avg_n_Wrap(u_djf(:,{-5:5},{150:210}),(/1,2/))
    u_mean_alt=dim_avg_n_Wrap(u_djf(:,{-5:5},{320:340}),(/1,2/))    
    end if

   ; u_mean_id=(/runave_n(u_mean_id, 3, 0, 0)/)
   ; u_mean_pa=(/runave_n(u_mean_pa, 3, 0, 0)/)
   ; u_mean_alt=(/runave_n(u_mean_alt, 3, 0, 0)/)



   delete(u_djf)
    if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc" then
    u10_rea_id(j,0:32)=tofloat(u_mean_id)
    u10_rea_pa(j,0:32)=tofloat( u_mean_pa)
    u10_rea_alt(j,0:32)=tofloat(u_mean_alt)
    else
    u10_rea_id(j,:)=tofloat( u_mean_id)
    u10_rea_pa(j,:)=tofloat( u_mean_pa)
    u10_rea_alt(j,:)=tofloat(u_mean_alt)
  end if
   delete(u_mean_id)
   delete(u_mean_pa)
   delete(u_mean_alt)
  end  do


;;;;读入WASwind
 diri1="/mnt/f/share/WASwind/waswind_v1_0_1.nc"
    id=ploy(u10_rea_id(0:5,:))
    pa=ploy(u10_rea_pa(0:5,:))    
    al=ploy(u10_rea_alt(0:5,:))
    
    id_inten_avg=id(0,:)
    pa_inten_avg=pa(0,:)
    al_inten_avg=al(0,:)
   
    was_id=u10_rea_id(6,:)
    was_pa=u10_rea_pa(6,:)
    was_alt=u10_rea_alt(6,:)



    ;;;计算趋势
    year=ispan(1, 34, 1)
    rc1=regCoef(ispan(1,34, 1),id_inten_avg)
    rc2=regCoef(ispan(1,34, 1),pa_inten_avg)
    rc3=regCoef(ispan(1,34, 1),al_inten_avg)
    

    paindex_rc =rc2*year+rc2@yintercept
    idindex_rc =rc1*year+rc1@yintercept
    alindex_rc =rc3*year+rc3@yintercept


    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/纬向沃克环流的u10m的变化")
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
    res@trYMaxF=4.0
    res@trYMinF=-8.0

    res@tmXBValues=ispan(1,34,5)
    res@tmXBLabels=ispan(1979,2012,5)
    res@tmXBMinorValues=ispan(1,34,1)
    res@vpHeightF=0.4
    res@vpWidthF=0.5
    res@xyLineColors=(/"green","blue","red","green","blue","red","green","blue","red"/)
    res@xyLineThicknesses=(/3,3,3,2,2,2/)
    res@xyDashPatterns=(/0,0,0,1,1,1/)
    res@tiYAxisString="surface zonal wind(m s~S~-1~N~)"
    res@tiXAxisFontHeightF = 0.02

    ; mid_plot=gsn_csm_xy(wks1, year, (/id(0,:),pa(0,:),al(0,:),idindex_rc,paindex_rc,alindex_rc,was_id,\
    ;   was_pa,was_alt/), res)
    
    mid_plot=gsn_csm_xy(wks1, year, (/id(0,:),pa(0,:),al(0,:),was_id,\
      was_pa,was_alt/), res)



    delete(res@xyLineColors)
    res@gsnXYFillColors = "green"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    res@gsnXYFillOpacities=0.1

    id_plot=gsn_csm_xy(wks1,year,id(1:2,:),res)       ; create line plot
    res@gsnXYFillColors = "blue"
    pa_plot=gsn_csm_xy(wks1,year,pa(1:2,:),res)       ; create line plot
    res@gsnXYFillColors = "red"
    al_plot=gsn_csm_xy(wks1,year,al(1:2,:),res)       ; create line plot
  
  ;;;WASwind画上去





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
    textres@lgLabels = (/"IWC ENS","PWC ENS","AWC ENS"/)
    textres@lgLabelFontHeights = (/0.012,0.012,0.012/)                           ; label font heights


   
    plot = simple_legend(wks1,mid_plot,gres,lineres,textres)   
    
    gres@YPosPercent = 96.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 50    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.
    
    lineres@lgDashIndexes=(/1,1,1/)
    lineres@lgLineThicknesses = 2.0                       ; line thicknesses
    textres@lgLabels = (/"IWC WASWind","PWC WASWind","AWC WASWind"/)

    plot2 = simple_legend(wks1,mid_plot,gres,lineres,textres)   

    overlay(mid_plot, id_plot)
    overlay(mid_plot, pa_plot)
    overlay(mid_plot, al_plot)



    draw(mid_plot)
    frame(wks1)
   
 end




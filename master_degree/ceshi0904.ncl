
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

;************************************************
; 画20cr,ncep,era的msf的趋势
;************************************************
function draw_all_msf(title,title2,wks,LeftString)
begin

    diri="/mnt/f/share/"
   print("***********************"+title+"**************************")
  
    title_uwnd=(/"ncep/20CR/uwnd.10m.mon.mean.nc","ncep/ncep1/uwnd.10m.mon.mean.nc",\
      "ncep/ncep2/uwnd.10m.mon.mean.nc","era/data/data_series/era_sur_global.nc","ICOADS/uwnd.mean.nc",\
      "WASwind/waswind_v1_0_1.monthly.nc"/)

     title_vwnd=(/"ncep/20CR/vwnd.10m.mon.mean.nc","ncep/ncep1/vwnd.10m.mon.mean.nc",\
      "ncep/ncep2/vwnd.10m.mon.mean.nc","era/data/data_series/era_sur_global.nc","ICOADS/vwnd.mean.nc",\
      "WASwind/waswind_v1_0_1.monthly.nc"/)

    if any(title.eq.title_uwnd) then
    fils1=diri+title 
    fils2= diri+ title_vwnd(ind(title.eq.title_uwnd))
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    timeind2=ind(date2.le.201212.and.date2.ge.197901)
    if title.eq."era/data/data_series/era_sur_global.nc" then
    u_r=f1->u10(timeind,:,:)
    v_r=f1->v10(timeind,:,:)
    else if  title.eq."ncep/ncep2/uwnd.10m.mon.mean.nc" then 
    u_r=f1->uwnd(timeind,0,:,:)      
    v_r=f2->vwnd(timeind2,0,:,:)      
    else if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc" then
      u_r= f1->u(timeind,:,:)
      v_r= f1->v(timeind,:,:)
    else
    u_r=f1->uwnd(timeind,:,:)
    v_r=f2->vwnd(timeind2,:,:)
    end if
  end if
  end if

    if isatt(u_r, "add_offset") then 
      u=u_r@scale_factor*u_r + u_r@add_offset
       v=v_r@scale_factor*v_r + v_r@add_offset
      print("**************"+"jius"+"***")
    else 
     u=u_r
     v=v_r
    end if
  
    copy_VarCoords(u_r, u)
    copy_VarAtts(u_r, u)
    u!0="time"
    u!1="lat"
    u!2="lon"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)

    copy_VarCoords(v_r, v)
    copy_VarAtts(v_r, v)
    v!0="time"
    v!1="lat"
    v!2="lon"
    delete(v&time)
    v&time=date(timeind)
    printVarSummary(v)
    delete(f1)
    delete(fils1)
    delete(date)
    delete(timeind)
    delete(v_r)
    delete(u_r)
    delete(f2)
    delete(date2)
    delete(timeind2)
    delete(fils2)
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
    v=f1[:]->V10M(timeind,:,:)
    u!0="time"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)
     v!0="time"
    delete(v&time)
    v&time=date(timeind)
    printVarSummary(v)
    delete(f1)
    delete(fils)
    delete(date)
    end if

    if title.eq."JRA-55/uwnd.10m.mon.mean.nc" then

    fils1=diri+title 
    fils2=diri+ "JRA-55/vwnd.10m.mon.mean.nc" 
    fils3=diri+"ncep/20CR/uwnd.10m.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f3=addfile(fils3, "r")

    date=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)

      u_r=f1->uwnd(252:659,:,:)
      v_r=f2->vwnd(252:659,:,:)
    if isatt(u_r, "add_offset") then 
      u=u_r@scale_factor*u_r + u_r@add_offset
      v=v_r@scale_factor*v_r + v_r@add_offset
    else 
      u=u_r
      v=v_r
    end if
    copy_VarCoords(u_r, u)
    copy_VarAtts(u_r, u)
    u!0="time"
    u!1="lat"
    u!2="lon"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)
    
     copy_VarCoords(v_r, v)
    copy_VarAtts(v_r, v)
    v!0="time"
    v!1="lat"
    v!2="lon"
    delete(v&time)
    v&time=date(timeind)
    printVarSummary(v)
    delete(f1)
    delete(fils1)
    delete(u_r)
    delete(v_r)
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
    
    if title.eq."WASwind/waswind_v1_0_1.monthly.nc".or.title.eq."ICOADS/uwnd.mean.nc"
    ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
    opt=(/0,5/)
    ;;先插值lon方向的
    u=(/linmsg_n(u, opt, 2)/)
    v=(/linmsg_n(v, opt, 2)/)
    ;;插值lat方向的
    opt=(/0,1/)
    u=(/linmsg_n(u, opt, 1)/)
    v=(/linmsg_n(v, opt, 1)/)
    delete(opt)
    fils2=systemfunc("ls /mnt/f/share/shixiong/landsea.nc")

    a=addfile(fils2, "r")
    lsdata = a->LSMASK
    lsm = landsea_mask(lsdata,u&lat,u&lon)   ;;得到的结果，1表示land   
    lsm1 = landsea_mask(lsdata,v&lat,u&lon)   ;;得到的结果，1表示land  
    u = mask(u,lsm.eq.1,False)
    u = mask(u,lsm.eq.3,False) 
    v = mask(v,lsm1.eq.1,False)
    v = mask(v,lsm1.eq.3,False) 
    end if


    ; u_pi=lonPivot(u, 30)
    ; printVarSummary(u_pi)
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    u_djf=month_to_season(u, "DJF")
    printVarSummary(u_djf)
    delete(u_djf&time)
    u_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    ;;;
    delete(u)
    ; delete(u_pi)
    ;;;计算线性趋势
    u_rc=regCoef_n(ispan(1,dimsizes(u_djf&time),1), u_djf, 0, 0)
    copy_VarCoords(u_djf(2,:,:), u_rc)
  
    ; v_pi=lonPivot(v, 30)
    ; printVarSummary(v_pi)
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    v_djf=month_to_season(v, "DJF")
    printVarSummary(v_djf)
    delete(v_djf&time)
    v_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    delete(v)
    ; delete(v_pi)
    ;;;计算线性趋势
    v_rc=regCoef_n(ispan(1,dimsizes(v_djf&time),1), v_djf, 0, 0)
    copy_VarCoords(v_djf(2,:,:), v_rc)
    delete(v_djf)




    if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc".or.title.eq."ICOADS/uwnd.mean.nc"   then
    u_rc=where(u_djf(1,:,:).eq.u_djf@_FillValue,u_rc@_FillValue , u_rc)
    end if
    delete(u_djf)

   porb_u=rc_test(u_rc)
   porb_v=rc_test(v_rc)
     if  title.eq. "WASwind/waswind_v1_0_1.monthly.nc" then      
   u_rc=where(porb_u.le.0.2.or.porb_v.le.0.1, u_rc, u_rc@_FillValue)
 else
   u_rc=where(porb_u.le.0.1, u_rc, u_rc@_FillValue)
end if
  u_rc=(/u_rc*10/)
  v_rc=(/v_rc*10/)

    ;;;;
    ;;;读入SLP的数据
    ;;;;

    
    title_uwnd=(/"ncep/20CR/prmsl.mon.mean.nc","ncep/ncep1/slp.mon.mean.nc",\
      "ncep/ncep2/mslp.mon.mean.nc","era/data/data_series/era_sur_global.nc","ICOADS/std/slp.mean.nc",\
      "Hadley/slp.mnmean.real.nc"/)

    if any(title2.eq.title_uwnd) then
    fils1=diri+title2 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    if title2.eq."ncep/ncep2/mslp.mon.mean.nc" then
    u_r=f1->mslp(timeind,:,:)
    else if  title2.eq."ncep/20CR/prmsl.mon.mean.nc"then 
    u_r=f1->prmsl(timeind,:,:)      
    else if title2.eq."era/data/data_series/era_sur_global.nc"
     u_r=f1->msl(timeind,:,:)
    else
    u_r=f1->slp(timeind,:,:)
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
    



    if title2.eq."ICOADS/std/slp.mean.nc"
    ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
    opt=(/0,5/)
    ;;先插值lon方向的
    u=(/linmsg_n(u, opt, 2)/)
    ;;插值lat方向的
    opt=(/0,1/)
    u=(/linmsg_n(u, opt, 1)/)
    print("********"+"ok"+"***********")
    delete(opt)
    ; fils2=systemfunc("ls /mnt/f/share/shixiong/landsea.nc")

    ; a=addfile(fils2, "r")
    ; lsdata = a->LSMASK
    ; lsm = landsea_mask(lsdata,u&lat,u&lon)   ;;得到的结果，1表示land   
    ; u = mask(u,lsm.eq.1,False)
    ; u = mask(u,lsm.eq.3,False) 
    end if

    delete(f1)
    delete(fils1)
    delete(timeind)
    delete(u_r)
    delete(date)
    end if
    
    if title2.eq."MERRA/MERRA_10m/slp.mean.nc"

    ;;;读入MERRA的数据
    m=""+ispan(1979,2012,1)+""
    fils=new((/34,12/),string)
    do i=0,33
    if i.le.13 then
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_slp/MERRA100.prod.assim.tavgM_2d_slv_Nx."+m(i)+"*"+".SUB.nc4")
    end if
    if i.ge.14.and.i.le.21 then 
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_slp/MERRA200.prod.assim.tavgM_2d_slv_Nx."+m(i)+"*"+".SUB.nc4")
    end if
    if i.ge.22 then
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_slp/MERRA300.prod.assim.tavgM_2d_slv_Nx."+m(i)+"*"+".SUB.nc4")
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
    u=f1[:]->SLP(timeind,:,:)
  
    u!0="time"
    delete(u&time)
    u&time=date(timeind)
    printVarSummary(u)
    delete(f1)
    delete(fils)
    delete(date)
    end if


    if title2.eq."JRA-55/prmsl.mon.mean.nc" then
   delete(date)
    fils1=diri+title2 
    fils2=diri+"ncep/20CR/prmsl.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
  
      u_r=f1->prmsl(252:659,:,:)
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
    end if
    

    year1=1979
    year2=2012
    yearend=201212
    
     if u@units.eq."Pa".or.u@units.eq."Pascals" then
       u=(/u/100.0/)
    end if
    ; u_pi=lonPivot(u, 30)
    ; printVarSummary(u_pi)
  
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    u_djf=month_to_season(u, "DJF")
    printVarSummary(u_djf)
    delete(u_djf&time)
    u_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    delete(u)
    ;;;计算线性趋势
    slp_rc=regCoef_n(ispan(1,dimsizes(u_djf&time),1), u_djf, 0, 0)
    copy_VarCoords(u_djf(2,:,:), slp_rc)
    
    prob=rc_test(slp_rc)
    delete(u_djf)
    slp_rc=(/slp_rc*10/)
     if title2.eq."ICOADS/std/slp.mean.nc"
        opt = 0
        wgt  = (/1.0, 3.0, 5.0, 3.0, 1.0/)
        wgt  = wgt/sum(wgt)             ; normalize
        slp_rc  = wgt_runave_n_Wrap(slp_rc, wgt, opt, 0)
        wgt2  = (/1.0, 2.0, 1.0/)
        wgt2  = wgt2/sum(wgt2)             ; normalize
        slp_rc   = wgt_runave_n_Wrap(slp_rc, wgt2, opt, 1)
        delete(opt)
      end if



    ;;;
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftString=LeftString

  ; res@cnFillPalette="GMT_panoply"
  res@cnFillPalette="BlueWhiteOrangeRed"
  res@cnFillColors=ispan(33,237,12)
  ; res@cnFillPalette="posneg_2"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.05
  res@lbLabelBarOn=False
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -20
  res@mpMaxLatF             = 20
  res@mpMinLonF             =0
  res@mpMaxLonF             = 360
  res@cnLineLabelsOn=False
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  res@mpShapeMode="FreeAspect"
  res@vpHeightF=0.3
  res@vpWidthF=0.6
  


  if title2.ne."Hadley/slp.mnmean.real.nc" then

  res@gsnScalarContour         = True
  res@mpFillOn              =  False          ; turn off map fill
  res@vcRefMagnitudeF         = 0.8           ; define vector ref mag
  res@vcRefLengthF            = 0.025            ; define length of vec ref
  res@vcGlyphStyle            = "LineArrow"    ; turn on curly vectors

  res@vcRefAnnoOrthogonalPosF  = -1.0
  res@vcRefAnnoPerimOn=False
  res@vcRefAnnoString2On=False
  if title2.eq."ICOADS/std/slp.mean.nc" then 
  res@vcMinDistanceF          = 0.01
   res@mpFillOn              =  True          ; turn off map fill
   res@cnFillMode = "CellFill"
 else
  res@vcMinDistanceF          = 0.025
  end if

  res@vcLineArrowHeadMaxSizeF =0.01
  res@vcLineArrowHeadMinSizeF =0.01
  plot=gsn_csm_vector_scalar_map(wks, u_rc,v_rc,slp_rc ,res)
  else 
  plot=gsn_csm_contour_map(wks, slp_rc ,res)
end if

  delete(u_rc)
  delete(v_rc)

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
  ;;;添加区域
   ;;indian
   pgx_indian=(/70,100,100,70,70/)
   pgy_indian=(/-5,-5,5,5,-5/)
   ;;;
   ;;;pacific
   pgx_pa=(/150,210,210,150,150/)
   pgy_pa=(/-5,-5,5,5,-5/)
   ;;;
   ;;;altantic
   pgx_alt=(/320,340,340,320,320/)
   pgy_alt=(/-5,-5,5,5,-5/)
   ;;;
   
   str1=unique_string("dum1")
   str2=unique_string("dum2")
   str3=unique_string("dum3")

   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   
    if title2.ne."Hadley/slp.mnmean.real.nc" then

  ; plot@$str1$ = gsn_add_polyline(wks, plot, pgx_indian, pgy_indian, pgres)
  plot@$str2$ = gsn_add_polyline(wks, plot, pgx_pa, pgy_pa, pgres)
  ; plot@$str3$ = gsn_add_polyline(wks, plot, pgx_alt, pgy_alt, pgres)
 end if
    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(1,graphic)

    iplot(0)=gsn_csm_contour(wks,prob,res_pro)
 delete(prob)

     opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.6         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger


    iplot= gsn_contour_shade(iplot, 0,0.05,opt)

    ;;;  

overlay(plot, iplot)

 return (plot)
 delete(res)
end 
;;;;
;;主程序
;;;

begin
  rea= (/"WASwind/waswind_v1_0_1.monthly.nc","WASwind/waswind_v1_0_1.monthly.nc","ncep/20CR/uwnd.10m.mon.mean.nc","era/data/data_series/era_sur_global.nc",\
 "JRA-55/uwnd.10m.mon.mean.nc", "MERRA/MERRA_10m/uwnd.10m.mon.mean.nc",\
     "ncep/ncep1/uwnd.10m.mon.mean.nc" , "ncep/ncep2/uwnd.10m.mon.mean.nc"/)

    rea2=(/"ICOADS/std/slp.mean.nc","Hadley/slp.mnmean.real.nc","ncep/20CR/prmsl.mon.mean.nc",\
       "era/data/data_series/era_sur_global.nc","JRA-55/prmsl.mon.mean.nc",\
      "MERRA/MERRA_10m/slp.mean.nc","ncep/ncep1/slp.mon.mean.nc",\
      "ncep/ncep2/mslp.mon.mean.nc"/)


   wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/太平洋U10M趋势_vector加上SLP")
   plot=new(8,graphic)
   title=(/"(a)ICOADS", "(b)HadSLP2r" ,"(c)20CR","(d)ERAIM","(e)JRA55","(f)MERRA","(g)NCEP1","(h)NCEP2"/)
   do i=0,7
    plot(i)=draw_all_msf(rea(i),rea2(i),wks,title(i))
  end do

    pres=True
    pres@gsnPanelLabelBar=True
    pres@gsnMaximize=False
  
    gsn_panel(wks, plot, (/4,2/), pres)
    
 end




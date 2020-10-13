
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
function draw_all_msf(title,wks,LeftString)
begin

    diri="/mnt/f/share/"
   print("***********************"+title+"**************************")
  


    
    title_uwnd=(/"ncep/20CR/prmsl.mon.mean.nc","ncep/ncep1/slp.mon.mean.nc",\
      "ncep/ncep2/mslp.mon.mean.nc","era/data/data_series/era_sur_global.nc","ICOADS/std/slp.mean.nc",\
      "Hadley/slp.mnmean.real.nc"/)

    if any(title.eq.title_uwnd) then
    fils1=diri+title 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    if title.eq."ncep/ncep2/mslp.mon.mean.nc" then
    u_r=f1->mslp(timeind,:,:)
    else if  title.eq."ncep/20CR/prmsl.mon.mean.nc"then 
    u_r=f1->prmsl(timeind,:,:)      
    else if title.eq."era/data/data_series/era_sur_global.nc"
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
    



    if title.eq."ICOADS/std/slp.mean.nc"
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
         ;;;区域平滑


    end if

    delete(f1)
    delete(fils1)
    delete(date)
    delete(timeind)
    delete(u_r)
    end if
    
    if title.eq."MERRA/MERRA_10m/slp.mean.nc"

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

    if title.eq."JRA-55/prmsl.mon.mean.nc" then

    fils1=diri+title 
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
    u_pi=lonPivot(u, 30)
    printVarSummary(u_pi)
  
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    u_djf=month_to_season(u_pi, "DJF")
    printVarSummary(u_djf)
    delete(u_djf&time)
    u_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    ; u_djf(0:dimsizes(u_djf&time)-2,:,:)=(/u_djf(1:,:,:)/)
    ; u_djf(dimsizes(u_djf&time)-1,:,:)=(/u({yearend},:,:)/)

    ; u_djf=(/month_to_annual(u_pi, 1)/)
    ; u_djf=(/month_to_season(u_pi, "JJA")/)
    ;;;
    delete(u)
    ;;;计算线性趋势
    u_rc=regCoef_n(ispan(1,dimsizes(u_djf&time),1), u_djf, 0, 0)
    copy_VarCoords(u_djf(2,:,:), u_rc)
    
    prob=rc_test(u_rc)
    delete(u_djf)
    u_rc=(/u_rc*10/)
     if title.eq."ICOADS/std/slp.mean.nc"
         opt = 0
        wgt  = (/1.0, 3.0, 5.0, 3.0, 1.0/)
        wgt  = wgt/sum(wgt)             ; normalize
        u_rc  = wgt_runave_n_Wrap(u_rc, wgt, opt, 0)
        wgt2  = (/1.0, 2.0, 1.0/)
        wgt2  = wgt2/sum(wgt2)             ; normalize
        u_rc   = wgt_runave_n_Wrap(u_rc, wgt2, opt, 1)
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
  ; res@cnFillPalette="BlueDarkRed18"
  res@cnFillPalette="GMT_panoply"
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
  res@mpMinLonF             =30
  res@mpMaxLonF             = 390
  res@cnLineLabelsOn=False
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  res@mpShapeMode="FreeAspect"
  res@vpHeightF=0.3
  res@vpWidthF=0.6
  ; res@mpFillDrawOrder="PostDraw"

    if title.eq."ICOADS/std/slp.mean.nc"
    res@cnFillMode="CellFill"
  else 
    res@cnFillMode="AreaFill"
  end if
  
  plot=gsn_csm_contour_map(wks, u_rc, res)
  delete(u_rc)

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
   
   
  plot@$str1$ = gsn_add_polyline(wks, plot, pgx_indian, pgy_indian, pgres)
  plot@$str2$ = gsn_add_polyline(wks, plot, pgx_pa, pgy_pa, pgres)
  plot@$str3$ = gsn_add_polyline(wks, plot, pgx_alt, pgy_alt, pgres)

    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.15       ; set contour spacing

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


    iplot= gsn_contour_shade(iplot, 0,0.15,opt)

    ;;;  

overlay(plot, iplot)

 return (plot)
end 
;;;;
;;主程序
;;;

begin

  rea=(/"ICOADS/std/slp.mean.nc","Hadley/slp.mnmean.real.nc","ncep/20CR/prmsl.mon.mean.nc",\
       "era/data/data_series/era_sur_global.nc","JRA-55/prmsl.mon.mean.nc",\
      "MERRA/MERRA_10m/slp.mean.nc","ncep/ncep1/slp.mon.mean.nc",\
      "ncep/ncep2/mslp.mon.mean.nc"/)

   wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料DJF的slp趋势")
   plot=new(8,graphic)
   title=(/"(a)ICOADS","(b)HadSLP","(c)20CR","(d)ERAIM","(e)JRA55","(f)MERRA","(g)NCEP1","(h)NCEP2"/)
   do i=0,7
    plot(i)=draw_all_msf(rea(i),wks,title(i))
  end do

    pres=True
    pres@gsnPanelLabelBar=True
    pres@lbLabelStride=2
    pres@gsnMaximize=False
    pres@pmLabelBarWidthF=0.6
    pres@pmLabelBarHeightF=0.07
    ; pres@lbTitleString="10~S~9~N~Kg s~S~-1"
    ; pres@lbTitleFontThicknessF=10.0
  
    ; pres@lbTitlePosition  = "Bottom"              ; title position
    ; pres@lbTitleFontHeightF= .02               ; make title smaller
             ; title direction
    gsn_panel(wks, plot, (/4,2/), pres)
    
 end




;#############################################
; 基于w最大值最小值定义的太平洋沃克环流的位置指数
;##############################################


function read_data(diri) ;;;读入数据
 begin 
  

  if diri.eq."era/data/data_series/w_0.75/".or.diri.eq."era/data/data_series/w_1.5/" then 

    if diri.eq."era/data/data_series/w_0.75/" then 
     fils1 = systemfunc("ls /mnt/f/share/era/data/data_series/w_0.75/era_interim_moda_"+"*"+".nc")
    else
     fils1 = systemfunc("ls /mnt/f/share/era/data/data_series/w_1.5/era_interim_moda_"+"*"+".nc")
    end if
   
    f1 = addfiles(fils1, "r")
    date=cd_calendar(f1[:]->time, 1)
    varname = getfilevarnames(f1[1])
    timeind=ind(date.le.201212.and.date.ge.197901)
    w_r = f1[:]->w(timeind,{500},{10:-10},:) 
   
  else if diri.eq."era/data/data_series/era_pre_global.nc" then 
    f1 = addfile("/mnt/f/share/"+diri , "r")
    date=cd_calendar(f1->time, 1)
    varname = getfilevarnames(f1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    w_r = f1->w(timeind,{500},{10:-10},:) 
       end if
  end if  

  if diri.eq."era/data/data_series/era_m_zonal_meri_new.nc" then 
    f1 = addfile("/mnt/f/share/"+diri , "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    w_r = f1->m_zonal(timeind,{500},{10:-10},:) 
    w_r = (/-1*w_r/)
  end if
   
  if isatt(w_r, "add_offset") then 
     w = w_r@scale_factor*w_r + w_r@add_offset
     copy_VarCoords(w_r, w)
  else 
     w = w_r
  end if
  
  if any(getVarDimNames(w).eq."longitude") then 
     w!1 = "lat"
     w!2 = "lon"
  end if
  printVarSummary(w)
  return w
  end 


function find_location(diri2,w) ;;;找出位置 diri2 : output directory
 begin 
 
  w_avg = dim_avg_n_Wrap(w, 1)

  w_avg_xi=w_avg(:,{120:240})
  w_avg_dong=w_avg(:,{180:280})


  location1=new((/dimsizes(w_avg&time)/),float)
  location2=new((/dimsizes(w_avg&time)/),float)

  do i=0,dimsizes(w_avg&time)-1,1
  location1(i)=w_avg_xi&lon(minind(w_avg_xi(i,:)))
  location2(i)=w_avg_dong&lon(maxind(w_avg_dong(i,:)))   
  end do
  
  copy_VarCoords(w_avg&time, location1)
  copy_VarCoords(w_avg&time, location2)
  
  ;;; save data
    opt=diri2+"wlocationxi.ascii"

    asciiwrite(opt, sprintf ("%5.2f,", location1))

    opt=diri2+"wlocationdong.ascii"

    asciiwrite(opt, sprintf ("%5.2f,", location2))
  ;;;

  return [/location1, location2/]
  end



;;**********************************************************
;;main 
;;********************************************************

begin 
  
   ;;;
   diri = (/"era/data/data_series/era_pre_global.nc",\
     "era/data/data_series/era_m_zonal_meri_new.nc"/)
  
   w_era = read_data(diri(0))
   w_ncep1 = read_data(diri(1))

   
   

   location_era = find_location("/mnt/f/share/rea_msf/era_",w_era)
   location_ncep1 = find_location("/mnt/f/share/rea_msf/era_",w_ncep1)
   

   location_era_xi =rmAnnCycle1D(location_era[0])  
   location_ncep1_xi =rmAnnCycle1D(location_ncep1[0])


   
   location_era_year = month_to_annual(location_era_xi, 1)
   location_ncep1_year = month_to_annual(location_ncep1_xi, 1)
   
   ; location_era_year  = (/month_to_seaDJF(location_era_xi, "DJF")/)
   ; location_ncep1_year  = (/month_to_seaDJF(location_ncep1_xi, "DJF")/)
   ; location_ncep2_year  = (/month_to_seaDJF(location_ncep2_xi, "DJF")/)
    
   location_era_year = (/dim_standardize_n(location_era_year, 1, 0)/)
   location_ncep1_year = (/dim_standardize_n(location_ncep1_year, 1, 0)/)
  
  
   location_era_xi = dim_standardize_n((/runave(location_era_xi, 3, 1)/),1,0)
   location_ncep1_xi =dim_standardize_n((/runave(location_ncep1_xi, 3, 1)/),1,0)
   ;;;
   
   ;;;read SST data

      f2 = addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
      date = cd_calendar(f2->time, 1)
      timeind = ind(date.ge.197901.and.date.le.201212)
      sst = lonFlip(f2->sst(timeind,{20:-20},:))
      ssta = (/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)
      ; ssta_de = dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
      ; copy_VarCoords(ssta, ssta_de)
      nino3_4 = dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))   
      nino3_4_std = dim_standardize_n((/runave(nino3_4, 3, 1)/),1,0)
      ; nino3_4=(/dtrend(nino3_4, True)/)
      nino3_4_anu=month_to_annual(nino3_4, 1)
      ; nino3_4_anu = (/month_to_seaDJF(nino3_4, "DJF")/)

      
      nino3_4_anu_std = dim_standardize_n(nino3_4_anu, 1, 0)
    ;;;

    ;;*************************************************** 超前滞后相关
    mxlag    = 120
    x_Lead_y = esccr(location_ncep1_xi,nino3_4,mxlag)
    y_Lead_x = esccr(nino3_4,location_ncep1_xi,mxlag)    ; switch the order of the series

    ccr = new ( 2*mxlag+1, double)    
    ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
    ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
    x=ispan(-120,120,1)     
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/w分解的部分位置与nino3.4相关性")
    resp=True
    resp@gsnXRefLine = 0
    plotp = gsn_csm_xy(wks,x,ccr,resp)
   
    ;*******************************************************************

  

   ;;;***************************************************
   ;;draw plot
   ;;;***************************************************

    time1 = ispan(1, 34*12, 1)
    time2 = ispan(1, 34, 1)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/ERA2.5与local得到的对比得到的垂直速度的对比")
    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnYRefLine=0.0
    ; res@gsnYRefLineDashPattern=1
    res@tmXTOn=False
    res@tmYROn=False

    res@tmXBMode="Explicit"
    res@trXMaxF=34*12
    res@trXMinF=1
    ; res@trYMaxF=2.0
    ; res@trYMinF=-2.0

    res@tmXBValues=ispan(1,34*12,10*12)
    res@tmXBLabels=ispan(1979,2012,10)
    res@tmXBMinorValues=ispan(1,34*12,12)


    res@vpHeightF=0.3
    res@vpWidthF=0.7


    res@xyLineColors=(/"black","blue","red"/)
    ; res@xyLineColors=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
    res@xyDashPattern=(/0,0,0,0/)
    res@xyLineThicknesses=(/3,3,3,3/)

    ;;;
    ;;;画图
    plot=new(2,graphic)
    res@gsnLeftString = "(a)w index month"
    plot(0)=gsn_csm_xy(wks, time1, (/location_era_xi,location_ncep1_xi,nino3_4_std/), res)
    res@gsnLeftString="(b)w index year"
    res@trXMaxF=34
    res@trXMinF=1
    
    res@tmXBValues=ispan(1,34,10)
    res@tmXBLabels=ispan(1979,2012,10)
    res@tmXBMinorValues=ispan(1,34,1)


    plot(1)=gsn_csm_xy(wks, time2,  (/location_era_year,location_ncep1_year,nino3_4_anu_std/) , res)

    text=(/"ERA2.5","ERA2.5_m","Nino3_4"/)
 
    
    ;;添加legend 
    gres = True
    gres@YPosPercent = 90.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 80    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
    gres@ItemSpacePercent = 10
    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"black","blue","red"/)
    lineres@lgLineThicknesses = 3.0                       ; line thicknesses
    lineres@LineLengthPercent = 5.                         ; expressed as %, 0->100, length of line
    


    textres = True
    textres@lgLabels = (/"ERA2.5","ERA2.5_m","Nino3_4"/)
    textres@lgLineLabelFontHeights= 2.0
    textres@lgLabelFontHeights = (/0.015,0.015,0.015/)                           ; label font heights
    
    do i=0,1
    plot(i) = simple_legend(wks,plot(i),gres,lineres,textres)   
    end do
    pres=True
    gsn_panel(wks, plot, (/2,1/), pres)
  ;;;

 
end
 
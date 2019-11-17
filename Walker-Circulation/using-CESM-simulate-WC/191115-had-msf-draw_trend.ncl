  ;;;************************************************************
  ;; 计算哈德来环流质量流函数
  ;; linshuheng 2019.11.17
  ;; 
  ;;;********************************************************
function read_rea_data(diri,filename,timestart,timeend)
 begin
    fils1 = diri + filename
    
    f1  =  addfile(fils1, "r")
    vnames =  getfilevarnames(f1)
    
    if any(vnames.eq."v") 
      varid = ind(vnames.eq."v")
      vatt = getfilevaratts(f1, "v")
    end if 

    if any(vnames.eq."vwnd") 
      varid = ind(vnames.eq."vwnd")
      vatt = getfilevaratts(f1, "vwnd")
    end if 
    
    date  =  cd_calendar(f1->time, 1)
    timeind  =  ind(date.le.timeend.and.date.ge.timestart)

    startyear_r = str_split_by_length(timestart, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(timeend, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    if any(vatt.eq."scale_factor") 
      v  =  short2flt(f1->$vnames(varid)$(timeind,::-1,:,:))
    else
      v  =  f1->$vnames(varid)$(timeind,::-1,:,:)
    end if

    level = f1->level(::-1)
 return [/v,level/]
 end 

function read_model_data(diri,filename,timestart,timeend)
 begin
    fils1 = diri + filename
    
    f1  =  addfile(fils1, "r")
    date  =  cd_calendar(f1->time, 1)
    timeind  =  ind(date.le.timeend.and.date.ge.timestart)

    startyear_r = str_split_by_length(timestart, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(timeend, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    v  =  f1->V(timeind,:,:,:)
    level = f1->lev
  

 return [/v,level/]
 end 

function get_had(v,level)   ;;;get Hadley mass streamfunction
 begin

    v_a = dim_avg_n_Wrap(v,3)
    a = 6371000
    g = 9.8
    pi = get_pi("double")
    dfi = pi

    if any(getVarDimNames(v_a).eq."lat") 
      cosfi = cos(v_a&lat/180.0)
    else    
      cosfi = cos(v_a&latitude/180.0)
    end if
    
    cosfi_new = conform_dims(dimsizes(v_a), cosfi, 2)


   ;;; vibeta(p, x, linlog, psfc, pbot, ptop)积分纬向平均的结果，得到mass
    msf = v_a
    psfc = v_a(:,1,:)
    psfc = 1013.0
   
    do i = 0,dimsizes(level) -1
      if  any(getVarDimNames(v_a).eq."lev") 
       msf(:,i,:)  = (/vibeta( level, v_a(time|:,lat|:,lev|:), 1, psfc, level(i), 0)/)
      else  
       msf(:,i,:)  = (/vibeta( level, v_a(time|:,latitude|:,level|:), 1, psfc, level(i), 0)/)
      end if
    end do

    msf_fi = v_a
    msf_fi := (/((2*a*dfi*cosfi_new)/g)*msf*100/)
    msf_fi = (/msf_fi/10^10/)
    copy_VarCoords(v_a, msf_fi)

   ;;;
  return msf_fi
 end 

function get_trend(fi_annual)
 begin
  ;;;trend 
  rc1 = regCoef_n(ispan(1, dimsizes(fi_annual(:,1,1)), 1), fi_annual, 0, 0)
   
  rc1 = (/rc1*10/)  ;;;convert to decadal trend mean

  copy_VarCoords(fi_annual(0,:,:), rc1)

  ;;;
  return rc1
 end 

function  draw_msf_trend(diri,filename,timestart,timeend,wks,LeftString,mode)
  begin 
    

    if mode.eq.0 then 
      data_CT = read_model_data(diri,filename,timestart,timeend)
    else
      data_CT = read_rea_data(diri,filename,timestart,timeend)   
    end if

    msf_CT = get_had(data_CT[0],data_CT[1])
   
    msf_CT_annual = month_to_annual(msf_CT, 1)

    ;;;; linear trend
    rc_ct = get_trend(msf_CT_annual)
    ;;;;

    ;;;get clim
    msf_ct_clm = dim_avg_n_Wrap(msf_CT, 0)
    ;;;
    
    printVarSummary(rc_ct)
    printVarSummary(msf_ct_clm)

    if any(getVarDimNames(rc_ct).eq."lev")
       rc_ct&lev@units = "hPa" 
       rc_ct&lev@long_name = "pressure_level"
    else 
       rc_ct&level@units = "hPa"
       rc_ct&level@long_name = "pressure_level"
       
    end if
  
   printVarSummary(rc_ct&level)
   




    ;;;;;;画图
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlRe"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.4
    res@cnMaxLevelValF = 0.4
    res@cnLevelSpacingF = 0.1

    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@lbLabelBarOn = True
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
  
    ; res@tmXBValues=(/0,90,180,270,357.5/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-90,90,10)

    ; res@tmXBValues=ispan(100,310,40)
    ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","60~S~o~N~W"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(100,310,10)
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7

    
    ; res@cnFillPalette = "MPL_seismic"

    res@gsnLeftString = LeftString
    plot_diff = gsn_csm_pres_hgt(wks,rc_ct,res)
    exit()
    
    ;;;;添加气候态的等值线
    res2 =res
    res2@gsnRightString=""
    res2@gsnLeftString=""
    res2@cnLevelSelectionMode = "ExplicitLevels"
    res2@cnLinesOn=True
    res2@cnLevels = ispan(-10,10,2)
    res2@cnLineLabelsOn=True
    res2@cnLineLabelBackgroundColor="Transparent"
    res2@cnLineThicknessF = 2.
    res2@cnFillOn = False
    res2@gsnContourNegLineDashPattern=1
    res2@gsnContourZeroLineThicknessF=4.0
    
    
    plot_clm = gsn_csm_pres_hgt(wks, msf_ct_clm, res2)
    
    overlay(plot_diff, plot_clm)

 return plot_diff
 delete(plot_diff)
 delete(plot_clm)
 end 
;;*****************************************
;;;;main program 
;;********************************************
;;;;;;
 begin
    diri1 =  "/home/ys17-19/lsh/data/era/"
    diri2 = "/home/ys17-19/lsh/CESM-data/"
    
    filename1 = "F_1850-PDAY/F_1850-PDAY.cam.h0.197801-201212_chazhi.nc"
    filename2 = "F_1850-PDAY-rmindian/F_1850-PDAY-rmindian.cam.h0.197801-201212_chazhi.nc"
    filename3 = "F_1850-PDAY-rmatlantic/F_1850-PDAY-rmatlantic.cam.h0.197801-201212_chazhi.nc"
    filename4 = "F_1850-PDAY-rmpacific/F_1850-PDAY-rmpacific.cam.h0.197801-201212_chazhi.nc"  
    filename5 = "F_1850-PDAY-rafixed/F_1850-PDAY-rafixed.cam.h0.197801-201212_chazhi.nc"
    filename6 = "era_pre_global.nc"

    timestart = 197901
    timeend = 201212

    wks = gsn_open_wks("pdf", "CESM各试验模拟的哈德来环流的趋势")

    plot_era = draw_msf_trend(diri1,filename6,timestart,timeend,wks,"ERA-Interim",1)
    plot_CT = draw_msf_trend(diri2,filename1,timestart,timeend,wks,"Control run",0)
    plot_rmid = draw_msf_trend(diri2,filename2,timestart,timeend,wks,"indian warming",0)
    plot_rmal = draw_msf_trend(diri2,filename3,timestart,timeend,wks,"altlantic warming",0)
    plot_rmpa = draw_msf_trend(diri2,filename4,timestart,timeend,wks,"pacific SST trend",0) 
    plot_rafixed = draw_msf_trend(diri2,filename5,timestart,timeend,wks,"GHG",0)
    

    
    pres = True
    pres@gsnPanelLabelBar = True
    
    gsn_panel(wks, plots, (/3,2/), pres)
  
 end 
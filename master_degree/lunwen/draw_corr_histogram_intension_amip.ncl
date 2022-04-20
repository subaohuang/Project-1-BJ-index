function  ploy(locationxi_rea,locationxi_amip)
begin
mnmx      = new ((/4,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
dim1=dimsizes(locationxi_amip&time)
mnmx(0,:) = dim_avg_n_Wrap(locationxi_rea, 0)
mnmx(1,0:dim1-1) = dim_avg_n_Wrap(locationxi_amip, 0)
mnmx(2,0:dim1-1) = dim_max_n_Wrap(locationxi_amip, 0)
mnmx(3,0:dim1-1) = dim_min_n_Wrap(locationxi_amip, 0)
return mnmx
end 
;;***********************************************************************
;;主程序
;;;**********************************************************************
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"intension_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    inten_rea_pa=f1->intension_pa(:,timeind)
    inten_rea_id=f1->intension_id(:,timeind)
    inten_rea_alt=f1->intension_alt(:,timeind)
   
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;;read amip
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    f1=addfile(diri+"intension_amip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)

    inten_amip_pa=f1->intension_pa(:,timeind)
    inten_amip_id=f1->intension_id(:,timeind)
    inten_amip_alt=f1->intension_alt(:,timeind)
    
    
     inten_amip_id = (/-1*inten_amip_id/)


    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;
printVarSummary(inten_amip_id)
printVarSummary(inten_amip_pa)
corr = escorc_n(inten_amip_id, inten_amip_pa , 1 ,1)
corr(18) = 0.75



;;;************************************************************
;画图
;;;**********************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/印度洋强度与大西洋强度相关系数的模式比较")

;********************************************************   
  res                       = True                 ; plot mods desired         
  res@gsnFrame              = False                ; don't advance frame yet
  res@gsnXYBarChart         = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  res@gsnXYBarChartColors= (/"SkyBlue"/)            ; choose colors
  
  ; res@tmXBOn                = False                ; turn off tickmarks at bot
  res@trYMinF               = min(corr)                    ; bring bars down to zero
  res@trXMinF               = 0                    ; adds space on either end
  res@trXMaxF               =  28                   ; of the 1st and last barsr
  res@tmXTOn  = False
  res@tmYROn = False
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(1, 27, 1)
  res@tmXBLabelFont = 0
  

  res@tmXBLabels = (/"M","a","b","c","d","e","f","g","h","i","z","k","l","m","n","o","p","q","r","s","t",\
                     "u","v","w","x","y","z"/)


  res@tmXBLabelFontHeightF = 0.015
  res@tmXBMinorOn = False
  
  res@tiMainString          = "correlation coefficient "
  res@vpHeightF = 0.3
  res@vpWidthF = 0.7



   x= ispan(1, 27, 1)
  plot = gsn_csm_xy (wks,x,corr,res)                  ; create plot

end
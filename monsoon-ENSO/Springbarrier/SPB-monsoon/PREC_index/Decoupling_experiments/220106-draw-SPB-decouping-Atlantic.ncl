;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
function spb_cal(data)  ;;  data:two dimension 12×yearnum array i.e.12×10
;;data1 : the monsoon anomous year  data2: next year
begin 
  dim = dimsizes(data)
  yearnum = dim(1)
  ;;; cal ACF, want to get a two dimensional array 12×13,initial month start
     acr = new((/12,13/),float)
     do i=0,11
         do j=0,12
            k=i+j 
            if k.le.11 then 
              acr(i,j) = escorc(data(i,0:yearnum-1),data(i+j,0:yearnum-1))
            else 
              acr(i,j) = escorc(data(i,0:yearnum-2),data(i+j-12,1:yearnum-1))
            end if 
         end do 
     end do
   ;;;
 ;;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   ; acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_positive = -1*acr_gra
   acr_gra_max = dim_max_n(acr_gra_positive, 1) 
 
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"

   maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   SPB = dim_sum_n_Wrap(acr_gra_max, 0)
   ; acr@maind = maind
   acr@SPB = SPB
; acr_gra_mar2aug = new(12,float)
  
;   do i=0,11
;     if i.le.2 then 
;     acr_gra_mar2aug(i) = sum(acr_gra(i,6-i+3:12))
;     end if

;     if i.le.9.and.i.gt.2 then 
;     acr_gra_mar2aug(i) = sum(acr_gra(i,6-i+3:6-i+3+5))
;     end if
;     if i.gt.9 then 
;     acr_gra_mar2aug(i) = sum(acr_gra(i,0:5+9-i))
;     end if
;   end do 
;    SPB = dim_avg_n_Wrap(acr_gra_mar2aug, 0)
   acr_gra@SPB = SPB
return acr
end 

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin


  startmon = 25601 
  endmon = 28012
;;;读入控制实验的数据
    ; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CP_CTRL/"
    ; fils1=systemfunc("ls "+diri+"B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc")

     diri="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CP_CTRL/"
    fils1=systemfunc("ls "+diri+"B2000_f19g16_CP_CTRL.TEMP.h1.0251-0340.nc")
    
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    ; lat = f1->lat
    ; lon = f1->lon

    ; timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 海温数据从1977年开始
    timeind:=ind(date.le.34012.and.date.ge.25601)  ;;; 海温数据从1977年开始


    sst_ctrl_r = f1->TEMP(timeind,0,:,:)

     

    sst_ctrl = PopLatLon(sst_ctrl_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_ctrl)


    
    ssta_ctrl = (/rmMonAnnCycTLL(sst_ctrl)/)
    copy_VarCoords(sst_ctrl, ssta_ctrl)

    ;;detrend 
    ssta_ctrl = (/dtrend_msg_n(ispan(1,dimsizes(ssta_ctrl&time),1),ssta_ctrl,False,True,0)/)
    ;;;

    delete(date)
    delete(timeind)
    
    nino34_ctrl =  dim_avg_n_Wrap(ssta_ctrl(:,{5:-5},{190:240}), (/1,2/))
;;;;
  


;;;split to 12 month nino34_12 is a array size (12,yearnum)

    nino34_ctrl_12 = split_12(nino34_ctrl)

;;;;



;;;读入decoupling实验的数据
    ; diri="/home/yangsong3/data-model/lsh/CESM/B/B2000_alt_north_year_CTRL/"
    ; fils1=systemfunc("ls "+diri+"B2000_alt_north_year_CTRL.TEMP.fvgrid.0251-0280.nc")

    diri="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/"
    fils1=systemfunc("ls "+diri+"B2000_alt_north_year_CTRL.TEMP_taux_tauy.h1.0251-0280.nc")

    ; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year/"
    ; fils1=systemfunc("ls "+diri+"B2000_alt_north_year.TEMP_taux_tauy.h1.0251-0280.nc")
    

    ; diri="/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_db_year_80_20/"
    ; fils1=systemfunc("ls "+diri+"lsh_B2000_alt_north_db_year_80_20.TEMP.h1.0251-0280.nc")

    
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    ; lat = f1->lat
    ; lon = f1->lon

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 海温数据从1977年开始

    sst_r = f1->TEMP(timeind,0,:,:)
    
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    ssta = (/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;

    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

;;;;
  
;;;split to 12 month nino34_12 is a array size (12,yearnum)
  nino34_12 = split_12(nino34)
;;;;



;;;
  spb_ctrl = spb_cal(nino34_ctrl_12)
  spb_decoupling = spb_cal(nino34_12)
;;;


 print("ctrl="+spb_ctrl@SPB)
 print("decoupling_atlantic="+spb_decoupling@SPB)
 
  plot = new(2,graphic)


wks = gsn_open_wks("pdf", "SPB_atlantic_decouping")

;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=False
    res@gsnAddCyclic = False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn= True
    ; res@pmLabelBarOrthogonalPosF = 0.2
    ; res@lbLabelStride=1
    res@lbOrientation="Vertical"
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmBorderThicknessF=1.0
    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF = 1.5

    res@tmXBValues = ispan(0,12,1)
    res@tmXBLabels = ispan(0,12,1) +""

    res@tmYLValues = ispan(0,11,1)
    ; res@tmYLLabels = (/"Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May"/) 
    res@tmYLLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/) 

    ; res@cnFillPalette = "BlueDarkOrange18"
    ; gsn_define_colormap(wks,"BlueDarkRed18")

     ; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

     gsn_define_colormap(wks,"BlueDarkRed18")


    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels

    res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 


    ; res@cnLevels = (/-0.25,-0.2,-0.15,-0.1,-0.05,0.0/)
    ; 
    ; res@cnFillColors = (/200,185,176,161,145,0,0/) 
    ; res@cnFillColors = (/17,33,49,65,81,0,0,129/)   
    res@lbLabelStride=1

   res@gsnLeftString = "CTRL"
   res@gsnRightString = ""+spb_ctrl@SPB*-1
   plot(0) = gsn_csm_contour(wks,spb_ctrl,res)

   res@gsnLeftString = "Decoupling Atlantic"
   res@gsnRightString = ""+spb_decoupling@SPB*-1
   plot(1) = gsn_csm_contour(wks,spb_decoupling,res)
  
   ; resp = True
   ; resp@gsMarkerColor = "black"
   ; resp@gsMarkerSizeF = 0.1
   ; resp@gsMarkerIndex = 1

   ; plot_gs = gsn_add_polymarker(wks, plot(0), spb_strmon@maind, ispan(0,11,1), resp)
   ; plot_gs1 = gsn_add_polymarker(wks, plot(1), spb_weakmon@maind, ispan(0,11,1), resp)
   ; plot_gs2 = gsn_add_polymarker(wks, plot(2), spb_normon@maind, ispan(0,11,1), resp)

   x1 = ispan(0,2,1)
   y1 = x1 ;;; line for mar

   x2 = ispan(0,7,1)
   y2 = x2 ;;; line for AUG
   
   x3 = ispan(3,12,1)
   y3 = ispan(2,11,1) ;;; line for MAR


   resp = True
   resp@gsLineColor = "green"
   resp@gsLineThicknessF = 6.0
   
  
   plot_gs = gsn_add_polyline(wks, plot(0), x1(::-1), y1, resp)
   plot_gs1 = gsn_add_polyline(wks, plot(1),  x1(::-1), y1, resp)


   ; plot_gs2 = gsn_add_polyline(wks, plot(2),  x1(::-1), y1, resp)
   resp@gsLineColor = "blue"
   plot_gs_aug = gsn_add_polyline(wks, plot(0), x2(::-1), y2, resp)
   plot_gs_aug1 = gsn_add_polyline(wks, plot(1),  x2(::-1), y2, resp)
   ; plot_gs_aug2 = gsn_add_polyline(wks, plot(2),  x2(::-1), y2, resp)

   resp@gsLineColor = "green"
   plot_gs_augs = gsn_add_polyline(wks, plot(0), x3, y3(::-1), resp)
   plot_gs_augs1 = gsn_add_polyline(wks, plot(1),  x3, y3(::-1), resp)
   ; plot_gs_augs2 = gsn_add_polyline(wks, plot(2),  x3, y3(::-1), resp)

   rest = True
   rest@txFontColor = "black"
   rest@txFontHeightF = 0.02
   
    plot_gt = gsn_add_text(wks, plot(0),"MAR", 2 ,1, rest)
    plot_gt1 = gsn_add_text(wks, plot(1),"MAR", 2 ,1, rest)
    ; plot_gt2 = gsn_add_text(wks, plot(2),"MAR", 5 ,5, rest)
  
    
    plot_gx = gsn_add_text(wks, plot(0),"AUG", 4.2 ,5.5, rest)
    plot_gx1 = gsn_add_text(wks, plot(1),"AUG", 4.2 ,5.5, rest)
    ; plot_gx2 = gsn_add_text(wks, plot(2),"AUG", 8.2 ,6.5, rest)
    

    plot_gxs = gsn_add_text(wks, plot(0),"MAR", 8 ,6.5, rest)
    plot_gxs1 = gsn_add_text(wks, plot(1),"MAR", 8 ,6.5, rest)
    ; plot_gxs2 = gsn_add_text(wks, plot(2),"AUG", 2 ,1, rest)
   pres = True
   gsn_panel(wks, plot, (/2,1/), pres)
   ; draw(plot)
   ; frame(wks)
   ;;;;;;; draw root_mean_square

    ;   resxy = True
    ;   resxy@gsnDraw = False
    ;   resxy@gsnFrame = False
    ;   resxy@trXMinF = 0
    ;   resxy@trXMaxF = 12

    ;   resxy@trYMaxF = 1.0
    ;   resxy@trYMinF = 0.0
      
    ;   resxy@xyMonoLineColor = False
    ;   resxy@xyLineColors = (/"red","blue","black"/)
    ;   resxy@xyDashPatterns= (/0,0,0/)
    ;   resxy@xyMonoLineThickness = True
    ;   resxy@xyLineThicknessF = 1.5

    ;   resxy@gsnLeftString = "R "
    ;   resxy@vpHeightF = 0.4
    ;   resxy@vpWidthF = 0.6
    ;   plotxy = new(1,graphic)
    ;   plotxy(0) = gsn_csm_xy(wks1, ispan(0,12,1), (/R_strmon,R_weakmon,R_normon/) , resxy) 
    ;   ;;添加legend 
    ;   gres = True
    ;   gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;     ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ;   gres@XPosPercent = 60    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ;   lineres = True
    ;   ; lineres@lgLineColors = (/"black"，/) ; line colors
    ;   lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
    ;   lineres@lgLineThicknesses = resxy@xyLineThicknessF                       ; line thicknesses
    ;   lineres@LineLengthPercent = 9.
    ;   lineres@lgDashIndexes   = resxy@xyDashPatterns                     ; expressed as %, 0->100, length of line

    ;   textres = True
    ;   textres@lgLabels = (/"Str monsoon","Weak monsoon","Normal monsoon"/)
    ;   ; textres@lgLineLabelFontHeights = (/0.05,
    ;   ; textres@lgLabelOpacity = (/0/) 

    ;   plotxy(0) = simple_legend(wks1,plotxy(0),gres,lineres,textres)
     
    ; draw(plotxy)
    frame(wks)
 
   end 

;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
function spb_cal(data,data2)  ;;  data:two dimension 12×yearnum array i.e.12×10
;;data1 : the monsoon anomous year  data2: next year
begin 
  dim = dimsizes(data)
  yearnum = dim(1)
  ;;; cal ACF, want to get a two dimensional array 12×13,initial month start June 
     acr = new((/12,13/),float)
     do i=0,11
         do j=0,12
            k=i+j 
            if k.le.11 then 
              acr(i,j) = escorc(data(i,0:yearnum-1),data(i+j,0:yearnum-1))
            else 
              acr(i,j) = escorc(data(i,0:yearnum-1),data2(i+j-12,0:yearnum-1))
            end if 
         end do 
     end do
   ;;;
 ;;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"
   maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   SPB = dim_sum_n_Wrap(acr_gra_max, 0)
   acr@maind = maind
   acr@SPB = SPB
return acr
end 
function split_12(nino34)
begin
    yearnum = dimsizes(nino34)/12
    nino34_12 = new((/12,yearnum/),float)

    ; do i=0,11   
    ; if i.le.6 then 
    ; nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
    ; else
    ; nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
    ; end if
    ; end do 
  do i=0,11   
    nino34_12(i,:) = nino34(i:dimsizes(nino34)-1:12)
  end do 
return nino34_12
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin

 ;;;time span
    startmon =25601
    endmon = 28012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;


;;;读入控制实验组的海温数据
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl = dim_avg_n_Wrap(f1->TEMP(timeind,{0:40000},{-5:5},:),2)

    sst_ctrl&z_t = (/sst_ctrl&z_t/100.0/)


    ssta_ct = (/rmMonAnnCycTLL(sst_ctrl)/)
    copy_VarCoords(sst_ctrl, ssta_ct)

    ;;detrend 
    ssta_ct = (/dtrend_msg_n(ispan(1,dimsizes(ssta_ct&time),1),ssta_ct,False,True,0)/)
    ;;;
    nino34_ct =  dim_avg_n_Wrap(ssta_ct(:,{5:-5},{190:240}), (/1,2/))
 
;;;读入B22000 WNPCLM海温数据
    
    ; model_name = (/"HEAT_B_BP"/)
    model_name =(/"lsh_B2000_WNP_heating_Idl"/)
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.TEMP.fvgrid.0251-0280.nc"
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_IDLthree/lsh_B2000_WNP_heating_IDLthree.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"+model_name+"/"+model_name+".TEMP.fvgrid.0251-0280.nc"
    
    f2 = addfile(fils2, "r")
    time := f2->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)


    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    sst_wnp = dim_avg_n_Wrap(f2->TEMP(timeind,{0:40000},{-5:5},:) ,2)
    sst_wnp&z_t = (/sst_wnp&z_t/100.0/)
    
    ssta_wnp = (/rmMonAnnCycTLL(sst_wnp)/)
    copy_VarCoords(sst_wnp, ssta_wnp)

    ;;detrend 
    ssta_wnp = (/dtrend_msg_n(ispan(1,dimsizes(ssta_wnp&time),1),ssta_wnp,False,True,0)/)
    nino34_wnp =  dim_avg_n_Wrap(ssta_wnp(:,{5:-5},{190:240}), (/1,2/))
;;;;;


  nino34_ct_12 = split_12(nino34_ct)
  nino34_wnp_12 = split_12(nino34_wnp)
  dim = dimsizes(nino34_ct_12)
  spb_ct =  spb_cal(nino34_ct_12(:,0:dim(1)-2),nino34_ct_12(:,1:dim(1)-1))
  spb_wnp =  spb_cal(nino34_wnp_12(:,0:dim(1)-2),nino34_wnp_12(:,1:dim(1)-1))

 
plot = new(2,graphic)
wks = gsn_open_wks("pdf", "nino34-acr_hov")
; wks1 = gsn_open_wks("pdf", "ISM-after-root_mean_square-newmonsoon.pdf")


;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=True
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
    gsn_define_colormap(wks,"BlueDarkRed18")


    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels

    res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 
  ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

   res@gsnLeftString = "CTRL"
    res@gsnRightString = ""+spb_ct@SPB
   plot(0) = gsn_csm_contour(wks,spb_ct,res)
   res@gsnLeftString = "weakmonsoon"
   res@gsnRightString = ""+spb_wnp@SPB
   plot(1) = gsn_csm_contour(wks,spb_wnp,res)
  

   resp = True
   resp@gsMarkerColor = "black"
   resp@gsMarkerSizeF = 0.1
   resp@gsMarkerIndex = 1

   plot_gs = gsn_add_polymarker(wks, plot(0), spb_ct@maind, ispan(0,11,1), resp)
   plot_gs1 = gsn_add_polymarker(wks, plot(1), spb_wnp@maind, ispan(0,11,1), resp)
   ; plot_gs2 = gsn_add_polymarker(wks, plot(2), spb_normon@maind, ispan(0,11,1), resp)

   pres = True
   gsn_panel(wks, plot, (/2,1/), pres)
   ; draw(plot)
   ; frame(wks)
   ; ;;;;;;; draw root_mean_square

   ;    resxy = True
   ;    resxy@gsnDraw = False
   ;    resxy@gsnFrame = False
   ;    resxy@trXMinF = 0
   ;    resxy@trXMaxF = 12

   ;    resxy@trYMaxF = 1.0
   ;    resxy@trYMinF = 0.0
      
   ;    resxy@xyMonoLineColor = False
   ;    resxy@xyLineColors = (/"red","blue","black"/)
   ;    resxy@xyDashPatterns= (/0,0,0/)
   ;    resxy@xyMonoLineThickness = True
   ;    resxy@xyLineThicknessF = 1.5

   ;    resxy@gsnLeftString = "R "
   ;    resxy@vpHeightF = 0.4
   ;    resxy@vpWidthF = 0.6
   ;    plotxy = new(1,graphic)
   ;    plotxy(0) = gsn_csm_xy(wks1, ispan(0,12,1), (/R_strmon,R_weakmon,R_normon/) , resxy) 
   ;    ;;添加legend 
   ;    gres = True
   ;    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
   ;      ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
   ;    gres@XPosPercent = 60    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

   ;    lineres = True
   ;    ; lineres@lgLineColors = (/"black"，/) ; line colors
   ;    lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
   ;    lineres@lgLineThicknesses = resxy@xyLineThicknessF                       ; line thicknesses
   ;    lineres@LineLengthPercent = 9.
   ;    lineres@lgDashIndexes   = resxy@xyDashPatterns                     ; expressed as %, 0->100, length of line

   ;    textres = True
   ;    textres@lgLabels = (/"Str monsoon","Weak monsoon","Normal monsoon"/)
   ;    ; textres@lgLineLabelFontHeights = (/0.05,
   ;    ; textres@lgLabelOpacity = (/0/) 

   ;    plotxy(0) = simple_legend(wks1,plotxy(0),gres,lineres,textres)
     
   ;  draw(plotxy)
    ; frame(wks)
   end 

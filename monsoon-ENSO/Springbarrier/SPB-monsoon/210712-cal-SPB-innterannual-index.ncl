;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
function spb_cal(data)  ;;  data:two dimension 12×yearnum array i.e.12×10
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
              acr(i,j) = escorc(data(i,0:yearnum-2),data(i+j-12,1:yearnum-1))
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

return SPB
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 197901
    endmon = 201912

    startyear = 1979
    endyear = 2019
;;;;

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    ; printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
;;;

;;;split to 12 month nino34_12 is a array size (12,yearnum) monsoon year 
    yearnum = dimsizes(nino34)/12 -1
    nino34_12 = new((/12,yearnum/),float)
;;;;


;;;;
  do i=0,11   
    if i.le.6 then 
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
     else
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)   ;;;;first year 1900 JJA -1901 MAM
    end if
  end do 

;;;;


;;; total year SPB
  spb_total = spb_cal(nino34_12)
;;;



;;; cal Every year
spb_year = new(yearnum,float)

ind_all = ispan(0, yearnum-1, 1)
ind_all@_FillValue =  99

spb_year@_FillValue = nino34_12@_FillValue


do i=0,yearnum-1,1
  indx = get1Dindex_Exclude(ind_all, i)
  spb_year(i) = spb_cal(nino34_12(:,indx))
  spb_year(i) = (/spb_total - spb_year(i)/)
  delete(indx)
end do

spb_year = dim_standardize(spb_year, 0)

opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual.ascii"    
asciiwrite(opt1, sprintf("%6.2f,",spb_year))  


;;; plot 
;;------------------------
;;------------------------
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.9
  res@vpHeightF = 0.4
  res@vpXF      = 0.06
  res@vpYF      = 0.75


  res@tiXAxisOn      = False
  res@tiYAxisOn      = False
  res@tmXTLabelsOn   = False
  res@tmXTOn         = False

  res@trXMinF        = 0
  res@trXMaxF        = yearnum

  res@trYMinF        = -4.0
  res@trYMaxF        = 4.0


  res@tmXBLabelFontHeightF = 0.017
  res@tmYLLabelFontHeightF = 0.017
  res@gsnLeftStringFontHeightF   = 0.023
  res@gsnRightStringFontHeightF  = 0.023

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(1, yearnum, 5)
  res@tmXBLabels = ispan(startyear+1, endyear, 5)

  res@tmXBMinorOn = True
  res@tmXBMinorValues = ispan(1, yearnum, 1)
  res@tmYLMode = "Explicit"
  res@tmYLValues = fspan(-5.0, 5.0, 11)
  res@tmYLLabels = fspan(-5.0, 5.0, 11)
  res@tmYLMinorOn = True
  res@tmYLMinorValues = fspan(-5.0, 5.0, 21)
  bres = res


  res@xyMarkLineModes = (/"MarkLines","Lines"/)
  res@xyMarkers = (/1,1/)
  res@xyMarkerSizes = (/0.03,0.001/)
  res@xyMarkerColors= (/"black","black"/)
  res@xyMarkerThicknesses = (/4.0,2.0/)
  res@xyLineThicknesses  = (/4.0,4.0/)
  res@xyLineColors       = (/"black","black"/)
  res@xyDashPatterns     = (/0,4/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018





  

  ;;------------------------ FIG.2 ----------------------------
  pltType = "pdf"
  pltName = "/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/SPB-interannual-index"
  wks = gsn_open_wks(pltType, pltName)
  res@xyLineColors  = (/"black","black"/)
  res@xyMarkerColors= (/"black","black"/)
  plot  = gsn_csm_xy (wks, ispan(1, yearnum, 1), spb_year, res)

    ;;---------------- 参考线 --------------
  x=ispan(1, yearnum, 1)
  y= x*0.1
  y = 0.0
   lnres   =    True
   lnres@gsLineColor="gray30"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=0
   lnres@gsnFrame = False
   dum000 = gsn_add_polyline(wks,plot,x,y,lnres)
  y= 0.5
   lnres@gsLineColor="red"
   lnres@gsLineDashPattern=11
   dum111 = gsn_add_polyline(wks,plot,x,y,lnres)
  y=-0.5
   lnres@gsLineColor="blue"
   dum222 = gsn_add_polyline(wks,plot,x,y,lnres)

  ;;---------------- 标题 --------------
  txres               = True
  txres@txFontHeightF = 0.019
  txres@txFont = "helvetica-bold"
  txres@txFontColor   = "black"
  gsn_text_ndc(wks, "Spring prediction barrier index",  0.5, 0.775, txres)
  
  draw(plot)
  frame(wks)

print("****done******")
;;;;

end 


   

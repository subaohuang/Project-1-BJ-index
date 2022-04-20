;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
function acr_cal(data,data2)  ;;  data:two dimension 12×yearnum array i.e.12×10
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
return acr
end 

function spb_cal(acr,method)
begin
   
   if method.eq.1 then 
   ;;method1 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
   spb =  dim_sum_n_Wrap(acr_gra_max, 0)
   end if

   if method.eq.2 then 
   ;;method2  3-6月梯度
    acr_gra = acr(:,1)
    acr_gra(7) = (acr(7,2) - acr(7,5))/3  ;;3月-6月 7是对应1月
    acr_gra(8) = (acr(8,1) - acr(8,4))/3            ;;2月
    acr_gra(9) = (acr(9,0) - acr(9,3))/3               ;;3月
    acr_gra(10) = (acr(10,0) - acr(10,2))/2            ;;4月      
    acr_gra(11) = (acr(11,11) - acr(11,12))/1            ;;5月  次年4月-5月
    do i =0,6
     acr_gra(i) = (acr(i,11-5-i+3) - acr(i,11-5-i+6))/2
    end do 
    spb =  dim_sum_n_Wrap(acr_gra, 0)
   end if

   if method.eq.3 then 
   ;;;method3 
      acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
      acr_gra = -1*acr_gra  ;;;change to positive

      acr_gra_max = acr_gra(:,1)

      do i=0,9
      acr_gra_max(i) = dim_max_n(acr_gra(i,11-5-i+4:11-5-i+6), 0)
      end do 

      acr_gra_max(10) = dim_max_n(acr_gra(10,0:2), 0)           ;;4月      
      acr_gra_max(11) = dim_max_n(acr_gra(10,11:12), 0)

      spb =  dim_sum_n_Wrap(acr_gra_max, 0)
    end if
 return spb
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

dim = dimsizes(nino34_12)


;;; total year SPB
  acr_total = acr_cal(nino34_12(:,0:dim(0)-2), nino34_12(:,1:dim(0)-1))
  spb_total = spb_cal(acr_total,1)
;;;



;;; cal Every year
spb_year = new(yearnum,float)

ind_all1 = ispan(0, yearnum-1, 1)  
ind_all2 = ispan(0, yearnum-2, 1)  ;;;remain the last year 

ind_all1@_FillValue =  99
ind_all2@_FillValue =  99

spb_year@_FillValue = nino34_12@_FillValue


do i=0,yearnum-1,1

  if i.ne.yearnum-1 then 
    indx = get1Dindex_Exclude(ind_all2, i)    ;;;去掉某一年数据
  else 
    indx = get1Dindex_Exclude(ind_all1, i)    ;;;去掉某一年数据  
  end if 
  data1 = nino34_12(:,indx)
  data2 = nino34_12(:,indx+1)
  acr = acr_cal(data1,data2)
  spb_year(i) = spb_cal(acr,1)
  spb_year(i) = (/spb_total - spb_year(i)/)
  delete(indx)
  delete(acr)
  delete(data1)
  delete(data2)
end do



spb_year = dim_standardize(spb_year, 0)

opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual_1980-2019_method1.ascii"    
asciiwrite(opt1, sprintf("%6.2f,",spb_year))  
exit()

;;; plot 
;;------------------------
;;------------------------
  res = True
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@tmBorderThicknessF = 2.0

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

  res@trYMinF        = -5.0
  res@trYMaxF        = 5.0


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
  pltName = "/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/test"
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


   

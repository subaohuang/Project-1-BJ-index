;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;;
;;*************************************************************
begin
    startmon = 193101
    endmon = 194012

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)


    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    delete(date)
    delete(timeind)
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
;;;

;;;split to 12 month nino34_12 is a array size (12,yearnum)
  yearnum = dimsizes(nino34)/12
  nino34_12 = new((/12,yearnum/),float)

  do i=0,11
    nino34_12(i,:) = nino34(i:dimsizes(nino34)-1:12)
  end do 

 printVarSummary(nino34_12)


;;; cal ACF, want to get a two dimensional array 12×13
     acr = new((/12,13/),float)
     ; ; do i=0,11
     ; ;   acr(i,:) = esacr(nino34(i:dimt-1-12+i),12)
     ; ; end do 
     ;  do i=0,11
     ;     do j=0,11
     ;      acr(i,j) = escorc(nino34(i:dimt-1-24+i),nino34(i+j:dimt-1-24+i+j))
     ;     end do 
     ; end do 
     
     do i=0,11
         do j=0,12
            k=i+j 
            if k.le.11 then 
              acr(i,j) = escorc(nino34_12(i,0:yearnum-1),nino34_12(i+j,0:yearnum-1))
            else 
              acr(i,j) = escorc(nino34_12(i,0:yearnum-2),nino34_12(i+j-12,1:yearnum-1))
            end if 
         end do 
     end do 

     printVarSummary(acr)
;;;

;;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   acr_gra = -1*acr_gra
   printVarSummary(acr_gra)

   acr_gra_max = dim_max_n(acr_gra, 1) ;;因为是负的
   printVarSummary(acr_gra_max)
   
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"
   
   maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   print(maind)

   SPB = dim_sum_n_Wrap(acr_gra_max, 0) 
   
   print(SPB)

;;;;


wks = gsn_open_wks("pdf", "acr_hov-1936")
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

  
   plot = gsn_csm_contour(wks,acr,res)

   resp = True
   resp@gsMarkerColor = "black"
   resp@gsMarkerSizeF = 0.1
   resp@gsMarkerIndex = 1
   plot_gs = gsn_add_polymarker(wks, plot, maind, ispan(0,11,1), resp)

   draw(plot)
   frame(wks)

   end 

   resp@gsMarkerIndex = 1
   plot_gs = gsn_add_polymarker(wks, plot, maind, ispan(0,11,1), resp)

   draw(plot)
   frame(wks)

   end 

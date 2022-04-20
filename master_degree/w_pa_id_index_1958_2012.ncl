
begin
diri="/mnt/f/share/era/data/"
paindex=asciiread("/mnt/f/share/era/data/pa_index_1958_2012.ascii", 660, "float")
idindex=asciiread("/mnt/f/share/era/data/id_index_1958_2012.ascii", 660, "float")

; paindex=asciiread("/mnt/f/share/era/data/pa_index_1950_2010.ascii", 732, "float")
; idindex=asciiread("/mnt/f/share/era/data/id_index_1950_2010.ascii", 732, "float")

;;;处理成年数据
    paindex_anu=month_to_annual(paindex, 1)
    paindex_anu_std=paindex_anu
    paindex_anu_std=dim_standardize(paindex_anu, 1)
    paindex_anu_std!0="time" 
    paindex_anu_std&time=ispan(1, 55, 1)


    idindex_anu=month_to_annual(idindex, 1)
    idindex_anu_std=dim_standardize(idindex_anu, 1)
    idindex_anu_std!0="time" 
    idindex_anu_std&time=ispan(1, 55, 1)
;;;
 
   
; ;;;读入nino3.4 index
;     nino3_4=asciiread(diri+"Hadley_nino3_4_1950_2010_nodtrend.ascii", 732, "float")
    
;     ; nino3_4=(/dtrend(nino3_4, True)/)

;     nino3_4_anu=month_to_annual(nino3_4, 1)
;     nino3_4_anu_std=dim_standardize(nino3_4_anu, 1)
;     nino3_4_anu_std!0="time" 
;     nino3_4_anu_std&time=ispan(1, 55, 1)

;;;读入pdo index
  
  fils1=systemfunc("ls "+diri+"pdo.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->TSAXIS, 1)
  timeind=ind(date.ge.195801.and.date.le.201012)
  pdo=f1->PDO(timeind)
  pdo!0="time"
  printVarSummary(pdo)
  pdo_anu=month_to_annual(pdo, 1)
  pdo_anu&year=ispan(1, 53, 1)

;;;画趋势 pacific
    paindex_anu_std_rc=regCoef(paindex_anu_std&time,paindex_anu_std)
    pa_rc =paindex_anu_std_rc*paindex_anu_std&time+paindex_anu_std_rc@yintercept

    paindex_anu_std_rc1=regCoef(ispan(1,23,1),paindex_anu_std(0:22))
    pa_rc1 =paindex_anu_std_rc1*ispan(1,23,1)+paindex_anu_std_rc1@yintercept 

    paindex_anu_std_rc2=regCoef(ispan(24,53,1),paindex_anu_std(23:52))
    pa_rc2=paindex_anu_std_rc2*ispan(24, 53, 1)+paindex_anu_std_rc2@yintercept   
;;;

;;;indian
    idindex_anu_std_rc=regCoef(idindex_anu_std&time,idindex_anu_std)
    id_rc =idindex_anu_std_rc*idindex_anu_std&time+idindex_anu_std_rc@yintercept

    idindex_anu_std_rc1=regCoef(ispan(1,23,1),idindex_anu_std(0:22))
    id_rc1 =idindex_anu_std_rc1*ispan(1, 23, 1)+idindex_anu_std_rc1@yintercept 

    idindex_anu_std_rc2=regCoef(ispan(24,53,1),idindex_anu_std(23:52))
    id_rc2=idindex_anu_std_rc2*ispan(24, 53, 1)+idindex_anu_std_rc2@yintercept   
;;;

;;;画图
;;;
   wks=gsn_open_wks("pdf", ""+diri+"w定义的太平洋位置指数_1958_2012")
   wks1=gsn_open_wks("pdf", ""+diri+"w定义的印度洋位置指数_1958_2012")
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    ; resxy@gsnLeftString="(c)"
    resxy@gsnXYBarChart = True            ; Create bar plot
    ; resxy@tmYLMode="Explicit"
    ; resxy@trYMaxF=2.0
    ; resxy@trYMinF=-2.0
    
    ; resxy@tmYLValues=fspan(-2, 2, 9)
    ; resxy@tmYLLabels=fspan(-2, 2, 9)
    resxy@gsnYRefLine  = 0.   
    resxy@gsnBelowYRefLineBarColors="blue"
    resxy@gsnAboveYRefLineBarColors="red"
    resxy@tiYAxisString="index"
    resxy@tmYLMinorOn=False
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(paindex_anu_std&time)
    resxy@tmXBValues=ispan(1,53,5)
    resxy@tmXBLabels=ispan(1958,2010,5)
    resxy@tmXBMinorValues=ispan(1,53,1)
    
    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9  
    plot = gsn_csm_xy (wks,paindex_anu_std&time , paindex_anu_std,resxy)

    plot1 = gsn_csm_xy (wks1,idindex_anu_std&time , idindex_anu_std,resxy)

;;;将nino3_4指数画上去
 
   res=True
   res@gsLineColor="black"
   res@gsLineThicknessF=2.0
   ; Polyline=gsn_add_polyline(wks, plot, paindex_anu_std&time, nino3_4_anu_std, res)
   
   ; Polyline0=gsn_add_polyline(wks1, plot1, idindex_anu_std&time, nino3_4_anu_std, res)

    res@gsLineColor="red"
   Polyline_0=gsn_add_polyline(wks, plot, paindex_anu_std&time, pdo_anu, res)
   
   Polyline_1=gsn_add_polyline(wks1, plot1, idindex_anu_std&time, pdo_anu, res)


   res@gsLineColor="blue"
   Polyline1=gsn_add_polyline(wks, plot, paindex_anu_std&time, pa_rc, res)
    res@gsLineColor="green"
   Polyline2=gsn_add_polyline(wks, plot, ispan(1, 23, 1), pa_rc1, res)
    res@gsLineColor="green"
   Polyline3=gsn_add_polyline(wks, plot, ispan(24, 53, 1), pa_rc2, res)
  
   
   res@gsLineColor="blue"
   Polyline4=gsn_add_polyline(wks1, plot1, idindex_anu_std&time, id_rc, res)
    res@gsLineColor="green"
   Polyline5=gsn_add_polyline(wks1, plot1, ispan(1, 23, 1), id_rc1, res)
    res@gsLineColor="green"
   Polyline6=gsn_add_polyline(wks1, plot1, ispan(24,53, 1), id_rc2, res)
  

   draw(plot)
   frame(wks)


   draw(plot1)
   frame(wks1)

   end
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

return SPB
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 190001
    endmon = 201912
    startyear = 1900
    endyear = 2019
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

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
;;;


;;;读入ERA wind data

   ; fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
   fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
   f1 = addfile(fils1, "r")
   time1 = f1->time
   date1=cd_calendar(time1, 1)

   ; endmon = endmon - 200
   endmon = 201212
   timeind1=ind(date1.ge.startmon.and.date1.le.endmon)
 
  
;;;;U,V
    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_200  = f1->uwnd(timeind1,{200},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

   ; fils2 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
   ; f2 = addfile(fils2, "r")
   ; time2 = f2->time
   ; date2 = cd_calendar(time2, 1)
   ; timeind2 = ind(date2.ge.197901.and.date2.le.201912)

;;;;U,V
    ; v_200  = short2flt(f2->v(timeind2,{200},:,:))
    ; v_850  = short2flt(f2->v(timeind2,{850},:,:))
    ; v_850 = (/rmMonAnnCycTLL(v_850)/)
    ; v_200 = (/rmMonAnnCycTLL(v_200)/)

     u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
;;;




;;;ISM index  850hPA ZONAL wind diff 5-15N 40-80E - 20-30N,70-90E
    ISM = u_850(:,1,1)

    ISM = (/dim_avg_n_Wrap(u_850(:,{5:15},{40:80}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(:,{20:30},{70:90}), (/1,2/))/)

    ISM_JJA = month_to_season(ISM, "JJA")  ;;;jja 少一年，从25601开始
;;;;


;;;strong and weak monsoon year 
   total_year = ispan(startyear,endyear,1)
   std = dim_stddev_Wrap(ISM_JJA)

   inds = ind(ISM_JJA.ge.std)
   indw = ind(ISM_JJA.le.-1*std)

   indnor = ind(ISM_JJA.gt.-0.5*std.and.ISM_JJA.lt.0.5*std)

   monstryear = total_year(inds)
   monwekyear = total_year(indw)
;;;;
print(monstryear)
print(monwekyear)


;;;;;  
  


;;;split to 12 month nino34_12 is a array size (12,yearnum) monsoon year 
  yearnum = dimsizes(nino34)/12 -1
  nino34_12 = new((/12,yearnum/),float)

  do i=0,11   
    if i.le.6 then 
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
     else
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
    end if
  end do 

 printVarSummary(nino34_12)
;;;


;;; strong and weak monsoon year ssta 

  nino34_12_strmon = nino34_12(:,inds)
  nino34_12_weakmon = nino34_12(:,indw)
  nino34_12_normon = nino34_12(:,indnor)

  nino34_12_strmon2 = nino34_12(:,inds+1)
  nino34_12_weakmon2 = nino34_12(:,indw+1)
  nino34_12_normon2 = nino34_12(:,indnor+1)

  spb_strmon = spb_cal(nino34_12_strmon,nino34_12_strmon2)
  spb_weakmon = spb_cal(nino34_12_weakmon,nino34_12_weakmon2)
  spb_normon = spb_cal(nino34_12_normon,nino34_12_normon2)
;;;;

print("spb_strmon="+spb_strmon)
print("spb_weakmon="+spb_weakmon)
print("spb_normon="+spb_normon)


end 


   

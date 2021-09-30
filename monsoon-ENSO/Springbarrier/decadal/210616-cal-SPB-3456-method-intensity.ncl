;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
function spb_cal(windows,data)  ;; windows: 10 years data:two dimension 12×yearnum array i.e.12×10
begin 
  yearnum = windows
  ;;; cal ACF, want to get a two dimensional array 12×13
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
 ;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive

   acr_gra_max = acr(:,1)

   

   acr_gra_max(0) = dim_max_n(acr_gra(0,2:5),0)
   acr_gra_max(1) = dim_max_n(acr_gra(1,1:4),0)
   acr_gra_max(2) = dim_max_n(acr_gra(2,0:3),0)
   acr_gra_max(3) = dim_max_n(acr_gra(3,0:2),0)
   acr_gra_max(4) = dim_max_n(acr_gra(4,:),0)
   acr_gra_max(5) = dim_max_n(acr_gra(5,:),0)

   do i =6,11
     ; acr_gra(i) = acr(i,11-i+2) - acr(i,11-i+5)
     acr_gra_max(i) = dim_max_n(acr_gra(i,12-i+3:12-i+6),0)
   end do 

   SPB = dim_avg_n_Wrap(acr_gra_max, 0) 
return SPB
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 189901
    endmon = 201912
    startyear = 1899
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


;;;split to 12 month nino34_12 is a array size (12,yearnum)
  yearnum = dimsizes(nino34)/12
  nino34_12 = new((/12,yearnum/),float)

  do i=0,11
    nino34_12(i,:) = nino34(i:dimsizes(nino34)-1:12)
  end do 

 printVarSummary(nino34_12)
;;;


;;;get 10years running windows SPB strength 

    length = yearnum
    nave = 10
    SPB_series = new((/length-nave+1/),float)

  do t = nave,length  
    SPB  = spb_cal(nave, nino34_12(:,t-nave:t-1))
    SPB_series(t-nave) = SPB 
    delete(SPB)
  end do 

 print(SPB_series)

 SPB_series!0 = "start_time_tenyears"
 SPB_series&start_time_tenyears = ispan(startyear,endyear-nave+1,1)

 printVarSummary(SPB_series)
 
 system("rm SPB_series_1899-2019-rmtrend-3456-method.nc")
 f2 = addfile("SPB_series_1899-2019-rmtrend-3456-method.nc", "c")

 ; system("rm SPB_series_1979-2019.nc")
 ; f2 = addfile("SPB_series_1979-2019.nc", "c")
 f2->SPB_series =  SPB_series
end 


   
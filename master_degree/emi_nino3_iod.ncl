;******************************************************************
 ;IOD & IOD modoki sst
;***************************************************************
;;选择月份的function
function iod_select(iod_year,month)                                                          ;
begin
 iod_month=new((/dimsizes(iod_year)/),string)
   if month.le.9 then
  do i=0,dimsizes(iod_year)-1
  iod_month(i)=""+iod_year(i)+""+0+""+month+""
  end do
   else
  do i=0,dimsizes(iod_year)-1
  iod_month(i)=""+iod_year(i)+""+""+month+""
  end do
  end if
  iod_month_new= stringtointeger(iod_month)
return iod_month_new
end 
;;----

begin
;;;读入
	diri="/mnt/f/share/ersst/"
	diri_pi="/mnt/f/share/era/"
	f0=addfile("/mnt/f/share/ersst/HadISST_sst_chazhi_2012.nc", "r")
	date=f0->time
	lon=f0->longitude
	lat=f0->latitude

	time=cd_calendar(date, 1)
	timeind=ind(time.le.201212.and.time.ge.195001)	
    time_new=time(timeind)
	sst=f0->sst_int(:,:,timeind)
    
 printVarSummary(sst)
 ssta=dtrend_msg_n(ispan(0,dimsizes(sst&time)-1,1),rmMonAnnCycLLT(sst(latitude|:,longitude|:,time|:)),\
            False,False,2)


 copy_VarCoords(sst, ssta)
 ssta=(/runave_n_Wrap(ssta, 3, -1, 2 )/)
;;;

;;;定义nino3 150w-90w , 5n-5s
	
	nino3=dim_avg_n_Wrap(ssta({-5:5},{210:270},:), (/0,1/))
    printVarSummary(nino3)
;;;

;;;定义EMI (ssta)(165E-140w,10s-10N) - 0.5*ssta(110w-70w,15s-5n)-0.5*ssta(125E-145E,10s-20N)


   area1=dim_avg_n_Wrap(ssta({-10:10},{165:220},:), (/0,1/))                 
   area2=dim_avg_n_Wrap(ssta({-15:5},{250:290},:), (/0,1/))
   area3=dim_avg_n_Wrap(ssta({-10:20},{125:145},:), (/0,1/))

   EMI=area1-0.5*area2-0.5*area3
;;;
 delete(area1)
 delete(area2)
 delete(area3)

;;;iod 50E-70E ,10s-10n  减去 90E-110E 10s-0N
    
   area1=dim_avg_n_Wrap(ssta({-10:10},{50:70},:), (/0,1/))                 
   area2=dim_avg_n_Wrap(ssta({-10:0},{90:110},:), (/0,1/))
   iod=area1-area2
;;;
 printVarSummary(iod)

;;;对各指数进行标准化处理
   nino3_std=dim_standardize_n_Wrap(nino3, 1, 0)
   EMI_std=dim_standardize_n_Wrap(nino3, 1, 0)
   iod_std=dim_standardize_n_Wrap(nino3, 1, 0)
;;; 挑选IOD事件 传统IOD：1961，1972,1982,1987,1997,2011,2012 正事件 ; IOD modoki :1963,1967,1976,1991,1994,2003,2008
  ;;;取JASON
    iod_year=(/1961,1972,1982,1987,1997,2011,2012/)
    iod_modoki_year=(/1963,1967,1976,1991,1994,2003,2008/)
    month=ispan(9,12,1)   ;;;取9,10,11,12
    iod_yr_mon=new((/dimsizes(iod_year),dimsizes(month)/),string)
    iod_modoki_yr_mon=new((/dimsizes(iod_year),dimsizes(month)/),string)

    do i=0,dimsizes(iod_year)-1
  do j=0,dimsizes(month)-1
    if (j+9).le.9 then 
     iod_yr_mon(i,j)=""+iod_year(i)+""+0+""+month(j)+""
     iod_modoki_yr_mon(i,j)=""+iod_modoki_year(i)+""+0+""+month(j)+""
   else 
     iod_yr_mon(i,j)=""+iod_year(i)+""+""+month(j)+""
     iod_modoki_yr_mon(i,j)=""+iod_modoki_year(i)+""+month(j)+""
    end if
   end do
  end do
;;;

;;;取iod和iodmodoki事件发生年份时，JASON，各太平洋的指数情况，取它们标准化的平均
   nino3_iod=new((/dimsizes(iod_year)/),typeof(nino3_std))
   EMI_iod=nino3_iod
 
   do i=0,dimsizes(iod_year)-1
     nino3_iod(i)=dim_avg_n(nino3_std(get1Dindex(time_new,iod_yr_mon(i,:))),0)
     EMI_iod(i)=dim_avg_n(EMI_std(get1Dindex(time_new,iod_yr_mon(i,:))),0)
   end do

   nino3_iod_modoki  =new((/dimsizes(iod_modoki_year)/),typeof(nino3_std))
   EMI_iod_modoki=nino3_iod_modoki
 
   do i=0,dimsizes(iod_modoki_year)-1
     nino3_iod_modoki(i)=dim_avg_n(nino3_std(get1Dindex(time_new,iod_modoki_yr_mon(i,:))),0)
     EMI_iod_modoki(i)=dim_avg_n(EMI_std(get1Dindex(time_new,iod_modoki_yr_mon(i,:))),0)
   end do
   print(nino3_iod)


end
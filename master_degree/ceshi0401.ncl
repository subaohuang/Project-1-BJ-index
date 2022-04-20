begin
model="CMCC-CM"
diri="/mnt/f/share/CMIP5/msf_history/"
fils1=systemfunc("ls "+diri+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")
 f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    ; yearstart=197901
    ; yearend=200512  
    ; year1=1979
    ; year2=2005
    ; timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date)
    end
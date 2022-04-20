;;;插值速度势资料到等压面上
;;;
begin
diri="/mnt/f/share/ncep/ncep1/"
fils=systemfunc("ls "+diri+"chi.mon.mean.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
print(date)
end
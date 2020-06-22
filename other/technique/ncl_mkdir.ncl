 begin
;;读入数据
  ; model="bcc-csm1-1"
   model="ACCESS1-0"
   DIR = " /mnt/f/share/CMIP5/msf_rcp45/ACCESS1-0"
   system("if ! test -d " + DIR +" ; then mkdir " + DIR + " ; fi")
   end 
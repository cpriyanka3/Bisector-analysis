pro call_bisector,starname

common_string='_G2_ord0to68' ; change common_string for mask and order selection

RV_path='/Users/priyanka/newpipeline/PARAS_data/RV/'+starname+'/'
RV_file=RV_path+starname+strmid(common_string,3,11)+'_frozen_RV.tbl'

readcol,RV_file,FILENAME1, HJD, RV_NO_DRIFTCORR, DRIFT, RV_FINAL, PHOTON_NOISE,CCF_FWHM,S2N_NEAR_550nm,PIPELINE_VERSION,WAVECAL_VERSION, ANALYSIS_VERSION, format=('a,d,d,d,d,d,d,d,d,d,d'),/silent


part1=strmid(filename1,0,56)
part2='/ccf/'
part3=strmid(filename1,62,5)
part4='_'+starname

part5='_ccf_ascii.tbl'

filename=part1+part2+part3+part4+common_string+part5

RV_final_rel=(RV_final-mean(RV_final))*1000.0

;filename=findfile('*'+starname+'*.tbl')

if ~file_test(starname) then file_mkdir,starname

openw,lun,starname+'/'+starname+'_bisec.op',/get_lun
printf,lun,'HJD      BIS (m/s)	  e_BIS (m/s)   RV (m/s)   e_RV (m/s) '    

n=n_elements(filename)


for i=0,n-1 do begin

readcol,filename(i),vel,ccf_data,ccf_fit,format=('d,d,d'),/silent


nel=n_elements(vel)

vel=vel[5:nel-1]
ccf_data=ccf_data[5:nel-1]
ccf_fit=ccf_fit[5:nel-1]


pix_arr=findgen(n_elements(vel))

;beg_vel=vel[0]+15.0
;end_vel=vel[n_elements(vel)-1]-10.0

beg_vel=vel[0]
end_vel=vel[n_elements(vel)-1]

;x1=closest(vel,beg_vel)
;x2=closest(vel,end_vel)

x_arr=dindgen(n_elements(ccf_data))

fit=mpfitpeak(x_arr,ccf_data,out,nterms=6,/negative)

core0=out[1]
x1=core0-4.0*out[2]
if x1 lt 0 then x1=0
x2 =core0+4.0*out[2]

norm_value=median([ccf_data[0:x1],ccf_data[x2:n_elements(vel)-1]])

ccf_data=ccf_data/norm_value


refit=mpfitpeak(x_arr[x1:x2],ccf_data[x1:x2],out,nterms=6,/negative)
contin=out[3] + (out[4]*x_arr) + out[5]

norm_value=min(contin)

ccf_data[x1:x2]=ccf_data[x1:x2]/norm_value

refit=mpfitpeak(x_arr[x1:x2],ccf_data[x1:x2],out,nterms=6,/negative)
contin=out[3] + (out[4]*x_arr) + out[5]

cont_l=median(contin[x1:x1+10])
cont_r=median(contin[x2-10:x2])

if cont_l le cont_r then cont=cont_l else cont=cont_r

prf=ccf_data

cont_x1=double(ccf_data[x1])
cont_x2=double(ccf_data[x2])

int=0.1

bisector,prf,x1,x2,cont,b_loc,b_int,vel,int,bis=bis,star=starname,name=filename(i),dv=dv,s2n=S2N_NEAR_550nm(i),/plot

printfile=strmid(filename(i),36,strlen(filename(i))-1)

printf,lun,hjd(i)+2400000.0,bis,dv,RV_final_rel(i),PHOTON_NOISE(i)*1000,format='(d14.6,3x,d8.2,3x,d8.2,3x,d8.2,3x,d8.2)'

endfor

close,lun

print,'Press .con to continue to plot the BIS vs RV plot'

stop

readcol,starname+'/'+starname+'_bisec.op',hjd,bis,e_bis,rv,e_rv,format=('d,d,d,d,d'),/silent

loadct,13,/silent
set_plot,'ps'
plotsym,0

y1=min(bis)-max(e_bis)
y2=max(bis)+max(e_bis)

device, filename=starname+'/'+starname+'.ps',/color
plot,rv,bis,xtitle='RV (m/sec)',ytitle='Bisector span velocity (m/sec)',psym=2,xthick=5,ythick=5,charthick=5,yr=[y1,y2]
oploterror,[rv],bis,[e_rv],e_bis,psym=2

device,/close
stop
end

!     3/3/2017 changed 1/1/2020
!     o:\ssmi\routines\readin_data_files_v8.f changed to o:\ssmis\routines\readin_data_files_v8.f
!     this change has no functional effect

!     this is the v8 version
!     compiled with compile64 


!     following comments pertain to v07 versiion dated 6/28/2016

!     8/28/2014 change 6/28/2016
!     iorbit added to call fd_l1b(1,iorbit)
!     call reset_iflag_cal is now done
!     following routines were updated to handling f17 37v anomaly
!     'o:\ssmis\routines\write_l1b_realtime.f'
!     'o:\ssmis\l2_processing\process_level2.f'
!    'o:\ssmis\routines\avg_cal_counts.f' 
!    'o:\ssmis\routines\ck_counts.f' 
!    'o:\ssmis\routines\ck_counts1.f' 
!    'o:\ssmis\routines\remove_moon_effect.f'
!    'o:\ssmis\routines\resample_ta.f'
!    'o:\ssmis\routines\fd_l1b.f' 


!    the 7/31/2014 version was changd on 8/28/2014 to run on Galaxy2 with l1b_inputs.txt file on Galaxy2.  Log file name and path changed.

!    the 4/21/2014 version was changed on 7/31/2014 for a new location of the l1b_inputs.txt file

!    The 9/26/2011 version was changed on 4/21/2014 for new location of where L1B files are written and new location of log files.

!    The 6/15/2011 automated code was changed on 9/26/2011 to the new code for v7 data.  for changes made, see non-automated version of this code 
!    and changes to make it run on ops1p instead of ssmiops

!     6/15/2011 the mk_batch_l1b_files.f program was changed to accept inputs from a file , d.smith  
!    the result is this program, mk_batch_l1b_files_automated.f
!    the .exe is located and run on \\ssmiops\f external drive where it is used in an automated fashion to make the quarterly files for NSIDC
!    program was tested to work on June 15, 2011 by Deborah and Scott

!     6/10/2010 version was recompiled on  7/24/2010. no changes were made, i just wanted to make sure this 
!     program was uptodate with ssmis_precessor

!     this program was updated on 6/10/2010 to correspond to the 6/10/2010 version of ssmis_processor.f
!     it should produce the same l1b files that the realtime processing produces.

     include 'o:\ssmis\routines\l2_module_ssmis_v8.f'            

     program mk_batch_l1b_files_v08_automated
     use l2_module                                            
     implicit none
 
     character(1) adrive
     character(100) filename
     character(100) infile
 
     real(8) start_time
     integer(4) iorbit,iorbit1,iorbit2,iexist,ilu
     integer(4) ibad,kbad
 
     real(8) secyr,secdy
     integer(4) lyear,idayjl,imon,idaymo
 
 c***********************************************************************************************
 CC!    Run this program on 
 !    OPEN INPUT FILE AND READ DRIVE, SATELLITE, ORBIT1 AND ORBIT2
 cold  infile='S:\SSMI\_programs\mk_batch_l1b_inputs.txt'
 cold  infile='S:\SSMI\_programs\mk_L1B\mk_batch_l1b_inputs.txt'
       infile='c:\process\ssmis_l1b\mk_l1b\mk_batch_l1b_inputs.txt'
       open(7,file=infile,status='old')
     read(7,'(a1)') adrive
     read(7,'(i2.2)') ksat
     read(7,'(i5.5)') iorbit1
     read(7,'(i5.5)') iorbit2
     close(7)
 
     call define_filenames
       call readin_data_files
 
     write(filename,9001) ksat
 cold 9001 format('s:\ssmi\_programs\mk_batch_l1b_files_v07_',i2.2,'_',i5.5,'_',i5.5,'.txt')
  9001 format('c:\process\ssmis_L1B\mk_L1B\fortran_mk_L1B_F',i2.2,'_log.txt')
     open(6,file=filename,position='append')              
      
       do iorbit=iorbit1,iorbit2
 
     kbad=0
     do ibad=1,nbad_orbits
     if(iorbit.eq.iorbit_bad(ibad)) kbad=1
     enddo
     if(kbad.eq.1) cycle
 
     call read_l1a_file(iorbit, iexist,start_time)
 
     if(iexist.eq.0) then
     write(*,1001) iorbit
     write(6,1001) iorbit                                         
  1001 format(' no l1a file for orbit ',i6)
       cycle
     endif
 
       call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
          write(*,'(i3,i7,i6,i5,2i3,f7.2)') ksat,numscan,iorbit,lyear,imon,idaymo,secdy/3600.
          write(6,'(i3,i7,i6,i5,2i3,f7.2)') ksat,numscan,iorbit,lyear,imon,idaymo,secdy/3600.
 
     if(numscan.lt.35) cycle !not enough scans to do scan averaging and resampling 
     
     call fd_l1b(1,iorbit)    !1 means do 92 ghz
     
       call resample_ta  !resample 19/22 to 37 ghz location
 
 !     all of the above routines depend on iflag_cal being defined as it has been in the past
 !     these routines use iflag_37v for the additional ac that is required for f17 37v
 !     now at this point in the process in prep for l2b outputing , iflag_cal is reset to include iflag_chn 
     call reset_iflag_cal
 
       call get_l1a_user_filename(adrive,ksat,iorbit, filename,iexist) 
 
       ilu=4
     call open_binary(ilu,filename,'new','write','denywr',60)
       call mk_l1b_file(ilu)  
     close(ilu)
 
     enddo  !iorbit                
      
     stop 'norm end'
     end                                              
      
                                                             
 
 
     include 'o:\ssmis\users\get_l1b_user_filename_v8.f' 
     
       include 'o:\ssmi\routines\filename_routines.f'         
     include 'o:\ssmi\routines\define_filenames_v8.f' 
 
     include 'o:\ssmis\routines\readin_data_files_v8.f' 
       include 'o:\algo99\mk_algo\mk_tables_v8\get_algo_delta_chi_tcos.f'  
     
     include 'o:\ssmis\routines\mk_l1b_file_v8.f'          
     include 'o:\algo99\mk_algo\mk_tables_v7\fd_22ghz_tb.f'
 
     include 'o:\ssmis\routines\read_l1a_file.f'
     include 'o:\ssmis\routines\avg_cal_counts.f' 
     include 'o:\ssmis\routines\ck_counts.f' 
     include 'o:\ssmis\routines\ck_counts1.f' 
     include 'o:\ssmis\routines\set_iflag_chn.f' 
     include 'o:\ssmis\routines\fd_arm_temp.f'
     include 'o:\ssmis\routines\remove_moon_effect.f'
     include 'o:\ssmis\routines\prep_group_v8.f' 
     include 'o:\ssmis\routines\shift_cal_counts.f' 
     include 'o:\ssmis\routines\fd_ta_adjusted_v8.f'
     include 'o:\ssmis\routines\resample_ta_v8.f'
     include 'o:\ssmis\routines\reset_iflag_cal.f' 
 
     include 'o:\ssmis\routines\fd_l1b.f'
      
     include 'o:\ssmi\routines\tatb_routines_v8.f'
     
     include 'o:\ssmi\routines\ssmi_geolocation.f'
     include 'o:\amsr_l2\v05\routines\pre_nut_routines.f'
      include 'o:\sun_moon\sun_vector.f'            
      include 'o:\sun_moon\moon_vector.f' 
 
      include 'x:\land\fdice.for'
       include 'x:\land\fdland3.f'
 
       include 'x:\mathlib\estima1x.f'
       include 'x:\mathlib\math_routines.f'
 
     include 'x:\syspgm\fd_date_2000.f'
       include 'x:\syspgm\openbig.for'
       include 'x:\syspgm\openbig_try.f'
     include 'x:\syspgm\open_file_routines.f' 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
!  This program interpolates the fields of a .h5 pipe file for new  !
!  values of my                                                     !
!                                                                   !
!  Interpolation in x: new values are filled as 0 at the end, or    !
!    final values are removed if mx is reduced.                     !
!  Interpolation in y: linear interpolation                         !
!  Interpolation in z: new values are filled as 0 in the middle, or !
!    middle values are removed if mz is reduced.                    !
!                                                                   !
!  SHC MAY/23                                                       !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main

use hdf5
use h5lt

implicit none
! Main memory
real(8), allocatable:: r_old(:),r_new(:),data_old(:,:), data_new(:,:)

! Input information
character(len=256) root,filinp,frout
integer tflag

! Previous and new field size
integer :: N_old, k_old,M_old
integer :: N_new, k_new,M_new, lv_i, lv

! Field to be written

real(8) ::  time,Re,alpha

character(len=256) frt  ! Input root and file

integer(hid_t) :: finp_id,fout_id          ! File identifier
integer(hid_t) :: dset_id,dspace_id
integer(hid_t) :: plis_id! Property list identifier

integer :: h5err, hdferr                   ! Error flag

! Dimensions of the datasets
integer(hsize_t),dimension(2):: dim_wr
integer(hsize_t),dimension(2):: dims2, dim_rd
integer(hsize_t),dimension(1):: hdims,gr_id(3)
integer(size_t), parameter :: siever = 4*1024*1024

character(len=256) text  ! To read files
character(len=6) :: dsetname='/Ur/Re' ! Uz

integer :: ii,jj,kk

double precision :: dt_old,dtcfl_old, dtcor_old, dt_new,dtcfl_new, dtcor_new

! Variables:: 


! Initialize FORTRAN interface.
call h5open_f(h5err) 

! Fields to be inbterpolated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!               PART 1                !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input file
open(19,file='hre.inter',status='old')
call avanza(text)
read(text,'(a256)') root
call avanza(text)
read(text,'(a256)') filinp
call avanza(text)
read(text,'(a256)') frout
call avanza(text)
read(text,*) N_new
call avanza(text)
read(text,*) k_new
call avanza(text)
read(text,*) M_new
call avanza(text)
read(text,*) tflag


frt = trim(root)//trim(filinp) ! Merge input root and input file


! Open input file
call h5fopen_f(frt,H5F_ACC_RDONLY_F,finp_id,h5err)



! Read the size of the four large arrays
call h5dopen_f(finp_id,dsetname, dset_id, h5err)



! Figure out the size of the array.

! Get the dataspace ID
call h5dget_space_f(dset_id,dspace_id,h5err)




! Getting dims from dataspace
call h5sget_simple_extent_dims_f(dspace_id, dims2, dim_rd, hdferr)


call h5sclose_f(dspace_id, h5err)



call h5dclose_f(dset_id, h5err)



! Create output file
call h5pcreate_f(H5P_FILE_ACCESS_F,plis_id,h5err)
call h5pset_sieve_buf_size_f(plis_id, siever, h5err)



call h5fcreate_f(frout,H5F_ACC_TRUNC_F,fout_id,h5err,H5P_DEFAULT_F,plis_id)
! call h5fcreate_f(frout,H5F_ACC_TRUNC_F,fout_id,h5err)



call h5pclose_f(plis_id,h5err)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!            HEADER VALUES            !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read datasets from input file
hdims = (/ 1 /)
if (tflag == 0) then
   call h5ltread_dataset_double_f (finp_id,"time",time,hdims,h5err)
else
   time = 0d0
endif



call h5ltread_dataset_double_f(finp_id,"Re",Re,hdims,h5err)



call h5ltread_dataset_double_f(finp_id,"alp" ,alpha,hdims,h5err)
call H5LTread_dataset_int_f(finp_id,"N",N_old,hdims,h5err)
call H5LTread_dataset_int_f(finp_id,"K",K_old,hdims,h5err)
call H5LTread_dataset_int_f(finp_id,"M",M_old,hdims,h5err)

! Temporal

call H5LTread_dataset_double_f(finp_id,"dt",dt_old,hdims,h5err)
call H5LTread_dataset_double_f(finp_id,"dtcfl",dtcfl_old,hdims,h5err)
call H5LTread_dataset_double_f(finp_id,"dtcor",dtcor_old,hdims,h5err)



! k_new = k_old
! M_new = m_old

! Necesario para que inicialice la sim bien
dt_new = dt_old
dtcfl_new = dtcfl_old
dtcor_new = dtcor_old


write(*,*) 'N, radial    ', N_old, ' ----------> ', '   N', N_new 
write(*,*) 'K, axial     ', K_old, ' ----------> ', '   K', K_new
write(*,*) 'M, azimuthal ', M_old, ' ----------> ', '   M', M_new





! Write datasets in output file

call h5ltmake_dataset_double_f(fout_id,"time",1,hdims,(/time/),h5err)
call h5ltmake_dataset_double_f(fout_id,"Re",1,hdims,(/Re/),h5err)
call h5ltmake_dataset_double_f(fout_id,"alpha",1,hdims,(/alpha/),h5err)

call h5ltmake_dataset_int_f(fout_id,"N" ,1,hdims,(/N_new/),h5err)
! call h5ltmake_dataset_int_f(fout_id,"M" ,1,hdims,(/M_old/),h5err)
! call h5ltmake_dataset_int_f(fout_id,"K" ,1,hdims,(/K_old/),h5err)

call h5ltmake_dataset_int_f(fout_id,"M" ,1,hdims,(/M_new/),h5err)
call h5ltmake_dataset_int_f(fout_id,"K" ,1,hdims,(/K_new/),h5err)

call h5ltmake_dataset_int_f(fout_id,"Mp",1,hdims,(/1/),h5err)

call h5ltmake_dataset_double_f(fout_id,"dt"   ,1,hdims,(/dt_new/),h5err)
call h5ltmake_dataset_double_f(fout_id,"dtcfl",1,hdims,(/dtcfl_new/),h5err)
call h5ltmake_dataset_double_f(fout_id,"dtcor",1,hdims,(/dtcor_new/),h5err)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!             MESH VALUES             !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read r

allocate(r_old(N_old),r_new(N_new))
hdims = (/ N_old /)
call h5ltread_dataset_double_f(finp_id,"r",r_old,hdims,h5err)
call newMesh(r_new,N_new)
hdims=(/N_new/)
call h5ltmake_dataset_double_f(fout_id,"r",1,hdims,r_new,h5err)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!                FIELDS               !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fields are read and written plane by plane

! Set global dimension

lv_i = (2*K_old-1)*(M_old-1) + K_old; !Da problemas

lv = (2*K_new-1)*(M_new-1) + K_new;
dim_wr = (/N_new,lv/)
dim_rd = dims2

allocate (data_old(dim_rd(1),dim_rd(2)))
allocate (data_new(dim_wr(1),dim_wr(2)))




call h5gcreate_f(fout_id, '/Ur', gr_id(1)   , h5err)
call h5gcreate_f(fout_id, '/Uz', gr_id(2)   , h5err)
call h5gcreate_f(fout_id, '/Ut', gr_id(3)   , h5err)




   ! Open it
   call H5LTread_dataset_double_f(finp_id,'/Ur/Re',data_old,dim_rd,h5err)
   call interpy(data_old,data_new,r_old,r_new,N_old,N_new,lv_i,lv)

               write(*,*) 'Hasta aqui bien'
            stop 
   call h5ltmake_dataset_double_f(gr_id(1),'Re',2,dim_wr,data_new,h5err)
   

   ! Open it
   call H5LTread_dataset_double_f(finp_id,'Ur/Im',data_old,dim_rd,h5err)
   call interpy(data_old,data_new,r_old,r_new,N_old,N_new,lv_i, lv)
   call h5ltmake_dataset_double_f(gr_id(1),'Im',2,dim_wr,data_new,h5err)



   ! Open it
   call H5LTread_dataset_double_f(finp_id,'/Uz/Re',data_old,dim_rd,h5err)
   call interpy(data_old,data_new,r_old,r_new,N_old,N_new,lv_i, lv)   
   call h5ltmake_dataset_double_f(gr_id(2),'Re',2,dim_wr,data_new,h5err)


   ! Open it
   call H5LTread_dataset_double_f(finp_id,'Uz/Im',data_old,dim_rd,h5err)
   call interpy(data_old,data_new,r_old,r_new,N_old,N_new,lv_i, lv)
   call h5ltmake_dataset_double_f(gr_id(2),'Im',2,dim_wr,data_new,h5err)

   ! Open it
   call H5LTread_dataset_double_f(finp_id,'/Ut/Re',data_old,dim_rd,h5err)
   call interpy(data_old,data_new,r_old,r_new,N_old,N_new,lv_i, lv)   
   call h5ltmake_dataset_double_f(gr_id(3),'Re',2,dim_wr,data_new,h5err)



   ! Open it
   call H5LTread_dataset_double_f(finp_id,'Ut/Im',data_old,dim_rd,h5err)
   call interpy(data_old,data_new,r_old,r_new,N_old,N_new,lv_i, lv)
   call h5ltmake_dataset_double_f(gr_id(3),'Im',2,dim_wr,data_new,h5err)



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Close files

call h5gclose_f(gr_id(1),h5err)
call h5gclose_f(gr_id(2),h5err)
call h5gclose_f(gr_id(3),h5err)

call h5fclose_f(finp_id,h5err)
call h5fclose_f(fout_id,h5err)

call h5close_f(h5err)




end program


subroutine interpy(y,ye,x,xe,my,mye,lv_i,lv)

! Interpolation in y. Linear interpolation of the field

implicit none

integer, intent(in) :: my,mye,lv, lv_i
real(8), intent(in) :: y(my,lv_i),x(my),xe(mye)  ! Initial field
real(8), intent(out):: ye(mye,lv)  ! Final   field (interpolated)


integer ii  ! Counter of the points in xe
integer jj  ! Counter of the points in x
integer kk


! Set jj counter at the first point of x
jj = 1
do ii = 1,mye
   ! Check if ye(ii) is between x(jj) and x(jj+1)
   do while (xe(ii) > x(jj+1))
   !   write(*,*) xe(ii),x(jj+1),ii,jj
      jj = jj + 1
   enddo


   ! Make linear interpolation
   ye(ii,:) = (y(jj+1,:)      -y(jj,:)        )/(x(jj+1)-x(jj))*xe(ii) - &
              (y(jj+1,:)*x(jj)-y(jj,:)*x(jj+1))/(x(jj+1)-x(jj))

   
enddo

endsubroutine

subroutine newMesh(r,i_N)
implicit none
integer ::n, N_, i_N
double precision :: r(i_N), dr
double precision, parameter :: d_pi = 3.141592653589793d0

N_ = i_N+int(dsqrt(dble(i_N)))
do n = N_-i_N+1, N_
   r(n-N_+i_N) = 0.5d0*( 1d0+dcos(d_PI*(N_-n)/dble(N_)) )
end do
do n = 1, 10
   dr = 1.5d0*r(1)-0.5d0*r(2)
   r = r*(1d0+dr) - dr
end do

end subroutine newMesh



subroutine avanza(text)
implicit none
character(len=256) text
do
   read(19,'(a)') text
   if(text(1:2)/='CC') exit
enddo
end subroutine

subroutine wrtfile_Non_MPIO(dat,rnorm,filenum, irec_start,ny1m,nz1m)
 
  use parameter_mod

  integer:: num_sdat

  integer*8 filenum,irec_start,iry1,iry2,irz1,irz2
  double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: dat
  real*4, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: sdat
  integer:: ip, iry, irz, i, j, k, recnum, ii
  integer*8 keg, kbg, jeg, jbg, icount,ny1m,nz1m
  double precision rnorm
!
!HXV
!
      icount=0
!

!
!  begin by converting data to REAL*4
!
   if (filenum > 1000) then
      print *, 'large file num: ', filenum
   endif

   do k = kb-1, ke+1
     do j = jb-1,je+1
       do i = 1, nxmax
         sdat(i,j,k) = rnorm*dat(i,j,k)
       enddo
     enddo
   enddo
!
!  send data to myid = 0
!HXV
!
  num_sdat = nxmax*(nylmax+2)*(nzlmax + 2)
!
  if(myid.ne.0) then
!HXV
!
    call MPI_ISEND(sdat(1,jb-1,kb-1),num_sdat,MPI_REAL4,&
         0   ,1,MPI_COMM_WORLD,req(1),IERR)
    call MPI_WAITALL(1,req,status_array,IERR)
  else
    do ip = 0, npes-1
      keg = keglobal(ip)
      kbg = kbglobal(ip)
      jeg = jeglobal(ip)
      jbg = jbglobal(ip)
!
!
      if (jbg.eq.1.and.jeg.eq.ny1m) then
        iry1 = 0        +1
       iry2 = ny1m+1    -1
        goto 10
      endif
!
!
      if(jbg.eq.1) then
	iry1 = 0        +1
	iry2 = jeg
	j    = jb - 2
!
!      else if(jeg.eq.nymax) then
!
!
!HXV
!
      else if(jeg.eq.ny1m) then
!
	iry1 = jbg
!
!	iry2 = nymax
!
!HXV
!
	iry2 = jeg+1    -1
!
	j = jb - 1
      else
	iry1 = jbg
	iry2 = jeg
	j = jb - 1
      endif
!
!
 10   continue
!
!
      if (kbg.eq.1.and.keg.eq.nz1m) then
        irz1=0          +1
        irz2=nz1m+1     -1
        goto 20
      endif
!
!
      if(kbg.eq.1) then
	irz1 = 0        +1
	irz2 = keg
	k = kb - 2
!
!      else if(keg.eq.nzmax) then
!
!
!HXV
!
      else if(keg.eq.nz1m) then
!
	irz1 = kbg
!
!	irz2 = nzmax
!
!HXV
!
	irz2 = keg+1    -1
!
	k = kb - 1
      else
        irz1 = kbg
	irz2 = keg
	k = kb - 1
      endif
 20   continue
!
!
      if(ip.ne.0) then
!       print *,'waiting for ip = ',ip
!
        call MPI_IRECV(sdat(1,jb-1,kb-1),num_sdat,MPI_REAL4,&
             IP,1,MPI_COMM_WORLD,req(1),IERR)
!
        call MPI_WAITALL(1,req,status_array,IERR)
      endif
      
      do irz = irz1, irz2
	k = k+1
	do iry = iry1, iry2
	  j = j+1
!
!HXV
!
      icount=icount+1
!
!	  recnum = irec_start + iry + irz*nymax
	  recnum = irec_start + iry -1 + (irz-1)*(nymax - 2)
!
!	  write(filenum,rec=recnum) (sdat(i,j,k),i=1,nxmax)
!
!HXV
!
	  write(filenum,rec=recnum) (sdat(i,iry-jbg+1,irz-kbg+1),i=1   +1,nxmax    -1)
!	  write(filenum,*) (sdat(i,iry-jbg+1,irz-kbg+1),i=1,nxmax)
!
        enddo
      enddo
    enddo
!
!HXV
!
    write(6,*) " number of records written = ",icount
!
  endif
  return
end subroutine wrtfile_Non_MPIO

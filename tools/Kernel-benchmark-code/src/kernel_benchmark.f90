PROGRAM kernel_benchmark
  USE cell_types,                      ONLY: cell_type,&
                                             cell_create
  USE kinds,                           ONLY: dp
  USE qmmm_gpw_forces,                 ONLY: qmmm_forces_with_gaussian_LG,  qmmm_forces_with_gaussian_LR
  
  IMPLICIT NONE
  
  ! variables that will be passed to kernel
  INTEGER                                            :: pgfs_size  
  REAL(KIND=dp), DIMENSION(3)                        :: cgrid_pw_grid_dr
  REAL(KIND=dp), DIMENSION(3)                        :: cgrid_pw_grid_npts
  REAL(KIND=dp)                                      :: cgrid_pw_grid_dvol
  INTEGER, DIMENSION(2, 3)                           :: cgrid_pw_grid_bounds
  INTEGER, DIMENSION(2, 3)                           :: cgrid_pw_grid_bounds_local
!  REAL(KIND=dp), DIMENSION(-20:20,-20:20,-20:20), TARGET     :: cgrid_pw_grid_cr3d
  REAL(KIND=dp), DIMENSION(:,:,:), POINTER           :: cgrid_pw_grid_cr3d
  INTEGER                                            :: num_mm_atoms
  REAL(KIND=dp), DIMENSION(:), POINTER               :: mm_charges 
  INTEGER, DIMENSION(:), POINTER                     :: mm_atom_index
  REAL(KIND=dp), DIMENSION(:,:), POINTER             :: mm_particles_r  
  INTEGER                                            :: para_env_num_pe
  INTEGER                                            :: para_env_mepos
  REAL(KIND=dp), DIMENSION(:, :), POINTER            :: Forces
  REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: per_pot_cr3d
  INTEGER, DIMENSION(3)                              :: per_pot_npts
  REAL(KIND=dp), DIMENSION(3)                        :: per_pot_dr
  INTEGER, DIMENSION(:), POINTER                     :: per_pot_mm_atom_index
  TYPE(cell_type), POINTER                           :: mm_cell
  REAL(KIND=dp), DIMENSION(3)                        :: dOmmOqm
  INTEGER                                            :: iw, par_scheme
  REAL(KIND=dp), DIMENSION(2)                        :: qmmm_spherical_cutoff
  LOGICAL                                            :: shells
  LOGICAL                                            :: periodic
  REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pot_pot0_2
  
  ! variables used in wrapper only
  INTEGER                                            :: mm_cell_id_nr, mm_cell_ref_count, mm_cell_symmetry_id
  LOGICAL                                            :: mm_cell_orthorhombic
  REAL(KIND=dp)                                      :: mm_cell_deth, pot_dx
  INTEGER, DIMENSION(3)                              :: mm_cell_perd
  REAL(KIND=dp), DIMENSION(3, 3)                     :: mm_cell_hmat, mm_cell_h_inv
  INTEGER                                            :: myunit
  INTEGER                                            :: i,j,k, ib1, ib2, jb1, jb2, kb1, kb2, il, iu, jl, ju
  INTEGER                                            :: imax, jmax, kmax
  INTEGER                                            :: t1, t2, t3, t4, count_rate, count_max
  
  WRITE(*,*) "Reading input files in ./data and initialising data structures needed by kernel..."
    
  CALL system_clock (t1, count_rate, count_max)
  
  ! read pgfs_size
  OPEN(newunit=myunit, file='./data/pgfs_size.dat', action='READ', status='OLD')
  READ(myunit, *) pgfs_size
  CLOSE(myunit)
  
   ! read periodic
  OPEN(newunit=myunit, file='./data/periodic.dat', action='READ', status='OLD')
  READ(myunit, *) periodic
  CLOSE(myunit)
  
  ! read cgrid
  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%dr.dat', action='READ', status='OLD')
  READ(myunit, *) cgrid_pw_grid_dr
  CLOSE(myunit)
  
  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%dvol.dat', action='READ', status='OLD')
  READ(myunit, *) cgrid_pw_grid_dvol
  CLOSE(myunit)

  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%bounds.dat', action='READ', status='OLD')
  DO i=1,2
     DO j=1,3
        READ(myunit, *) cgrid_pw_grid_bounds(i,j)
     END DO
  END DO
  CLOSE(myunit)
  
  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%bounds_local.dat', action='READ', status='OLD')
  DO i=1,2
     DO j=1,3
        READ(myunit, *) cgrid_pw_grid_bounds_local(i,j)
     END DO
  END DO
  CLOSE(myunit)

  OPEN(newunit=myunit, file='./data/cgrid%cr3d.dat', action='READ', status='OLD')
!  DO k=-20,20
!     DO j=-20,20
!        DO i=-20,20
!           READ(myunit, *) cgrid_pw_grid_cr3d (i,j,k)
!        END DO
!     END DO
!  END DO
!  CLOSE(myunit)
!  cgrid_pw_grid_cr3d_PTR => cgrid_pw_grid_cr3d

     READ(myunit, *) ib1
     READ(myunit, *) ib2
     READ(myunit, *) jb1
     READ(myunit, *) jb2
     READ(myunit, *) kb1
     READ(myunit, *) kb2
     ALLOCATE (cgrid_pw_grid_cr3d(ib1:ib2,jb1:jb2,kb1:kb2))
     DO k=kb1,kb2
        DO j=jb1,jb2
           DO i=ib1,ib2
              READ(myunit, *) cgrid_pw_grid_cr3d(i,j,k)
           END DO
        END DO
     END DO
     CLOSE(myunit)
  
  ! read num_mm_atoms
  OPEN(newunit=myunit, file='./data/num_mm_atoms.dat', action='READ', status='OLD')
  READ(myunit, *) num_mm_atoms
  CLOSE(myunit)

  ! read mm_charges
  OPEN(newunit=myunit, file='./data/mm_charges.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (mm_charges(imax))
  READ(myunit, *) mm_charges
  CLOSE(myunit)

  ! read mm_atom_index
  OPEN(newunit=myunit, file='./data/mm_atom_index.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (mm_atom_index(imax))
  READ(myunit, *) mm_atom_index
  CLOSE(myunit)
  
  ! read mm_particles
  OPEN(newunit=myunit, file='./data/mm_particles%r.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (mm_particles_r(imax,3))
  DO i=1,imax
     READ(myunit, *) mm_particles_r(i,:)
  END DO
  CLOSE(myunit)
  
  ! read para_env
  OPEN(newunit=myunit, file='./data/para_env.dat', status='OLD', action='READ')
  READ(myunit, *) para_env_num_pe
  READ(myunit, *) para_env_mepos
  CLOSE(myunit)
  
  ! read Forces
  OPEN(newunit=myunit, file='./data/Forces.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  READ(myunit, *) jmax
  ALLOCATE (Forces(imax,jmax))
  DO j=1,jmax
     DO i=1,imax
        READ(myunit, *) Forces(i,j)
     END DO
  END DO
  CLOSE(myunit)


  ! read per_pot

  IF (periodic) THEN

     OPEN(newunit=myunit, file='./data/per_pot%TabLR%cr3d.dat', status='OLD', action='READ')
     READ(myunit, *) imax
     READ(myunit, *) jmax
     READ(myunit, *) kmax
     ALLOCATE (per_pot_cr3d(imax,jmax,kmax))
     DO k=1,kmax
        DO j=1,jmax
           DO i=1,imax
              READ(myunit, *) per_pot_cr3d(i,j,k) 
           END DO
        END DO
     END DO
     CLOSE(myunit)



     OPEN(newunit=myunit, file='./data/per_pot%TabLR%pw_grid%npts.dat', status='OLD', action='READ')
     READ(myunit, *) per_pot_npts
     CLOSE(myunit)
  
     OPEN(newunit=myunit, file='./data/per_pot%TabLR%pw_grid%dr.dat', status='OLD', action='READ')
     READ(myunit, *) per_pot_dr
     CLOSE(myunit)

     OPEN(newunit=myunit, file='./data/per_pot%mm_atom_index.dat', status='OLD', action='READ')
     READ(myunit, *) imax
     ALLOCATE (per_pot_mm_atom_index(imax))
     READ(myunit, *) per_pot_mm_atom_index
     CLOSE(myunit)
  ELSE
     OPEN(newunit=myunit, file='./data/pot%mm_atom_index.dat', status='OLD', action='READ')
     READ(myunit, *) imax
     ALLOCATE (per_pot_mm_atom_index(imax))
     READ(myunit, *) per_pot_mm_atom_index
     CLOSE(myunit)

     OPEN(newunit=myunit, file='./data/pot%pot0_2.dat', status='OLD', action='READ')
     READ(myunit, *) imax
     READ(myunit, *) jmax
     ALLOCATE (pot_pot0_2(imax,jmax))
        DO j=1,jmax
           DO i=1,imax
              READ(myunit, *) pot_pot0_2(i,j)
           END DO
        END DO
     CLOSE(myunit)

     OPEN(newunit=myunit, file='./data/pot%dx.dat', status='OLD', action='READ')
     READ(myunit, *) pot_dx
     CLOSE(myunit)

     OPEN(newunit=myunit, file='./data/cgrid%pw_grid%npts.dat', action='READ', status='OLD')
     READ(myunit, *) cgrid_pw_grid_npts
     CLOSE(myunit)
  
  END IF
  
  ! read mm_cell
  OPEN(newunit=myunit, file='./data/mm_cell.dat', status='OLD', action='READ')
  READ(myunit, *) mm_cell_id_nr
  READ(myunit, *) mm_cell_ref_count
  READ(myunit, *) mm_cell_symmetry_id
  READ(myunit, *) mm_cell_orthorhombic
  READ(myunit, *) mm_cell_deth
  READ(myunit, *) mm_cell_perd
  READ(myunit, *) mm_cell_hmat
  READ(myunit, *) mm_cell_h_inv
  CLOSE(myunit)
  
  CALL cell_create(mm_cell,&
       mm_cell_id_nr, mm_cell_ref_count, mm_cell_symmetry_id,&
       mm_cell_orthorhombic,&
       mm_cell_deth, mm_cell_hmat, mm_cell_h_inv, mm_cell_perd)
  
  
  ! read d0mm0qm
  OPEN(newunit=myunit, file='./data/dOmmOqm.dat', status='OLD', action='READ')
  READ(myunit, *) dOmmOqm
  CLOSE(myunit)
  
  
  ! read iw
  OPEN(newunit=myunit, file='./data/iw.dat', status='OLD', action='READ')
  READ(myunit, *) iw
  CLOSE(myunit)
  
  
  ! read par_scheme
  OPEN(newunit=myunit, file='./data/par_scheme.dat', status='OLD', action='READ')
  READ(myunit, *) par_scheme
  CLOSE(myunit)
  
  
  ! read qmmm_spherical_cutoff
  OPEN(newunit=myunit, file='./data/qmmm_spherical_cutoff.dat', status='OLD', action='READ')
  READ(myunit, *) qmmm_spherical_cutoff
  CLOSE(myunit)
  
  
  ! read shells
  OPEN(newunit=myunit, file='./data/shells.dat', status='OLD', action='READ')
  READ(myunit, *) shells
  CLOSE(myunit)


  CALL system_clock (t2, count_rate, count_max)
  WRITE(*,'(A23,F6.2,A8)') "Done initialising, took", REAL(t2 - t1)/REAL(count_rate), " seconds"
  WRITE(*, *) "Calling kernel to compute QM/MM forces..."
  
 ! call kernel    
 IF (periodic) THEN  
 CALL qmmm_forces_with_gaussian_LG  (pgfs_size,&
      cgrid_pw_grid_dr, cgrid_pw_grid_dvol, cgrid_pw_grid_bounds, cgrid_pw_grid_bounds_local, &
      cgrid_pw_grid_cr3d,&
      num_mm_atoms,&
      mm_charges,& 
      mm_atom_index,&
      mm_particles_r,&
      para_env_num_pe, para_env_mepos,&
      Forces,&
      per_pot_cr3d, per_pot_npts, per_pot_dr, per_pot_mm_atom_index,&
      mm_cell,&
      dOmmOqm,&
      iw,&
      par_scheme,&
      qmmm_spherical_cutoff,&
      shells)
  ELSE
      CALL qmmm_forces_with_gaussian_LR  (pgfs_size,&
      cgrid_pw_grid_dr, cgrid_pw_grid_npts, cgrid_pw_grid_dvol, cgrid_pw_grid_bounds, cgrid_pw_grid_bounds_local, &
      cgrid_pw_grid_cr3d,&
      num_mm_atoms,&
      mm_charges,&
      mm_atom_index,&
      mm_particles_r,&
      para_env_num_pe, para_env_mepos,&
      Forces,&
      per_pot_mm_atom_index, pot_dx, pot_pot0_2,&
      mm_cell,&
      dOmmOqm,&
      iw,&
      par_scheme,&
      qmmm_spherical_cutoff,&
      shells)
 END IF
 
 CALL system_clock (t3, count_rate, count_max)
 WRITE(*,'(A26,F6.2,A8)') "Done running kernel, took", REAL(t3 - t2)/REAL(count_rate), " seconds"
 WRITE(*, *) "Writing resulting forces.."
 
 ! write updated forces as a check to compare to reference values
 OPEN(newunit=myunit, file='Forces.out', status='replace', action='WRITE')
 DO j = 1, SIZE(Forces,2)
    DO i = 1, SIZE(Forces,1)
       WRITE(myunit, '(I8, I8, ES14.6)') i, j, Forces(i,j)
    END DO
 END DO
 CLOSE(myunit)

 CALL system_clock (t4, count_rate, count_max)
 WRITE(*,'(A25,F6.2,A8)') "Done writing forces, took", REAL(t4 - t3)/REAL(count_rate), " seconds"
 WRITE(*,*) "Kernel benchmark completed"
 
 
END PROGRAM

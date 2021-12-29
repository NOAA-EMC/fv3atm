module module_cap_cpl
!
!*** this module contains the debug subroutines for fv3 coupled run
!
! revision history
!  12 Mar 2018: J. Wang       Pull coupled subroutines from fv3_cap.F90 to this module
!
  use ESMF

  implicit none

  private
  public diagnose_cplFields
!
  contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

    subroutine diagnose_cplFields(gcomp, clock_fv3, fcstpe, &
                                  statewrite_flag, stdiagnose_flag, state_tag)

      type(ESMF_GridComp), intent(in)       :: gcomp
      type(ESMF_Clock),intent(in)           :: clock_fv3
      logical, intent(in)                   :: fcstpe
      logical, intent(in)                   :: statewrite_flag
      integer, intent(in)                   :: stdiagnose_flag
      character(len=*),         intent(in)  :: state_tag                        !< Import or export.

      character(len=*),parameter :: subname='(module_cap_cpl:diagnose_cplFields)'
      type(ESMF_Time) :: currTime
      type(ESMF_State) :: state
      character(len=240) :: timestr
      integer :: timeslice = 1
      character(len=160) :: nuopcMsg
      character(len=160) :: filename
      integer :: rc
!
      call ESMF_ClockGet(clock_fv3, currTime=currTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_TimeGet(currTime, timestring=timestr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ClockPrint(clock_fv3, options="currTime", preString="current time: ", unit=nuopcMsg)
      call ESMF_LogWrite(trim(subname)//' '//trim(state_tag)//' '//trim(nuopcMsg), ESMF_LOGMSG_INFO)

      if(trim(state_tag) .eq. 'import')then
        call ESMF_GridCompGet(gcomp, importState=state, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(stdiagnose_flag > 0 .and. fcstpe)then
         call state_diagnose(state, ':IS', rc=rc)
        end if

        ! Dump Fields out
        if (statewrite_flag) then
          write(filename,'(A)') 'fv3_cap_import_'//trim(timestr)//'_'
          call State_RWFields_tiles(state,trim(filename), timeslice, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if
      end if

      if(trim(state_tag) .eq. 'export')then
        call ESMF_GridCompGet(gcomp, exportState=state, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(stdiagnose_flag > 0 .and. fcstpe)then
         call state_diagnose(state, ':ES', rc=rc)
        end if

        ! Dump Fields out
        if (statewrite_flag) then
          write(filename,'(A)') 'fv3_cap_export_'//trim(timestr)//'_'
          call State_RWFields_tiles(state,trim(filename), timeslice, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if
      end if

    end subroutine diagnose_cplFields

  !-----------------------------------------------------------------------------

  ! This subroutine requires ESMFv8 - for coupled FV3
    subroutine State_RWFields_tiles(state,filename,timeslice,rc)

      type(ESMF_State), intent(in)          :: state
      character(len=*), intent(in)          :: fileName
      integer, intent(in)                   :: timeslice
      integer, intent(out)                  :: rc

      ! local
      type(ESMF_Field)                       :: firstESMFFLD
      type(ESMF_Field),allocatable           :: flds(:)
      type(ESMF_GridComp) :: IOComp
      type(ESMF_Grid) :: gridFv3

      character(len=256) :: msgString
      integer                                :: i, icount, ifld
      integer                                :: fieldcount, firstfld
      character(64), allocatable             :: itemNameList(:), fldNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: typeList(:)

      character(len=*),parameter :: subname='(module_cap_cpl:State_RWFields_tiles)'

      ! local variables

      rc = ESMF_SUCCESS
      !call ESMF_LogWrite(trim(subname)//trim(filename)//": called", ESMF_LOGMSG_INFO, rc=rc)

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! find first stateitem that is a field and the count of fields
      firstfld = 0; fieldcount = 0
      do i = icount,1,-1
        if(typeList(i) == ESMF_STATEITEM_FIELD) firstfld = i
        if(typeList(i) == ESMF_STATEITEM_FIELD) fieldcount = fieldcount + 1
      enddo
      !write(msgString,*) trim(subname)//' icount = ',icount," fieldcount =
      !",fieldcount," firstfld = ",firstfld
      !call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

      allocate(flds(fieldCount),fldNameList(fieldCount))
      ifld = 1
      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          fldNameList(ifld) = itemNameList(i)
          ifld = ifld + 1
        endif
      enddo

      call ESMF_LogWrite(trim(subname)//": write "//trim(filename)//"tile1-tile6", ESMF_LOGMSG_INFO, rc=rc)
      ! get first field
      call ESMF_StateGet(state, itemName=itemNameList(firstfld), field=firstESMFFLD, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

      call ESMF_FieldGet(firstESMFFLD, grid=gridFv3, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

      IOComp = ESMFIO_Create(gridFv3, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

      do ifld=1, fieldCount
        call ESMF_StateGet(state, itemName=fldNameList(ifld), field=flds(ifld), rc=rc)
      enddo

      call ESMFIO_Write(IOComp, filename, flds, filePath='./', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

! -- Finalize ESMFIO
      deallocate(flds)
      deallocate(fldNameList)
      call ESMFIO_Destroy(IOComp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize()

      !call ESMF_LogWrite(trim(subname)//trim(filename)//": finished", ESMF_LOGMSG_INFO, rc=rc)

    end subroutine State_RWFields_tiles

  !-----------------------------------------------------------------------------

  subroutine state_diagnose(State,string, rc)
    ! ----------------------------------------------
    ! Diagnose status of state
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    character(len=*), intent(in), optional :: string
    integer, intent(out), optional  :: rc

    ! local variables
    integer                       :: i,j,n
    integer                       :: itemCount
    character(len=64) ,pointer    :: itemNameList(:)
    character(len=64)             :: lstring
    character(len=256)            :: tmpstr

    type(ESMF_Field)            :: lfield
    type(ESMF_StateItem_Flag)   :: itemType
    real(ESMF_KIND_R8), pointer :: dataPtr2d(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr3d(:,:,:)
    integer                     :: lrc, dimCount
    character(len=*),parameter  :: subname='(FV3: state_diagnose)'

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_StateGet(State, itemCount=itemCount, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(itemNameList(itemCount))

    call ESMF_StateGet(State, itemNameList=itemNameList, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    do n = 1, itemCount
       call ESMF_StateGet(State, itemName=trim(itemNameList(n)), itemType=itemType, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if(itemType == ESMF_STATEITEM_FIELD)then
       call ESMF_StateGet(State, itemName=trim(itemNameList(n)), field=lfield, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_FieldGet(lfield, dimCount=dimcount, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       if(dimcount == 2)then
         call ESMF_FieldGet(lfield, farrayPtr=dataPtr2d, rc=lrc)
         if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
           minval(dataPtr2d),maxval(dataPtr2d),sum(dataPtr2d)
         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
       else
         call ESMF_FieldGet(lfield, farrayPtr=dataPtr3d, rc=lrc)
         if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
           minval(dataPtr3d),maxval(dataPtr3d),sum(dataPtr3d)
         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
       end if
     end if
    enddo
     deallocate(itemNameList)

    if (present(rc)) rc = lrc
    call ESMF_LogWrite(subname//' exit', ESMF_LOGMSG_INFO)

  end subroutine state_diagnose

  !-----------------------------------------------------------------------------

end module module_cap_cpl

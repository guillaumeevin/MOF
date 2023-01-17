! compile with the command "R CMD SHLIB disag.f90" in a command prompt

!=========================== indexx ==========================
! Use the Heapsort algorithm to index an array arrin of length n.
! Output the array indx such that arrin(indx(j)) is in ascending
! order for j = 1,2,...,n.  The input quantities n and arrin are
! not changed.
!
! This is a Numerical Recipes routine, but modified by one
! line to work if n equals 1.
subroutine indexx(arr,n,indx)

  parameter (M=7,NSTACK=50)
  integer n,indx(n)
  double precision arr(n)
  integer i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
  double precision a

  do j=1,n
     indx(j)=j
  enddo

  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M) then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,1,-1
           if(arr(indx(i)).le.a) goto 2
           indx(i+1)=indx(i)
        enddo
        i=0
2       indx(i+1)=indxt
     enddo
     if(jstack.eq.0) return

     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2

  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp

     if(arr(indx(l+1)).gt.arr(indx(ir))) then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     endif

     if(arr(indx(l)).gt.arr(indx(ir))) then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     endif

     if(arr(indx(l+1)).gt.arr(indx(l))) then
        itemp=indx(l+1)
        indx(l+1)=indx(l)
        indx(l)=itemp
     endif

     i=l+1
     j=ir
     indxt=indx(l)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a) goto 3

4    continue
     j=j-1
     if(arr(indx(j)).gt.a) goto 4
     if(j.lt.i) goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3

5    indx(l)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     if(jstack.gt.NSTACK) return
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     endif
  endif
  goto 1

end subroutine indexx


!=========================== disagPrec_F ==========================
! For each daily simulations at n gauges, we find the closest field in observations for which
! daily precipitation structures are available.
subroutine disagPrec_F(Yobs24, Yobs1, mObs, cObs, Ysim24, mSim, cSim, nTobs,&
    &nStat, nTsim, nLagScore, Ysim1, codeDisag)
    implicit none
    !  Input/Output
    integer, intent(in) :: nTobs, nStat, nTsim, nLagScore
    integer, intent(in) :: mObs(nTobs) ! index of the month (1..12) for the matrix of 3-day observations
    integer, intent(in) :: mSim(nTsim) ! index of the month (1..12) for the matrix of 3-day simulations
    integer, intent(in) :: cObs(nTobs) ! class (=1,2,3,4) of precip for observations
    integer, intent(in) :: cSim(nTsim) ! class (=1,2,3,4) of precip for simulations
    double precision, intent(in), dimension(nTobs*24,nStat) :: Yobs1    ! Matrix of hourly observations: nTobs*24 [days] x nStat [stations]
    double precision, intent(in), dimension(nTobs,nStat) :: Yobs24  ! Matrix of observations amounts for daily periods: nTobs [days] x nStat [stations]
    double precision, intent(in), dimension(nTsim,nStat) :: Ysim24   ! Matrix of simulated amounts for daily periods: nTsim [days] x nStat [stations]
    double precision, intent(out), dimension(nTsim*24,nStat) :: Ysim1 ! Matrix of disag. simulated amounts: nTsim [days] x nStat [stations]
    double precision, intent(out), dimension(nTsim,nStat) :: codeDisag ! Matrix indicating how it has been disaggregated

    !  Locals
    integer :: nBestField,i,j,k,jRnd,iLag,ii,iBest
    double precision :: naVal,rmseIJ, r
    PARAMETER(nBestField=10,naVal=-9999.)
    double precision :: adimObs(nStat), adimSim(nStat)
    double precision, dimension(24) :: Yobs24H
    double precision, dimension(nTobs) :: rmseI
    integer, dimension(nTobs) :: indBestRMSEI, indBestrmseDayI
    integer, dimension(nBestField) :: indBestFieldI
	integer, dimension(24) :: indexi1, indexj1
    logical :: notSameClass

    ! For each simulated day
    do i = 1, nTsim
       
	   ! index hours for this day i 
	   indexi1 = (/ (ii, ii = ((i-1)*24+1), (i*24)) /)
	   
       !======= Step 1: Minimize Score for this selection
       ! sum of square differences for the closest fields on two time steps

       ! for the first time step
       if(i <= nLagScore) then
           do j = 1, nTobs
              ! same month and class
              notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
              ! if any observed value is missing in the observed prec. field or if the months do not correspond
              ! we discard this observed day
              if((any(Yobs24(j,:) == naVal)) .OR. notSameClass) then
                rmseI(j) = 1E30
              else
                ! absolute differences between adimensioned precipitation for this day
                if(sum(Ysim24(i,:))==0) then
                    adimSim = 0
                else
                    adimSim = Ysim24(i,:)/sum(Ysim24(i,:))
                end if

                 if(sum(Yobs24(j,:))==0) then
                    adimObs = 0
                else
                    adimObs = Yobs24(j,:)/sum(Yobs24(j,:))
                end if
                rmseI(j) = sum(abs(adimSim-adimObs))
              end if
           enddo
       else
         ! discard first elements
         do j = 1, nLagScore
            rmseI(j) = 1E30
         end do

         ! for the next days, compute score
         loopStationScore: do j = (nLagScore+1), nTobs
              ! same month and class
              notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
              ! if any observed value is missing in the observed prec. field or if the months do not correspond
              ! we discard this observed day
              if(any(Yobs24(j,:) == naVal) .OR. notSameClass) then
                rmseI(j) = 1E30
              else
                ! absolute differences between adimensioned precipitation for this day
                if(sum(Ysim24(i,:))==0) then
                    adimSim = 0
                else
                    adimSim = Ysim24(i,:)/sum(Ysim24(i,:))
                end if

                 if(sum(Yobs24(j,:))==0) then
                    adimObs = 0
                else
                    adimObs = Yobs24(j,:)/sum(Yobs24(j,:))
                end if
                rmseIJ = sum(abs(adimSim-adimObs))
                ! add differences for the previous days, just non na values
                do iLag = 1, nLagScore
                    if(any(Yobs24(j-iLag,:)==naVal)) then
                       rmseI(j) = 1E30
                       cycle loopStationScore
                    else
                        if(sum(Ysim24(i-iLag,:))==0) then
                            adimSim = 0
                        else
                            adimSim = Ysim24(i-iLag,:)/sum(Ysim24(i-iLag,:))
                        end if

                         if(sum(Yobs24(j-iLag,:))==0) then
                            adimObs = 0
                        else
                            adimObs = Yobs24(j-iLag,:)/sum(Yobs24(j-iLag,:))
                        end if
                        rmseIJ = rmseIJ + sum(abs(adimSim-adimObs))
                    end if
                end do
                rmseI(j) = rmseIJ
              end if
           enddo loopStationScore
       end if

       call indexx(rmseI,nTobs,indBestRMSEI)
       indBestFieldI = indBestRMSEI(1:nBestField)

       !======= Step 3: Look at the different case and disaggregate =====
       loopStationDisag: do k = 1, nStat
         ! initialise code
         codeDisag(i,k) = naVal

         !!!!! case 1: no occurrence for this period and this station, nothing to disaggregate.
         if(Ysim24(i,k)==0) then
            codeDisag(i,k) = 0
            Ysim1(indexi1,k) = 0

         !!!!! case 2: loop at the closest fields, if there is the same number of occurrences, we take the observed
         ! temporal structure (less restrictive than the exact same structure)
         else
            do j = 1, nBestField
				
              iBest = (indBestFieldI(j)-1)*24
			  indexj1 = (/ (ii, ii = (iBest+1), (iBest+24)) /)
			  Yobs24H = Yobs1(indexj1,k)
                
              ! check if there is an observed precitation for the selected observed field and for this station
              if(sum(Yobs24H)>0) then
                ! update code of disaggregation
                codeDisag(i,k) = dble(j)
                ! simulated precipitation for these 3 days are observed precipitation for the close field, rescaled
                ! by the simulated cumulated value for these 3 days
                Ysim1(indexi1,k) = Yobs24H*Ysim24(i,k)/sum(Yobs24H)
                ! codes to analyse how large daily intensities are disaggregated
                if(any(Ysim1(indexi1,k)>maxval(Yobs1(:,k)))) then
                    codeDisag(i,k) = codeDisag(i,k) + 10000.
                end if
                ! get out of the loop for this station
                cycle loopStationDisag
              end if
            end do
         end if

         !!!!! case 3: if we did not find similar structure then we find, for this station, days with similar amounts
         if(codeDisag(i,k)==naVal) then
               ! for the first time step
               if(i <= nLagScore) then
                   do j = 1, nTobs
                      ! same month and class
                      notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
                      ! if any observed value is missing in the observed prec. field or if the months do not correspond
                      ! we discard this observed day
                      if((Yobs24(j,k) == naVal) .OR. notSameClass) then
                        rmseI(j) = 1E30
                      else
                        rmseI(j) = abs(Ysim24(i,k)-Yobs24(j,k))
                      end if
                   enddo
               else
                 ! discard first elements
                 do j = 1, nLagScore
                    rmseI(j) = 1E30
                 end do

                 ! for the next days, compute score
                 do j = (nLagScore+1), nTobs
                      ! same month and class
                      notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
                      ! if any observed value is missing in the observed prec. field or if the months do not correspond
                      ! or if there is no observed precip for this day
                      ! we discard this observed day
                      if(((Yobs24(j,k) == naVal) .OR. (Yobs24(j,k) == 0)) .OR. notSameClass) then
                        rmseI(j) = 1E30
                      else
                        rmseIJ = abs(Ysim24(i,k)-Yobs24(j,k))
                        ! add differences for the previous days, just non na values
                        if(any(Yobs1(((j-nLagScore-1)*24+1):((j-1)*24),k)==naVal)) then
                            rmseIJ = 1E30
                        else
                            do iLag = 1, nLagScore
                               rmseIJ = rmseIJ + abs(Ysim24(i-iLag,k)-Yobs24(j-iLag,k))
                            end do
                        end if
                        rmseI(j) = rmseIJ
                      end if
                   enddo
               end if

            ! pick one of 100 days with similar occ structures and intensities at this station
            call indexx(rmseI,nTobs,indBestrmseDayI)
            call RANDOM_NUMBER(r)
            jRnd = indBestrmseDayI(int(r*100)+1)
            indexj1 =  (/ (ii, ii = ((jRnd-1)*24+1), (jRnd*24)) /)
            Yobs24H = Yobs1(indexj1,k)
			
			if(sum(Yobs24H)>0) then
				Ysim1(indexi1,k) = Yobs24H*Ysim24(i,k)/sum(Yobs24H)
				codeDisag(i,k) = 2000
				! codes to analyse how large 3-day intensities are disaggregated
				if(any(Ysim1(indexi1,k)>maxval(Yobs1(:,k)))) then
					codeDisag(i,k) = codeDisag(i,k) + 10000.
				end if
            else
				! last solution: uniform disaggregation
				Ysim1(indexi1,k) = Ysim24(i,k)/24
				codeDisag(i,k) = 3000.
            end if
        end if
    enddo loopStationDisag
  enddo
end subroutine disagPrec_F

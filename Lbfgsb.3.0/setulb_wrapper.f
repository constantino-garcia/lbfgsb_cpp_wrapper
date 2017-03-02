c
c Copyright Constantino Antonio Garcia 2017
c
c This Source Code Form is subject to the terms of the Mozilla Public
c License, v. 2.0. If a copy of the MPL was not distributed with this
c file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
      subroutine task_to_integer(task, itask)
        integer :: itask
        character * 60 :: task
           SELECT CASE (task)
           CASE ('START')
            itask = 0
           CASE ('NEW_X')
            itask = 1
           CASE ('FG_START')
            itask = 2
           CASE ('FG_LNSRCH')
            itask = 3
           CASE ('RESTART_FROM_LNSRCH')
            itask = 4
           CASE ('CONVERGENCE: '//
     +        'NORM_OF_PROJECTED_GRADIENT_<=_PGTOL')
            itask = 5
           CASE ('CONVERGENCE: '//
     +           'REL_REDUCTION_OF_F_<=_FACTR*EPSMCH')
            itask = 6
           CASE ('ABNORMAL_TERMINATION_IN_LNSRCH')
            itask = 7
           CASE ('ERROR: N .LE. 0')
            itask = 8
           CASE ('ERROR: M .LE. 0')
            itask = 9
           CASE ('ERROR: FACTR .LT. 0')
            itask = 10
           CASE ('ERROR: INVALID NBD')
            itask = 11
           CASE ('ERROR: NO FEASIBLE SOLUTION')
            itask = 12
           CASE DEFAULT
               itask = 13
           END SELECT
      end subroutine task_to_integer

      subroutine integer_to_task(task, itask)
        integer :: itask
        character * 60 :: task

        SELECT CASE (itask)
             CASE (0)
              task = 'START'
             CASE (1)
              task =  'NEW_X'
             CASE (2)
              task =  'FG_START'
             CASE (3)
              task =  'FG_LNSRCH'
             CASE (4)
              task =  'RESTART_FROM_LNSRCH'
             CASE (5)
              task =  'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
             CASE (6)
              task =  'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
             CASE (7)
              task =  'ABNORMAL_TERMINATION_IN_LNSRCH'
             CASE (8)
              task =  'ERROR: N .LE. 0'
             CASE (9)
              task =  'ERROR: M .LE. 0'
             CASE (10)
              task =  'ERROR: FACTR .LT. 0'
             CASE (11)
              task =  'ERROR: INVALID NBD'
             CASE (12)
              task =  'ERROR: NO FEASIBLE SOLUTION'
        END SELECT
      end subroutine integer_to_task

      subroutine csave_to_integer(csave, icsave)
        integer :: icsave
        character * 60 :: csave
        SELECT CASE(csave)
          CASE ('START')
           icsave = 0
          CASE ('FG')
           icsave = 1
          CASE ('CONVERGENCE')
           icsave = 2
          CASE ('WARNING: '//
     +           'ROUNDING ERRORS PREVENT PROGRESS')
           icsave = 3
          CASE ('WARNING: XTOL TEST SATISFIED')
           icsave = 4
          CASE ('WARNING: STP = STPMAX')
           icsave = 5
          CASE ('WARNING: STP = STPMIN')
           icsave = 6
          CASE ('ERROR: STP .LT. STPMIN')
           icsave = 7
          CASE ('ERROR: STP .GT. STPMAX')
           icsave = 8
          CASE ('ERROR: INITIAL G .GE. ZERO')
           icsave = 9
          CASE ('ERROR: FTOL .LT. ZERO')
           icsave = 10
          CASE ('ERROR: GTOL .LT. ZERO')
           icsave = 11
          CASE ('ERROR: XTOL .LT. ZERO')
           icsave = 12
          CASE ('ERROR: STPMIN .LT. ZERO')
           icsave = 13
          CASE ('ERROR: STPMAX .LT. STPMIN')
           icsave = 14
          CASE DEFAULT
c             This case should never happen
c             checked later in cpp using assert
             icsave = 15
          END SELECT
      end subroutine csave_to_integer

      subroutine integer_to_csave(csave, icsave)
        integer :: icsave
        character * 60 :: csave
                  SELECT CASE (icsave)
                    CASE (0)
                     csave = 'START'
                    CASE (1)
                     csave = 'FG'
                    CASE (2)
                     csave = 'CONVERGENCE'
                    CASE (3)
                     csave = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
                    CASE (4)
                     csave = 'WARNING: XTOL TEST SATISFIED'
                    CASE (5)
                     csave = 'WARNING: STP = STPMAX'
                    CASE (6)
                     csave = 'WARNING: STP = STPMIN'
                    CASE (7)
                     csave = 'ERROR: STP .LT. STPMIN'
                    CASE (8)
                     csave = 'ERROR: STP .GT. STPMAX'
                    CASE (9)
                     csave = 'ERROR: INITIAL G .GE. ZERO'
                    CASE (10)
                     csave = 'ERROR: FTOL .LT. ZERO'
                    CASE (11)
                     csave = 'ERROR: GTOL .LT. ZERO'
                    CASE (12)
                     csave = 'ERROR: XTOL .LT. ZERO'
                    CASE (13)
                     csave = 'ERROR: STPMIN .LT. ZERO'
                    CASE (14)
                     csave = 'ERROR: STPMAX .LT. STPMIN'
                    END SELECT
      end subroutine integer_to_csave

      subroutine setulb_wrapper(n, m, x, l, u, nbd, f, g, factr, pgtol,
     +                         wa, iwa, itask, iprint, icsave,
     +                         lsave0, lsave1, lsave2, lsave3,
     +                         isave, dsave) bind(c)
          use iso_c_binding
          integer(c_int) :: n, m, nbd(n), iwa(3 * n), iprint, isave(44),
     +      itask, icsave
          real(c_double) :: x(n), l(n), u(n), f, g(n), factr, pgtol,
     +                      wa(2 * m * n + 5 * n + 11 * m * m + 8 * m),
     +                      dsave(29)
          character * 60 :: task, csave
          logical(c_bool) :: lsave0, lsave1, lsave2, lsave3
          logical lsave(4)

          lsave(1) = lsave0
          lsave(2) = lsave1
          lsave(3) = lsave2
          lsave(4) = lsave3
          call integer_to_task(task, itask)
          call integer_to_csave(csave, icsave)

          call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol,
     +           wa, iwa, task, iprint, csave,
     +           lsave, isave, dsave)

          lsave0 = lsave(1)
          lsave1 = lsave(2)
          lsave2 = lsave(3)
          lsave3 = lsave(4)

          call task_to_integer(task, itask)
          call csave_to_integer(csave, icsave)


      end subroutine setulb_wrapper


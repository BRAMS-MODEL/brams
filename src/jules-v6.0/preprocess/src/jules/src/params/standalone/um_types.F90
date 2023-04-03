! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! JULES version of UM module um_types

MODULE um_types

USE precision_mod, ONLY: real32, real64,                                      &
                          integer32 => int32, integer64 => int64

IMPLICIT NONE

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we have tested.
INTEGER, PARAMETER :: logical64 = integer64
INTEGER, PARAMETER :: logical32 = integer32

! Explicit kind type for reals used in routines imported by the UM and LFRic
INTEGER, PARAMETER :: real_jlslsm = real32

END MODULE um_types


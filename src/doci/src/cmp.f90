module qtzlibkindsmodule
    integer, parameter :: qt_k_ulong = selected_int_kind(9) ! 4-byte integer
    integer, parameter :: qt_k_ulongf = selected_int_kind(9) ! 4-byte integer
end module qtzlibkindsmodule

module qtzlibinterfacemodule
    ! fortran interface to zlib1.dll
    ! calling convention according
    ! ...libsrczlib1.2.3zlib-1.2.3win32dll_faq.txt
    !
    ! * the exported symbols are exclusively defined in the source
    ! files "zlib.h" and "zlib.def", found in an official zlib
    ! source distribution.
    ! * the symbols are exported by name, not by ordinal.
    ! * the exported names are undecorated.
    ! * the calling convention of functions is "c" (cdecl).
    ! * the zlib1.dll binary is linked to msvcrt.dll.
    !
    ! source:
    ! http://www.zlib.net/
    !
    use qtzlibkindsmodule
    implicit none
    
    ! ---------
    ! constants
    
    integer, parameter :: z_no_flush = 0 ! #define z_no_flush 0
    integer, parameter :: z_partial_flush = 1 ! #define z_partial_flush 1
    ! /* will be removed, use z_sync_flush instead */
    integer, parameter :: z_sync_flush = 2 ! #define z_sync_flush 2
    integer, parameter :: z_full_flush = 3 ! #define z_full_flush 3
    integer, parameter :: z_finish = 4 ! #define z_finish 4
    ! /* allowed flush values ; see deflate() below for details */
    
    integer, parameter :: z_ok = 0 ! #define z_ok 0
    integer, parameter :: z_stream_end = 1 ! #define z_stream_end 1
    integer, parameter :: z_need_dict = 2 ! #define z_need_dict 2
    integer, parameter :: z_errno =(-1) ! #define z_errno (-1)
    integer, parameter :: z_stream_error =(-2) ! #define z_stream_error (-2)
    integer, parameter :: z_data_error =(-3) ! #define z_data_error (-3)
    integer, parameter :: z_mem_error =(-4) ! #define z_mem_error (-4)
    integer, parameter :: z_buf_error =(-5) ! #define z_buf_error (-5)
    integer, parameter :: z_version_error =(-6) ! #define z_version_error (-6)
    ! /*return codes for the compression/decompression functions. negative
    ! * values are errors, positive values are used for special but normal events.
    ! */
    
    integer, parameter :: z_no_compression = 0 ! #define z_no_compression 0
    integer, parameter :: z_best_speed = 1 ! #define z_best_speed 1
    integer, parameter :: z_best_compression = 9 ! #define z_best_compression 9
    integer, parameter :: z_default_compression =(-1) ! #define z_default_compression (-1)
    ! /* compression levels */
    
    integer, parameter :: z_filtered = 1 ! #define z_filtered 1
    integer, parameter :: z_huffman_only = 2 ! #define z_huffman_only 2
    integer, parameter :: z_default_strategy = 0 ! #define z_default_strategy 0
    ! /* compression strategy ; see deflateinit2() below for details */
    
    integer, parameter :: z_binary = 0 ! #define z_binary 0
    integer, parameter :: z_ascii = 1 ! #define z_ascii 1
    integer, parameter :: z_unknown = 2 ! #define z_unknown 2
    ! /* possible values of the data_type field */
    
    integer, parameter :: z_deflated = 8 ! #define z_deflated 8
    ! /* the deflate compression method (the only one supported in this version) */
    
    integer, parameter :: z_null = 0 ! #define z_null 0 /* for initializing zalloc, zfree, opaque */
    
    ! #define zlib_version zlibversion()
    ! /* for compatibility with versions less than 1.0.2 */
    
    ! ----------
    ! interfaces
    
    ! ________
    ! int compress (bytef *dest, ulongf *destlen, const bytef *source, ulong sourcelen);
    ! compresses the source buffer into the destination buffer. sourcelen is the byte
    ! length of the source buffer.
    ! upon entry, destlen is the total size of the destination buffer, which must be at
    ! least 0.1% larger than sourcelen plus 12 bytes.
    ! upon exit, destlen is the actual size of the compressed buffer.
    interface
        ! int compress (bytef *dest, ulongf *destlen, const bytef *source, ulong sourcelen);
        integer function compress( dest, destlen, source, sourcelen )
            !x dec$ attributes c, dllimport :: compress
            !dec$ attributes c :: compress
            !dec$ attributes reference :: dest, destlen, source
            use qtzlibkindsmodule
            character(1), intent(out) :: dest(*)
            integer (qt_k_ulongf), intent(inout) :: destlen
            character(1), intent(in) :: source(*)
            integer (qt_k_ulong), intent(in) :: sourcelen
        end function compress
    end interface
    
    ! __________
    ! int uncompress (bytef *dest, ulongf *destlen, const bytef *source, ulong sourcelen);
    ! decompresses the source buffer into the destination buffer. sourcelen is the byte
    ! length of the source buffer.
    ! upon entry, destlen is the total size of the destination buffer, which must be
    ! large enough to hold the entire uncompressed data. (the size of the uncompressed
    ! data must have been saved previously by the compressor and transmitted to the
    ! decompressor by some mechanism outside the scope of this compression
    ! library.)
    ! upon exit, destlen is the actual size of the compressed buffer.
    interface
    ! int uncompress (bytef *dest, ulongf *destlen, const bytef *source, ulong sourcelen);
        integer function uncompress( dest, destlen, source, sourcelen )
            !x dec$ attributes c, dllimport :: uncompress
            !dec$ attributes c :: uncompress
            !dec$ attributes reference :: dest, destlen, source
            use qtzlibkindsmodule
            character(1), intent(out) :: dest(*)
            integer (qt_k_ulongf), intent(inout) :: destlen
            character(1), intent(in) :: source(*)
            integer (qt_k_ulong), intent(in) :: sourcelen
        end function uncompress
    end interface

contains

    character(100) function qtzl_cf_getzliberrormessage( ierror ) result(cmsg)
        integer, intent(in) :: ierror
        integer is
        
        cmsg = 'zlib :: '
        is = len_trim(cmsg) + 2
        select case ( ierror )
        case (z_ok)
        cmsg = ' '
        case (z_stream_end)
        cmsg(is:) = 'stream end (z_stream_end).'
        case (z_need_dict)
        cmsg(is:) = 'need dictionary (z_need_dict).'
        case (z_errno)
        cmsg(is:) = 'there is an error reading or writing the files (z_errno).'
        case (z_stream_error)
        cmsg(is:) = 'invalid compression level is supplied (z_stream_error).'
        case (z_data_error)
        cmsg(is:) = 'input data is corrupted or invalid (z_data_error).'
        case (z_mem_error)
        cmsg(is:) = 'memory could not be allocated for processing (z_mem_error).'
        case (z_buf_error)
        cmsg(is:) = 'not enough room in the output buffer (z_buf_error).'
        case (z_version_error)
        cmsg(is:) = 'version of zlib.h and the version of the library linked do not match (z_version_error).'
        case default
        write (cmsg(is:), 10000) ierror
10000 format('unknown error; error code = ', i0)
        end select
    
    end function qtzl_cf_getzliberrormessage

end module qtzlibinterfacemodule

!and here is a simple test program to show how to use it. sorry for the german texts, but since language doesn't matter here, i leave it as it is:

program t_zlib1
    use qtzlibinterfacemodule
    implicit none
    
    integer, parameter :: islen = 80
    character(1) :: c1asource(islen)
    character(islen) csource
    
    integer, parameter :: idlen = islen * (1. + 1./1000) + 12
    character(1) :: c1adest(idlen)
    equivalence(csource,c1asource)
    
    integer (qt_k_ulongf) :: idestlen
    integer (qt_k_ulong) :: isourcelen
    
    integer iret, j
    
    write(*,"(/'quelltext eingeben: ')")
    read(*,'(a)') csource
    
    isourcelen = len_trim(csource)
    idestlen = idlen
    iret = compress( c1adest, idestlen, c1asource, isourcelen )
    
    write(*,6000) isourcelen, idestlen
6000 format(/'laenge des urspruenglichen textes: ', i0/ &
    'komprimierter text (hexadez.), laenge:', i0)
    do j = 1, idestlen
        write(*,6001) c1adest(j)
6001 format(z2.2, 'h ')
    end do
    
    isourcelen = idestlen
    c1asource(1:isourcelen) = c1adest(1:idestlen)
    idestlen = idlen
    iret = uncompress( c1adest, idestlen, c1asource, isourcelen )
    
    c1asource(1:idestlen) = c1adest(1:idestlen) ! -> csource
    write(*,6010) isourcelen, idestlen, csource(1:idestlen)
6010 format(/'laenge des komprimierten textes: ', i0/ &
    'de-komprimierter text (laenge: ', i0, '): ', a)

end program t_zlib1
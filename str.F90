FUNCTION STR(K)
!   "Convert an integer to string."
    INTEGER, INTENT(IN) :: K
    CHARACTER(LEN=20) :: STR
    write (STR, *) K
    STR = adjustl(STR)
END FUNCTION STR
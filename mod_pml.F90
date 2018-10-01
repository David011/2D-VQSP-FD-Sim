!_______________________________________________________________________

MODULE PML
  USE PRECISION, ONLY: PP
  IMPLICIT NONE

  REAL (KIND=PP), DIMENSION (:,:), ALLOCATABLE ::                      &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
                                          ODA_ZH ,  ODB_ZH ,  ODC_ZH , &
                                          ODD_ZH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                                         ODD2_ZH 

  REAL (KIND=PP), DIMENSION (:,:), ALLOCATABLE ::                      &
                                        PZXL , PZXH ,                  &
                                        PXXL , PXXH ,                  &
                                        OXXL , OXXH ,                  &
                                               PXZH ,                  &
                                               PZZH ,                  &
                                               OZZH ,                  &
                                          RXXL , RXXH , SXXL , SXXH ,  &
                                          RZXL , RZXH ,                &
                                          GXXL , GXXH , HXXL , HXXH ,  &
                                          GZXL , GZXH,                 &
                                          RXZH ,                       &
                                          GXZH ,                       &
                                          RZZH , SZZH,                 &
                                          GZZH , HZZH


END MODULE PML

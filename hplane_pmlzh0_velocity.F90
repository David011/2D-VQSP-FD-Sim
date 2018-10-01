!=======================================================================

SUBROUTINE HPLANE_PMLZH0_VELOCITY                                   

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MZ, TPML
  USE WAVEFIELD   , ONLY:                                              &
                          UM, WM, QX, QZ
  USE AUXIL       , ONLY:                                              &
                           X_MIN  , X_MAX
  USE NONREF_BOUND, ONLY:                                              &
                          NR_DATA_Z,                                   &
                          SZ_U,       SZ_W , SZ_UP,        SZ_WP,      &
                          SZ_QX,      SZ_QZ, SZ_QXP,       SZ_QZP,     &
                          ABW01,   ABW02,   ABW10,   ABW11,            &
                          ABW12,   ABW20,   ABW21,   ABW22,            &
                          ABH01,   ABH02,   ABH10,   ABH11,            &
                          ABH12,   ABH20,   ABH21,   ABH22

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: I, LST, LST1, LST2

!-----------------------------------------------------------------------

  LST  = MZ + TPML
  LST1 = LST-1
  LST2 = LST-2

    DO  I = MAX(1-TPML,X_MIN), MIN(MX+TPML,X_MAX)

      UM(I,  LST) =                  ABH01 *    UM(I,     LST1)        &       
                                   + ABH02 *    UM(I,     LST2)        &
       + ABH10 * SZ_U (I  )%TIKMZ0 + ABH11 * SZ_U (I  )%TIKMZ1         &
                                   + ABH12 * SZ_U (I  )%TIKMZ2         &
       + ABH20 * SZ_UP(I  )%TIKMZ0 + ABH21 * SZ_UP(I  )%TIKMZ1         &
                                   + ABH22 * SZ_UP(I  )%TIKMZ2

      WM(I,  LST) =                  ABW01 *    WM(I,     LST1)        &
                                   + ABW02 *    WM(I,     LST2)        &
       + ABW10 * SZ_W (I  )%TIKMZ0 + ABW11 * SZ_W (I  )%TIKMZ1         &
                                   + ABW12 * SZ_W (I  )%TIKMZ2         &
       + ABW20 * SZ_WP(I  )%TIKMZ0 + ABW21 * SZ_WP(I  )%TIKMZ1         &
                                   + ABW22 * SZ_WP(I  )%TIKMZ2
      
      QX(I,  LST) =                   ABH01 *    QX(I,     LST1)        &   
                                    + ABH02 *    QX(I,     LST2)        &
       + ABH10 * SZ_QX (I  )%TIKMZ0 + ABH11 * SZ_QX (I  )%TIKMZ1        &
                                    + ABH12 * SZ_QX (I  )%TIKMZ2        &
       + ABH20 * SZ_QXP(I  )%TIKMZ0 + ABH21 * SZ_QXP(I  )%TIKMZ1        &
                                    + ABH22 * SZ_QXP(I  )%TIKMZ2
      
      QZ(I,  LST) =                   ABW01 *    QZ(I,     LST1)        &
                                    + ABW02 *    QZ(I,     LST2)        &
       + ABW10 * SZ_QZ (I  )%TIKMZ0 + ABW11 * SZ_QZ (I  )%TIKMZ1        &
                                    + ABW12 * SZ_QZ (I  )%TIKMZ2        &
       + ABW20 * SZ_QZP(I  )%TIKMZ0 + ABW21 * SZ_QZP(I  )%TIKMZ1        &
                                    + ABW22 * SZ_QZP(I  )%TIKMZ2
    END DO


END SUBROUTINE HPLANE_PMLZH0_VELOCITY

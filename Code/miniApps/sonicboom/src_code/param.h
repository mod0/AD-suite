      integer nsmax,ntmax,nsegmax,nfacmax, ndegmax
      parameter (nsmax=30000, ntmax=6*nsmax, nsegmax=7*nsmax)
      parameter (nfacmax=10000, ndegmax=100)

      real*8 stifsmx(nsmax)
      real*8 stifsmy(nsmax)
      real*8 stifsmz(nsmax)

      real*8 ktet1_1(ntmax),ktet1_2(ntmax),ktet1_3(ntmax)
      real*8 ktet1_4(ntmax),ktet1_5(ntmax)
      real*8 ktet1_6(ntmax),ktet1_7(ntmax),ktet1_8(ntmax)
      real*8 ktet1_9(ntmax),ktet1_10(ntmax)
      real*8 ktet1_11(ntmax),ktet1_12(ntmax)
      real*8 ktet2_2(ntmax),ktet2_3(ntmax),ktet2_4(ntmax)
      real*8 ktet2_5(ntmax),ktet2_6(ntmax)
      real*8 ktet2_7(ntmax),ktet2_8(ntmax),ktet2_9(ntmax)
      real*8 ktet2_10(ntmax),ktet2_11(ntmax)
      real*8 ktet2_12(ntmax)
      real*8 ktet3_3(ntmax),ktet3_4(ntmax),ktet3_5(ntmax)
      real*8 ktet3_6(ntmax),ktet3_7(ntmax)
      real*8 ktet3_8(ntmax),ktet3_9(ntmax),ktet3_10(ntmax)
      real*8 ktet3_11(ntmax),ktet3_12(ntmax)
      real*8 ktet4_4(ntmax),ktet4_5(ntmax),ktet4_6(ntmax)
      real*8 ktet4_7(ntmax),ktet4_8(ntmax)
      real*8 ktet4_9(ntmax),ktet4_10(ntmax), ktet4_11(ntmax)
      real*8 ktet4_12(ntmax), ktet5_5(ntmax)
      real*8 ktet5_6(ntmax),ktet5_7(ntmax),ktet5_8(ntmax)
      real*8 ktet5_9(ntmax),ktet5_10(ntmax)
      real*8 ktet5_11(ntmax),ktet5_12(ntmax)
      real*8 ktet6_6(ntmax),ktet6_7(ntmax),ktet6_8(ntmax)
      real*8 ktet6_9(ntmax),ktet6_10(ntmax)
      real*8 ktet6_11(ntmax),ktet6_12(ntmax)
      real*8 ktet7_7(ntmax),ktet7_8(ntmax),ktet7_9(ntmax)
      real*8 ktet7_10(ntmax),ktet7_11(ntmax)
      real*8 ktet7_12(ntmax)
      real*8 ktet8_8(ntmax),ktet8_9(ntmax),ktet8_10(ntmax)
      real*8 ktet8_11(ntmax),ktet8_12(ntmax)
      real*8 ktet9_9(ntmax),ktet9_10(ntmax),ktet9_11(ntmax)
      real*8 ktet9_12(ntmax)
      real*8 ktet10_10(ntmax),ktet10_11(ntmax),ktet10_12(ntmax)
      real*8 ktet11_11(ntmax),ktet11_12(ntmax)
      real*8 ktet12_12(ntmax)


      COMMON/mvmshto0/stifsmx, stifsmy, stifsmz
      COMMON/tmymvmsh0/ktet1_1,ktet1_2,ktet1_3,ktet1_4,ktet1_5,
     $                 ktet1_6,ktet1_7,
     $                 ktet1_8,ktet1_9,ktet1_10,ktet1_11,ktet1_12
      COMMON/tmymvmsh1/ktet2_2,ktet2_3,ktet2_4,ktet2_5,ktet2_6,
     $                 ktet2_7,ktet2_8,
     $                 ktet2_9,ktet2_10,ktet2_11,ktet2_12
      COMMON/tmymvmsh2/ktet3_3,ktet3_4,ktet3_5,ktet3_6,ktet3_7,
     $                 ktet3_8,ktet3_9,
     $                 ktet3_10,ktet3_11,ktet3_12
      COMMON/tmymvmsh3/ktet4_4,ktet4_5,ktet4_6,ktet4_7,ktet4_8,
     $                 ktet4_9,ktet4_10,
     $                 ktet4_11,ktet4_12
      COMMON/tmymvmsh4/ktet5_5,ktet5_6,ktet5_7,ktet5_8,ktet5_9,ktet5_10,
     $                 ktet5_11,ktet5_12
      COMMON/tmymvmsh5/ktet6_6,ktet6_7,ktet6_8,ktet6_9,ktet6_10,
     $                 ktet6_11,ktet6_12
      COMMON/tmymvmsh6/ktet7_7,ktet7_8,ktet7_9,ktet7_10,ktet7_11,
     $                 ktet7_12
      COMMON/tmymvmsh7/ktet8_8,ktet8_9,ktet8_10,ktet8_11,ktet8_12
      COMMON/tmymvmsh8/ktet9_9,ktet9_10,ktet9_11,ktet9_12
      COMMON/tmymvmsh9/ktet10_10,ktet10_11,ktet10_12
      COMMON/tmymvmsh10/ktet11_11,ktet11_12
      COMMON/tmymvmsh11/ktet12_12

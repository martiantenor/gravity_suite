      subroutine comelp(hk,ck,ce)
c     ==================================================
c     Purpose: Compute complete elliptic integrals K(k)
c              and E(k)
c     Input  : K  --- Modulus k ( 0 < k < 1 )
c     Output : CK --- K(k)
c              CE --- E(k)
c
c     From Shanjie Zhang and Jianming Jin's code, as
c     attached to their book "Computation of Special
c     Functions." That code is explicitly licensed for
c     use in other programs, provided acknowledgement
c     that the work is copyrighted:
c
c     Copyright (c) 1996 John Wiley & Sons, Inc.
c     =================================================

cf2py intent(in) HK
cf2py intent(out) CK,CE

      implicit double precision (a-h,o-z)

      pk = 1.0d0 - hk*hk

      if (hk.eq.1.0) then
         ck=1.0d+300
         ce=1.0d0
      else
         ak=(((.01451196212d0*pk+.03742563713d0)*pk
     &      +.03590092383d0)*pk+.09666344259d0)*pk+
     &      1.38629436112d0
         bk=(((.00441787012d0*pk+.03328355346d0)*pk+
     &      .06880248576d0)*pk+.12498593597d0)*pk+.5d0
         ck=ak-bk*dlog(pk)
         ae=(((.01736506451d0*pk+.04757383546d0)*pk+
     &      .0626060122d0)*pk+.44325141463d0)*pk+1.0d0
         be=(((.00526449639d0*pk+.04069697526d0)*pk+
     &      .09200180037d0)*pk+.2499836831d0)*pk
         ce=ae-be*dlog(pk)
      endif

      return
      end

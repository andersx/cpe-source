C =======================================================================
C Here is a random number generator from an expert. It has passed more 
C stringent tests than those found in Knuth's book. There are also versions 
C of the generator in C, Pascal, Modula-2, Ada, and Basic.
C
C - David LaSalle
C  SCRI
C  Florida State University
C  Tallahassee, FL  32306-4052
C  (904)644-8532
C
C arpanet: minuit@scri1.scri.fsu.edu  (128.186.4.1)
C         or
C         minuit%fsu@nmfecc.arpa
C =======================================================================
C This random number generator originally appeared in "Toward a Universal 
C Random Number Generator" by George Marsaglia and Arif Zaman. 
C Florida State University Report: FSU-SCRI-87-50 (1987)
C 
C It was later modified by F. James and published in "A Review of Pseudo-
C random Number Generators" 
C 
C THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
C       (However, a newly discovered technique can yield 
C         a period of 10^600. But that is still in the development stage.)
C
C It passes ALL of the tests for random number generators and has a period 
C   of 2^144, is completely portable (gives bit identical results on all 
C   machines with at least 24-bit mantissas in the floating point 
C   representation). 
C 
C The algorithm is a combination of a Fibonacci sequence (with lags of 97
C   and 33, and operation "subtraction plus one, modulo one") and an 
C   "arithmetic sequence" (using subtraction).
C
C On a Vax 11/780, this random number generator can produce a number in 
C    13 microseconds. 
C======================================================================== 
C
C Initialization of TEST - the clean way. 
C Modified by Dirk Porezag (porezag@physik.tu-chemnitz.de)
C
      block data 
      logical TEST
      common /raset4/ TEST
      data TEST /.FALSE./
      end
C
      subroutine RMARIN(IJ,KL)
C
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient 
C length to complete an entire calculation with. For example, if sveral 
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random 
C number generator can create 900 million different subsequences -- with 
C each subsequence having a length of approximately 10^30.
C 
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0
C
      implicit none
      real*8 U(97), C, CD, CM
      real*8 s, t
      integer I97, J97, IJ, KL
      integer i, j, k, l, ii, jj, m
      logical TEST
      common /raset1/ U
      common /raset2/ C, CD, CM
      common /raset3/ I97, J97
      common /raset4/ TEST

      if ( IJ .lt. 0  .or.  IJ .gt. 31328  .or.
     &     KL .lt. 0  .or.  KL .gt. 30081 ) then
        print *,'First seed must be in (0,31328)'
        print *,'Second seed must be in (0,30081)'
        stop
      endif

      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169) 

      do 2 ii = 1, 97
        s = 0.0
        t = 0.5
        do 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          if (mod(l*m, 64) .ge. 32) then
            s = s + t
          endif
          t = 0.5 * t
3       continue
        U(ii) = s
2     continue

      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0

      I97 = 97
      J97 = 33

      TEST = .TRUE.
      return
      end
C
C =======================================================================
C
      function ranmar()
C
C This is the random number generator proposed by George Marsaglia in 
C Florida State University Report: FSU-SCRI-87-50
C It was slightly modified by F. James to produce an array of pseudorandom
C numbers.
C
      real*8 U(97), C, CD, CM
      real*8 ranmar, uni
      integer I97, J97
      logical TEST
      common /raset1/ U
      common /raset2/ C, CD, CM
      common /raset3/ I97, J97
      common /raset4/ TEST
 
      if ( .NOT. TEST ) then
        print *,'Call the init routine RMARIN before calling RANMAR'  
        stop
      endif

      uni = U(I97) - U(J97)
      if ( uni .lt. 0.0 ) uni = uni + 1.0
      U(I97) = uni
      I97 = I97 - 1
      if (I97 .eq. 0) I97 = 97
      J97 = J97 - 1
      if (J97 .eq. 0) J97 = 97
      C = C - CD
      if ( C .lt. 0.0 ) C = C + CM
      uni = uni - C
      if ( uni .lt. 0.0 ) uni = uni + 1.0
      RANMAR = uni

      return
      end

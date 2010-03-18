#include <math.h>
#include <stdio.h>





inline int     Nint (float  d) { return (d>0) ? (int)(d+0.5) : -(int)(-d+0.5); }
int     odd     (int x){return(x&1);}

double MpMax(double a,double b)
{
	if (a > b)
		return a;
	else
		return b;
}

double MpMin(double a,double b)
{
	if (a > b)
		return b;
	else
		return a;
}









void ThreeJSymbolM (double l1, double l2, double l3, double m1, 
		    double &m2min, double &m2max, double *thrcof, int ndim, 
		    int &errflag)
{
  const double zero = 0.0, eps = 0.01, one = 1.0, two = 2.0;

  int nfin, nlim, i, n, index, lstep, nfinp1, nfinp2, nfinp3, nstep2;
  double oldfac, dv, newfac, sumbac = 0.0, thresh, a1s, sumfor, sumuni, 
         sum1, sum2, x, y, m2, m3, x1, x2, x3, y1, y2, y3, cnorm, 
         ratio, a1, c1, c2, c1old = 0.0, sign1, sign2;

  // Parameter adjustments
  --thrcof;

  errflag = 0;

  // "huge" is the square root of one twentieth of the largest floating
  // point number, approximately.
  double huge   = sqrt(10000000000000000000000000000.0 / 20.0),
         srhuge = sqrt(huge),
         tiny   = one / huge,
         srtiny = one / srhuge;

  // lmatch = zero

  //  Check error conditions 1, 2, and 3. 
  if (l1 - fabs(m1) + eps < zero 
      || fmod(l1 + fabs(m1) + eps, one) >= eps + eps) {
    errflag = 1;    
    printf("l1-abs(m1) less than zero or l1+abs(m1) not integer.");
    return;
  } else if (l1+l2-l3 < -eps || l1-l2+l3 < -eps || -(l1) + l2+l3 < -eps) {
    errflag = 2;
    printf( "l1, l2, l3 do not satisfy triangular condition.");
    return;
  } else if (fmod(l1 + l2 + l3 + eps, one) >= eps + eps) {
      errflag = 3;
      printf( "l1+l2+l3 not integer.");
      return;
  }

  // limits for m2 
  m2min = MpMax(-l2,-l3-m1);
  m2max = MpMin(l2,l3-m1);

  // Check error condition 4. 
  if (fmod(m2max - m2min + eps, one) >= eps + eps) {
    errflag = 4;
    printf( "m2max-m2min not integer.");
    return;
  }
  if (m2min < m2max - eps) goto L20;
  if (m2min < m2max + eps) goto L10;

  //  Check error condition 5. 
  errflag = 5;
 printf("m2min greater than m2max.");
  return;

  // This is reached in case that m2 and m3 can take only one value. 
L10:
  // mscale = 0 
  thrcof[1] = ((int(fabs(l2-l3-m1)+eps)%2 == 1) ? -one : one) / sqrt(l1+l2+l3+one);
  return;

  // This is reached in case that M1 and M2 take more than one value. 
L20:
  // mscale = 0 
  nfin = int(m2max - m2min + one + eps);
  if (ndim - nfin >= 0) goto L23;

  // Check error condition 6. 

  errflag = 6;
  printf("Dimension of result array for 3j coefficients too small.");
  return;

  //  Start of forward recursion from m2 = m2min 

L23:
  m2 = m2min;
  thrcof[1] = srtiny;
  newfac = 0.0;
  c1 = 0.0;
  sum1 = tiny;

  lstep = 1;
L30:
  ++lstep;
  m2 += one;
  m3 = -m1 - m2;

  oldfac = newfac;
  a1 = (l2 - m2 + one) * (l2 + m2) * (l3 + m3 + one) * (l3 - m3);
  newfac = sqrt(a1);

  dv = (l1+l2+l3+one) * (l2+l3-l1) - (l2-m2+one) * (l3+m3+one) 
                                   - (l2+m2-one) * (l3-m3-one);

  if (lstep - 2 > 0) c1old = fabs(c1);

// L32:
  c1 = -dv / newfac;

  if (lstep > 2) goto L60;

  //  If m2 = m2min + 1, the third term in the recursion equation vanishes,    
  //  hence 

  x = srtiny * c1;
  thrcof[2] = x;
  sum1 += tiny * c1 * c1;
  if (lstep == nfin) goto L220;
  goto L30;

L60:
  c2 = -oldfac / newfac;

  // Recursion to the next 3j coefficient 
  x = c1 * thrcof[lstep-1] + c2 * thrcof[lstep-2];
  thrcof[lstep] = x;
  sumfor = sum1;
  sum1 += x * x;
  if (lstep == nfin) goto L100;

  // See if last unnormalized 3j coefficient exceeds srhuge 

  if (fabs(x) < srhuge) goto L80;

  // This is reached if last 3j coefficient larger than srhuge, 
  // so that the recursion series thrcof(1), ... , thrcof(lstep) 
  // has to be rescaled to prevent overflow 

  // mscale = mscale + 1 
  for (i = 1; i <= lstep; ++i) {
    if (fabs(thrcof[i]) < srtiny) thrcof[i] = zero;
    thrcof[i] /= srhuge;
  }
  sum1 /= huge;
  sumfor /= huge;
  x /= srhuge;

  // As long as abs(c1) is decreasing, the recursion proceeds towards 
  // increasing 3j values and, hence, is numerically stable.  Once 
  // an increase of abs(c1) is detected, the recursion direction is 
  // reversed. 

L80:
  if (c1old - fabs(c1) > 0.0) goto L30;

  //  Keep three 3j coefficients around mmatch for comparison later 
  //  with backward recursion values. 

L100:
  // mmatch = m2 - 1 
  nstep2 = nfin - lstep + 3;
  x1 = x;
  x2 = thrcof[lstep-1];
  x3 = thrcof[lstep-2];

  //  Starting backward recursion from m2max taking nstep2 steps, so 
  //  that forwards and backwards recursion overlap at the three points 
  //  m2 = mmatch+1, mmatch, mmatch-1. 

  nfinp1 = nfin + 1;
  nfinp2 = nfin + 2;
  nfinp3 = nfin + 3;
  thrcof[nfin] = srtiny;
  sum2 = tiny;

  m2 = m2max + two;
  lstep = 1;
L110:
  ++lstep;
  m2 -= one;
  m3 = -m1 - m2;
  oldfac = newfac;
  a1s = (l2-m2+two) * (l2+m2-one) * (l3+m3+two) * (l3-m3-one);
  newfac = sqrt(a1s);
  dv = (l1+l2+l3+one) * (l2+l3-l1) - (l2-m2+one) * (l3+m3+one)
                                   - (l2+m2-one) * (l3-m3-one);
  c1 = -dv / newfac;
  if (lstep > 2) goto L120;

  // if m2 = m2max + 1 the third term in the recursion equation vanishes 

  y = srtiny * c1;
  thrcof[nfin - 1] = y;
  if (lstep == nstep2) goto L200;
  sumbac = sum2;
  sum2 += y * y;
  goto L110;

L120:
  c2 = -oldfac / newfac;

  // Recursion to the next 3j coefficient 

  y = c1 * thrcof[nfinp2 - lstep] + c2 * thrcof[nfinp3 - lstep];

  if (lstep == nstep2) goto L200;

  thrcof[nfinp1 - lstep] = y;
  sumbac = sum2;
  sum2 += y * y;

  // See if last 3j coefficient exceeds SRHUGE 

  if (fabs(y) < srhuge) goto L110;

  // This is reached if last 3j coefficient larger than srhuge, 
  // so that the recursion series thrcof(nfin), ... , thrcof(nfin-lstep+1)    
  // has to be rescaled to prevent overflow. 

  // mscale = mscale + 1 
  for (i = 1; i <= lstep; ++i) {
    index = nfin - i + 1;
    if (fabs(thrcof[index]) < srtiny) thrcof[index] = zero;
    thrcof[index] /= srhuge;
  }
  sum2 /= huge;
  sumbac /= huge;

  goto L110;

  //  The forward recursion 3j coefficients x1, x2, x3 are to be matched 
  //  with the corresponding backward recursion values y1, y2, y3. 

L200:
  y3 = y;
  y2 = thrcof[nfinp2-lstep];
  y1 = thrcof[nfinp3-lstep];

  //  Determine now ratio such that yi = ratio * xi  (i=1,2,3) holds 
  //  with minimal error. 

  ratio = (x1*y1 + x2*y2 + x3*y3) / (x1*x1 + x2*x2 + x3*x3);
  nlim = nfin - nstep2 + 1;

  if (fabs(ratio) < one) goto L211;
  for (n = 1; n <= nlim; ++n)
    thrcof[n] = ratio * thrcof[n];
  sumuni = ratio * ratio * sumfor + sumbac;
  goto L230;

L211:
  ++nlim;
  ratio = one / ratio;
  for (n = nlim; n <= nfin; ++n) 
    thrcof[n] = ratio * thrcof[n];
  sumuni = sumfor + ratio * ratio * sumbac;
  goto L230;

L220:
  sumuni = sum1;

  // Normalize 3j coefficients 

L230:
  cnorm = one / sqrt((l1+l1+one) * sumuni);

  // Sign convention for last 3j coefficient determines overall phase 
  if (thrcof[nfin] > 0)
	sign1 = 1;
  else
	sign1 = -1;
  sign2 = (int(fabs(l2-l3-m1)+eps)%2 == 1) ? -one : one;
  if (sign1 * sign2 <= 0.0) goto L235;
  else goto L236;

L235:
  cnorm = -cnorm;

L236:
  if (fabs(cnorm) < one) goto L250;

  for (n = 1; n <= nfin; ++n)
    thrcof[n] = cnorm * thrcof[n];
  return;

L250:
  thresh = tiny / fabs(cnorm);
  for (n = 1; n <= nfin; ++n) {
    if (fabs(thrcof[n]) < thresh) thrcof[n] = zero;
    thrcof[n] = cnorm * thrcof[n];
  }
} 








double ClebschGordan (double l1, double m1, double l2, double m2, 
		      double l3, double m3, int &errflag)
{
  const double err = 0.01;
  double CG = 0.0, m2min, m2max, *cofp;

  // static array for calculation of 3-j symbols
  const int ncof = 100;
  static double cof[ncof];

  // reset error flag
  errflag = 0;
  
  // Check for physical restriction. 
  // All other restrictions are checked by the 3-j symbol routine.
  if ( fabs(m1 + m2 - m3) > err) {
    errflag = 7;
    printf( "m1 + m2 - m3 is not zero.");
    return 0;
  } 
  
  // calculate minimum storage size needed for ThreeJSymbolM()
  // if the dimension becomes negative the 3-j routine will capture it
  int njm = Nint(MpMin(l2,l3-m1) - MpMax(-l2,-l3-m1) + 1); 
  
  // allocate dynamic memory if necessary
  cofp = (njm > ncof) ? new double[njm] : cof;

  // calculate series of 3-j symbols
  ThreeJSymbolM (l1,l2,l3,m1, m2min,m2max, cofp,njm, errflag);

  // calculated Clebsch-Gordan coefficient
  if (! errflag)
    CG = cofp[Nint(m2-m2min)] * (odd(Nint(l1-l2+m3)) ? -1 : 1) * sqrt(2*l3+1); 

  // free dynamic memory if necessary
  if (njm > ncof) delete [] cofp;

  return CG;
}



;;; Scheme bindings to  Multidimensional Interpolation with Radial Basis Functions 
;;; by John Burkardt 
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.

;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.

;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; -------------------------------------------------------------------------------

(module
 rbf
 
 (
  make-interpolant
  eval-interpolant
  )

  (import scheme (chicken base) (chicken foreign) (chicken format) srfi-4)

  
  (define-record-type <RBFInterpolant>
    (make-RBFInterpolant rank basis scale npoints xpoints weights)
    RBFInterpolant?
    (rank  RBFInterpolant-rank)
    (basis RBFInterpolant-basis)
    (scale  RBFInterpolant-scale)
    (npoints RBFInterpolant-npoints)
    (xpoints RBFInterpolant-xpoints)
    (weights RBFInterpolant-weights))


  
  #>
  #include <rbf_interp_nd.hpp>
  
  <#

  (define cphi1 (foreign-value "phi1" (function void (int (c-pointer double) double (c-pointer double)))))
  (define cphi2 (foreign-value "phi2" (function void (int (c-pointer double) double (c-pointer double)))))
  (define cphi3 (foreign-value "phi3" (function void (int (c-pointer double) double (c-pointer double)))))
  (define cphi4 (foreign-value "phi4" (function void (int (c-pointer double) double (c-pointer double)))))

  
;;   Evaluates a radial basis function interpolant.
;;    Input, int M, the spatial dimension.
;;  
;;      Input, int ND, the number of data points.
;;  
;;      Input, double XD[M*ND], the data points.
;;  
;;      Input, double R0, a scale factor.  R0 should be larger than the typical
;;      separation between points, but smaller than the maximum separation.
;;      The value of R0 has a significant effect on the resulting interpolant.
;;  
;;      Input, void PHI ( int N, double R[], double R0, double V[] ), a 
;;      function to evaluate the radial basis functions.
;;  
;;      Input, double W[ND], the weights, as computed by RBF_WEIGHTS.
;;  
;;      Input, int NI, the number of interpolation points.
;;  
;;      Input, double XI[M*NI], the interpolation points.
;;  
;;      Output, double RBF_INTERP_ND[NI], the interpolated values.
;;  void rbf_interp_nd ( int m, int nd, double xd[], double r0, 
;;                       void phi ( int n, double r[], double r0, double v[] ), double w[], 
;;                       int ni, double xi[], double fi[] );
  (define crbf_interp_nd (foreign-lambda void "rbf_interp_nd"
                                         int int f64vector double
                                         (function void (int (c-pointer double) double (c-pointer double)))
                                         f64vector int f64vector f64vector))

;;   Computes weights for radial basis function interpolation. 
;;      Input, int M, the spatial dimension.
;;  
;;      Input, int ND, the number of data points.
;;  
;;      Input, double XD[M*ND], the data points.
;;  
;;      Input, double R0, a scale factor.  R0 should be larger than the typical
;;      separation between points, but smaller than the maximum separation.
;;      The value of R0 has a significant effect on the resulting interpolant.
;;  
;;      Input, void PHI ( int N, double R[], double R0, double V[] ), a 
;;      function to evaluate the radial basis functions.
;;  
;;      Input, double FD[ND], the function values at the data points.
;;  
;;      Output, double RBF_WEIGHT[ND], the weights.
;;  void rbf_weight ( int m, int nd, double xd[], double r0, 
;;                    void phi ( int n, double r[], double r0, double v[] ), 
;;                    double fd[], double fi[] );

  (define crbf_weight (foreign-lambda void "rbf_weight"
                                      int int f64vector double
                                      (function void (int f64vector double f64vector))
                                      f64vector f64vector))

  (define (make-interpolant basis m scale xpoints ypoints)
    (let ((npoints (quotient (f64vector-length xpoints) m)))
      (assert (= npoints (f64vector-length ypoints)))
      (let ((weights (make-f64vector npoints))
            (basis-ptr (case basis
                         ((mq multiquadric) cphi1)
                         ((imq inverse-multiquadric) cphi2)
                         ((tp thin-plate) cphi3)
                         ((ga gaussian) cphi4)
                         (else (error 'make-interpolant "unknown basis" basis)))))
        (crbf_weight m npoints xpoints scale basis-ptr ypoints weights)
        (make-RBFInterpolant m basis scale npoints xpoints weights))
      ))

  
  (define (eval-interpolant interpolant ipoints)
    (let* (
           (m (RBFInterpolant-rank interpolant))
           (basis (RBFInterpolant-basis interpolant))
           (xpoints (RBFInterpolant-xpoints interpolant))
           (nd (RBFInterpolant-npoints interpolant))
           (scale (RBFInterpolant-scale interpolant))
           (weights (RBFInterpolant-weights interpolant))
           (ni (quotient (f64vector-length ipoints) m))
           (fpoints (make-f64vector ni))
           (basis-ptr (case basis
                        ((mq multiquadric) cphi1)
                        ((imq inverse-multiquadric) cphi2)
                        ((tp thin-plate) cphi3)
                        ((ga gaussian) cphi4)
                        (else (error 'make-interpolant "unknown basis" basis))))
           )
      (crbf_interp_nd m nd xpoints scale basis-ptr weights ni ipoints fpoints)
      fpoints))


  
  
  ) ; end of library


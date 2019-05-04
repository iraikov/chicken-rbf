(import scheme (chicken base) (chicken format) srfi-1 srfi-4 rbf)

(define (to-4-dp f)
  (/ (round (* f 10000)) 10000))

(define (=4 n1 n2)
  (= (to-4-dp n1) (to-4-dp n2)))


;(define (f x) (+ (sin (/ x 10.)) (expt (/ x 50.) 2.)))

(let* ((xpoints (list-tabulate 301 (lambda (x) x)))
       (ypoints (map (lambda (x) (* 0.3 x x)) xpoints))
       (ip (make-interpolant 'imq 1 1.
                             (list->f64vector xpoints)
                             (list->f64vector ypoints))))
  (let ((ipoints (list-tabulate 9991 (lambda (x) (* 0.1 x)))))
    (for-each (lambda (x y)
                (printf "eval ~A = ~A~%" x (f64vector-ref (eval-interpolant ip (f64vector x)) 0))
                (assert (=4 (f64vector-ref (eval-interpolant ip (f64vector x)) 0) y))
                )
              '(1 51 101)
              '(0.3 780.3 3060.3))
    ))


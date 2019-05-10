(import scheme (chicken base) (chicken format) srfi-1 srfi-4 rbf)

(define (to-4-dp f)
  (/ (round (* f 10000)) 10000))

(define (=4 n1 n2)
  (= (to-4-dp n1) (to-4-dp n2)))


(define (meshgrid x y)
  (let ((lenx (length x))
        (leny (length y)))
    (let recur ((i 0) (j 0) (xlst x) (ylst y) (ax '()))
      (if (< i lenx)
          (if (< j leny)
              (recur i (+ 1 j)
                     xlst (cdr ylst)
                     (append (reverse (list (car xlst) (car ylst))) ax))
              (recur (+ 1 i) 0 (cdr xlst) y ax))
          (reverse ax)))
    ))
          
  
(let* ((spacing 40.) (peak 20.)
       (xdim 100.) (ydim 100.) (xs spacing) (ys spacing) 
       (nxsteps (inexact->exact (+ 1 (round (/ (* 2 xdim) xs)))))
       (nysteps (inexact->exact (+ 1 (round (/ (* 2 ydim) ys)))))
       (xpoints (list-tabulate nxsteps (lambda (x) (- (* x xs) xdim))))
       (ypoints (list-tabulate nysteps (lambda (y) (- (* y ys) ydim))))
       (coords (meshgrid xpoints ypoints))
       (rpoints (list-tabulate (/ (length coords) 2) (lambda (i) peak)))
       (ip (make-interpolant 'ga 2 (/ spacing 2.71828)
                             (list->f64vector coords)
                             (list->f64vector rpoints))))
  (let* ((dx 1.0) (dy 1.0)
         (ixsteps (inexact->exact (+ 1 (round (/ (* 2 xdim) dx)))))
         (iysteps (inexact->exact (+ 1 (round (/ (* 2 ydim) dy)))))
         (ixpoints (list-tabulate ixsteps (lambda (x) (- (* x dx) xdim))))
         (iypoints (list-tabulate iysteps (lambda (y) (- (* y dy) ydim))))
         (icoords (meshgrid ixpoints iypoints))
         (y (eval-interpolant ip (list->f64vector icoords)))
         (n (/ (length icoords) 2)))

    (let recur ((i 0))
      (if (< i n)
          (begin
            (printf "~A ~A ~A~%"
                    (list-ref icoords (* i 2))
                    (list-ref icoords (+ 1 (* i 2)))
                    (f64vector-ref y i))
            (recur (+ 1 i)))))
  ))
  

#!/bin/bash
makePath="align_project"
binPath=$makePath"/build/bin"
imgPath="tests"

case $1 in

--clean)
	for i in {1..11}; do rm -rf $imgPath/out/"$i"; done
	for i in {1..6}; do rm -rf $imgPath/out/"$i"b; done
	rm -rf $imgPath/out/test/*
;;

--help)
	$binPath/./align --help
;;

--align) 
	for i in {1..11}; do mkdir -p $imgPath/out/"$i" && $binPath/./align $imgPath/in/"$i".bmp $imgPath/out/"$i"/align.bmp --align ; done
;;

--align-subpixel) 
	for i in {1..11}; do mkdir -p $imgPath/out/"$i" && $binPath/./align $imgPath/in/"$i".bmp $imgPath/out/"$i"/align_subpixel_$2.bmp --align --subpixel $2 ; done
;;

--align-subpixel-bicubic) 
	for i in {1..11}; do mkdir -p $imgPath/out/"$i" && $binPath/./align $imgPath/in/"$i".bmp $imgPath/out/"$i"/align_subpixel_"$2"_bicubic.bmp --align --subpixel $2 --bicubic-interp ; done
;;

--gray-world)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/gray_world.bmp --gray-world; done
;;

--custom)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/custom.bmp --custom $2; done
;;

--unsharp)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/unsharp.bmp --unsharp; done
;;

--median)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/median_$2.bmp --median $2; done
;;

--median-linear)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/median_linear_$2.bmp --median-linear $2; done
;;

--median-const)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/median_const_$2.bmp --median-const $2; done
;;

--gaussian)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/gaussian_$2_$3.bmp --gaussian $2 $3; done
;;

--gaussian-separable)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/gaussian_separable_$2_$3.bmp --gaussian-separable $2 $3; done
;;


--sobel-x)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/sobel-x.bmp --sobel-x; done
;;

--sobel-y)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/sobel-y.bmp --sobel-y; done
;;

--resize)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/resize_$2_linear.bmp --resize $2; done
;;

--resize-bicubic)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/resize_$2_bicubic.bmp --resize-bicubic $2; done
;;

--white-balance)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/white-balance.bmp --white-balance; done
;;

--autocontrast)
	for i in {1..11}; do $binPath/./align $imgPath/out/"$i"/align.bmp $imgPath/out/"$i"/autocontrast_$2.bmp --autocontrast $2; done
;; 

*)
	echo "Unknown command"
esac

exit 0
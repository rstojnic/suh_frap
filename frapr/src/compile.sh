#!/bin/bash
export PKG_LIBS=`Rscript -e "Rcpp:::LdFlags()"`
export PKG_CXXFLAGS="`Rscript -e "RcppArmadillo:::CxxFlags()"` `Rscript -e "Rcpp:::CxxFlags()"`"

rm -f gillespieSim.o
rm -f gillespieSim.so
rm -f gillespieSim3D.o
rm -f gillespieSim3D.so

R CMD SHLIB gillespieSim.cpp
R CMD SHLIB gillespieSim3D.cpp


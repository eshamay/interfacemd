#!/bin/bash
./curvefit $1 '0.5*(p1+p2)-0.5*(p1-p2)*tanh((x-x0)/d)' 'p1=600000;p2=0.0;x0=30.0;d=5.0'

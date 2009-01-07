#!/bin/bash
./curvefit $1 '0.5*(p1+p2)-0.5*(p1-p2)*tanh((x-x0)/d)' 'p1=1.0;p2=0.0;x0=0.0;d=5.0'

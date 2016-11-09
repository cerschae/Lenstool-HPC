# ksi = (a*a-b*b)/(a*a+b*b) (cos(2theta) + i sin(2theta))
awk '/#/ {print} !/#/ {q=$5/$4;e=(1-q*q)/(1+q*q);t=$6*3.145927/180;e1=e*cos(2*t);e2=e*sin(2*t);$4=e1;$5=e2;$6=e;print}' $1

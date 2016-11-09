# chi = (a*a-b*b)/(a*a+b*b) (cos(2theta) + i sin(2theta))
awk '/#/ {print} !/#/ {theta=0.5*atan2($5,$4)*180/3.1415927;e=sqrt($4*$4+$5*$5);a=sqrt(1+e);b=sqrt(1-e);$4=a;$5=b;$6=theta; print}' $1

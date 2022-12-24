%研究第一个Bohachevsky 1 function
function z=boh(x,y)

z=x.^2 + 2*y.^2 - 0.3*cos(3*pi*x) - 0.4*cos(4*pi*y) + 0.7;  %这里相当于输入的横纵坐标是向量


end
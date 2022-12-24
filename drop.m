%研究第三个drop wave function
function z=drop(x1,x2)

z=-(1+cos(12*sqrt(x1.^2+x2.^2)))./(0.5*(x1.^2+x2.^2)+2);  %这里相当于输入的横纵坐标是向量


end
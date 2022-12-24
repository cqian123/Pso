%研究第四个Easom function
function z=easom(x)
    z = 0.;
    for i=1:30    
        z = z+x(:,i).^2-10*cos(2*pi*x(:,i))+10;  %这里相当于输入的横纵坐标是向量
    end
end
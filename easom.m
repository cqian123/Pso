%�о����ĸ�Easom function
function z=easom(x)
    z = 0.;
    for i=1:30    
        z = z+x(:,i).^2-10*cos(2*pi*x(:,i))+10;  %�����൱������ĺ�������������
    end
end
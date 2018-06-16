function v = randcir( N )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

v=-2.*rand(N,2)+1;
radius=sqrt(v(:,1).^2+v(:,2).^2);
v=complex(v(:,1)./radius,v(:,2)./radius);

end


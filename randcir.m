function v = randcir( N )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

v=-2.*rand(N,2)+1;
radius=sqrt(v(:,1).^2+v(:,2).^2);
v=complex(v(:,1)./radius,v(:,2)./radius);

end


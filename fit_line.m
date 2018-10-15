function [a,resnorm] = fit_line(a0,data)
% ����ά���ݽ���ֱ����ϣ���ֱ����ϳ�������ʽ��  
% x-a(1)    y-a(2)  z-a(3)   
% ------ = ------ = ------   
%  a(4)     a(5)     a(6)   
% ԭ������� d �Ǹ����ݵ㵽���ֱ�ߵľ��룬���� lsqnonlin �� d ������С����   
% a0 �ǳ�ʼֵ��data ����ά���ݣ���һ���� x���ڶ����� y���������� z 
[a,resnorm] = lsqnonlin(@fit_line_fun,a0);   
   function d=fit_line_fun(a)   
	% �������������Ӻ���   
	xdata=data(1,:);     
	ydata=data(2,:);      
	zdata=data(3,:);
          
	point=a(1:3);      
	v=a(4:6);      
	d(1:length(xdata))=0;       
	for n=1:length(xdata)           
		m=[xdata(n);ydata(n);zdata(n)]-point(:);
		d(n)=norm(cross(m,v(:)))/norm(v(:));% ���ý������ε�֪ʶ�����d  
	end  
   end 
end


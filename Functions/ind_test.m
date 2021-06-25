function res=ind_test(V)
	T=length(V);
	J=zeros(T,4);
	for i = 2:T
		J(i,1)=V(i-1)==0 & V(i)==0;
		J(i,2)=V(i-1)==0 & V(i)==1;
		J(i,3)=V(i-1)==1 & V(i)==0;
		J(i,4)=V(i-1)==1 & V(i)==1;
	end
	V_00=sum(J(:,1));
	V_01=sum(J(:,2));
	V_10=sum(J(:,3));
	V_11=sum(J(:,4));
	p_00=V_00/(V_00+V_01);
	p_01=V_01/(V_00+V_01);
    if V_10 == 0
        p_10 = 0;
    else
        p_10=V_10/(V_10+V_11);
    end
    if V_11==0
        p_11 = 0;
    else
	p_11=V_11/(V_10+V_11);
    end
	hat_p=(V_01+V_11)/(V_00+V_01+V_10+V_11);
    al = (1-hat_p)^(V_00+V_10)*(hat_p)^(V_01+V_11);
    bl = (p_00)^(V_00) *(p_01)^V_01 *(p_10)^V_10* (p_11)^V_11;
    res= -2*log(al/bl);
 
out=res;


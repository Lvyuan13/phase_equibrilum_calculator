#include<stdio.h>
#include<iostream>
#include<math.h>
using namespace std;
//double tc1=154.58,pc1=5043,tc2=126.15,pc2=3394,w1=0.019,w2=0.045,R=8.3144;
double t=295;//在t为110K的情况下，变动x。求取P,y;
double tc2=369.825,pc2=4247.66,tc1=304.1282 ,pc1=7377.3,w2=0.152,w1=0.22394;//数据取自refprp NIST 9.0
#define R 8.3144
double kij=0;//refpro NIST 9.0数据
int main() {
	FILE*fp;
	fp=fopen("satured mixture2.txt","w");
	double x1,x2,y1,y2;
//double ya1[100],ya2[100];
	double p;
	int i;
	double Y;
	double Y1,Y2;
	//double ddt;
	double paa,pbb;
	double faa,fbb;
	double L;
	double fcc;
	double sum_y(double,double,double);
	void calculate_y(double,double,double,double*);
	double y[2];
	double ddp=100;
	int it;//防止死循环
	for(i=0; i<=100;i++) {
		x1=0.01*i;
		x2=1-x1;
		for(p=100; p<=20000; p=p+ddp) {
			//	if(p<100)  ddt=1;
			//	if(p>1000) ddt=0.1;
			//	if(p>4000) ddt=0.01;
			//	if(p>4600) ddt=0.001;
			//	if(p>4680) ddt=0.001;

			faa=sum_y(t,x1,p);
			fbb=sum_y(t,x1,p+ddp);
			if((faa*fbb)<=0) {
				paa=p;
				pbb=p+ddp;
				it=0;
				do {
					faa=sum_y(t,x1,paa);
					fbb=sum_y(t,x1,pbb);
					L=(faa*pbb-fbb*paa)/(faa-fbb);
					fcc=sum_y(t,x1,L);
					if(fcc*faa<0)
					pbb=L;
					else
					paa=L;
					it++;
					if (it>=80) break;
				} while(fabs(fcc)>1e-10) ;

				break;
			}
			else {
                //do nothing
			}

		}
       //  cout<<"寻根区间为（"<<paa<<","<<pbb<<")"<<endl;
         calculate_y(t,x1,L,y);
         y1=y[0];
         y2=y[1];
         double F(double t,double x1,double x2,double y1,double y2,double p,int j);
         double fl1=F(t,x1,x2,y1,y2,L,1);
         double fv1=F(t,x1,x2,y1,y2,L,2);
         double fl2=F(t,x1,x2,y1,y2,L,3);
         double fv2=F(t,x1,x2,y1,y2,L,4);
         cout<<"逸度::phil1,phiv1,phil2,phiv2  "<<x1*fl1<<"  "<<y1*fv1<<"  "<<x2*fl2<<"  "<<y2*fv2<<"  "<<endl;
         // 如果逸度相等，就说明算法收敛
	    printf("x1=%lf,x2=%lf,y1=%lf,y2=%lf,p=%lf,t=%lf\n",x1,x2,y1,y2,L,t);
		fprintf(fp,"x1=%lf,x2=%lf,y1=%lf,y2=%lf,p=%lf,t=%lf\n",x1,x2,y1,y2,L,t);
		//i++;

	}

	fclose(fp);
	return 0;
}

//double sum_y(double t,double x1,double p)
double sum_y(double t,double x1,double p){
	void calculate_y(double,double,double,double*);

	double sum;
	//double y1,y2;
	double y[2];
	calculate_y(t,x1,p,y);
	if (y[0]<0) y[0]=fabs(y[0]);
	if (y[1]<0) y[1]=fabs(y[1]);
	sum=y[0]+y[1]-1;
	return sum;

}


void calculate_y(double t,double x1,double p,double*y){
int i,j;
	double ya1[100],ya2[100];
	double fl1,fl2,fv1,fv2;
	double F(double,double,double,double,double,double,int);
	double x2=1-x1;
	double y1,y2;
	y1=x1;
	y2=x2;
	ya1[0]=y1;
	ya2[0]=y2;
	for(i=1; i<100; i++) { //循环一百次，应该可以迭代出收敛的根，如果不能，则需要增大循环次数；

		j=1;    //j是控制返回值的变量；
		fl1=F(t,x1,x2,y1,y2,p,j);
//		cout<<"fl1为"<<fl1<<endl;
		j=2;
		fv1=F(t,x1,x2,y1,y2,p,j);
//		cout<<"fv1为"<<fv1<<endl;
		j=3;
		fl2=F(t,x1,x2,y1,y2,p,j);
//		cout<<"fl2为"<<fl2<<endl;
		j=4;
		fv2=F(t,x1,x2,y1,y2,p,j);
		ya1[i]=x1*fl1/fv1;
		ya2[i]=x2*fl2/fv2;
		y1=ya1[i]/(ya1[i]+ya2[i]);
		y2=ya2[i]/(ya1[i]+ya2[i]);
		if((fabs(ya1[i]-ya1[i-1])<=1e-10)&&(fabs(ya2[i]-ya2[i-1])<=1e-10))
			break;
	}
y[0]=ya1[i];
y[1]=ya2[i];

}
double F(double t,double x1,double x2,double y1,double y2,double p,int j) {//j为控制返回值的变量
	double  Al,Bl,Av,Bv;//PR方程中的量；
	double ac1,ac2;
	double alpha1,alpha2;
	double b1,b2;
	double m1,m2;
	double a1,a2;
	double bl,bv;
	double al,av;
	double zl,zv;
	double fl1,fl2,fv1,fv2;
	ac1=0.457235*pow((R*tc1),2)/pc1;
	ac2=0.457235*pow((R*tc2),2)/pc2;
	m1=0.37464+1.54226*w1-0.26992*pow(w1,2);
	m2=0.37464+1.54662*w1-0.26992*pow(w1,2);
	alpha1=pow((1+m1*(1-sqrt(t/tc1))),2.0);
	alpha2=pow((1+m2*(1-sqrt(t/tc2))),2.0);
	a1=ac1*alpha1;
	a2=ac2*alpha2;
	double root_l(double,double,double,double);
	double root_v(double,double,double,double);
	b1=0.077796*R*tc1/pc1;
	b2=0.077796*R*tc2/pc2;
	al=x1*x1*a1+2*x1*x2*sqrt(a1*a2)*(1-kij)+x2*x2*a2; //al
	av=y1*y1*a1+2*y1*y2*sqrt(a1*a2)*(1-kij)+y2*y2*a2;//av
	bl=x1*b1+x2*b2;
	bv=y1*b1+y2*b2;
	Al=al*p/pow((R*t),2.0);
	Av=av*p/pow((R*t),2.0);
	Bl=bl*p/(R*t);
	Bv=bv*p/(R*t);
	zl=root_l(1,Bl-1,Al-3*Bl*Bl-2*Bl,-Al*Bl+Bl*Bl+Bl*Bl*Bl);
	zv=root_v(1,Bv-1,Av-3*Bv*Bv-2*Bv,-Av*Bv+Bv*Bv+Bv*Bv*Bv);
	fl1=exp(b1/bl*(zl-1)-log(zl-Bl)+Al/(2*sqrt(2)*Bl)*(b1/bl-2/al*(x1*a1+x2*(1-kij)*sqrt(a1*a2)))*log((zl+(sqrt(2)+1)*Bl)/(zl-(sqrt(2)-1)*Bl)));
	fv1=exp(b1/bv*(zv-1)-log(zv-Bv)+Av/(2*sqrt(2)*Bv)*(b1/bv-2/av*(y1*a1+y2*(1-kij)*sqrt(a1*a2)))*log((zv+(sqrt(2)+1)*Bv)/(zv-(sqrt(2)-1)*Bv)));
	fl2=exp(b2/bl*(zl-1)-log(zl-Bl)+Al/(2*sqrt(2)*Bl)*(b2/bl-2/al*(x2*a2+x1*(1-kij)*sqrt(a1*a2)))*log((zl+(sqrt(2)+1)*Bl)/(zl-(sqrt(2)-1)*Bl)));
	fv2=exp(b2/bv*(zv-1)-log(zv-Bv)+Av/(2*sqrt(2)*Bv)*(b2/bv-2/av*(y2*a2+y1*(1-kij)*sqrt(a1*a2)))*log((zv+(sqrt(2)+1)*Bv)/(zv-(sqrt(2)-1)*Bv)));
//aij在这里为x1*a1+x2*sqrt(a1*a2)
	if(j==1)
		return fl1;
	if(j==2)
		return fv1;
	if(j==3)
		return fl2;
	if(j==4)
		return fv2;
}

double root_l(double a,double b,double c,double d) {
//double a,b,c,d; //作为主函数时需要显示，作为子函数时不用显示；
	double L,A,B,C; //L为判别式，A,B,C为判别符号
	double X1,X2,X3; //分别为三个解
	double Y1,Y2,Y3; //为在虚根情况下使用的符号，在运算中不体现，出于完整性考虑
	double P1,P2,P3; //在存在虚根的情况下使用的符号，分别表示虚根的实部和虚部，在运算中不体现，出于完整性考虑
	double GS1,GS2;  //表示三次根号下Y1，Y2,实际运算中不体现，出于完整性考虑
	double K;
	double T;
	double rad;
	double Xmax,Xmin;    //用到的变量声明
	double max(double a,double b,double c);
	double min(double a,double b,double c);   //两个子函数的声明
//printf("please input a b c d\n");         //作为子函数时隐藏
//scanf("%lf%lf%lf%lf",&a,&b,&c,&d);         //作为子函数时隐藏
	/*	printf("a=%lf b=%lf c=%lf d=%lf\n",a,b,c,d);  */
	A=b*b-3*a*c;
	B=b*c-9*a*d;
	C=c*c-3*b*d;
	L=B*B-4*A*C;// 判别式
//	printf("A=%lf,B=%lf,C=%lf\n",A,B,C);
//	printf("L=%lf\n",L);
	if(A==0&&B==0) {
		X1=X2=X3=-b/(3*a);
//		printf("X1=X2=X3=%.3f\n",X1);
		Xmax=X1;
		Xmin=X1;


	}
//-------------------第一种情况，三重根 ------------------------------------------------------------//
	else    {
		if(L>0) {
			Y1=A*b+(3*a*(-B+sqrt(L)))/2;
			Y2=A*b+(3*a*(-B-sqrt(L)))/2;
//			printf("Y1=%.2f Y2=%.2f\n\n",Y1,Y2);
			if(Y1>0)
				GS1=pow(Y1,1.0/3);        //GS1为三次根号Y1
			else
				GS1=-pow(fabs(Y1),1.0/3);
			if(Y2>0)
				GS2=pow(Y2,(1.0/3));      //GS2为三次根号Y2
			else     GS2=-pow(fabs(Y2),(1.0/3));
			P1=(-2*b+GS1+GS2)/(6*a);
			P2=sqrt(3)*(GS1-GS2)/(6*a);
			X1=(-b-GS1-GS2)/(3*a);
//			printf("X1=%.3f\n",X1);
//			printf("X2=%.3f+%.3fi\n",P1,P2);
//			printf("X3=%.3f-%.3fi\n  X1and X2 are flase roots which can be overlooked",P1,P2);
			return 0;      //此时不存在Xmax,Xmin,返回0即可；
		}
		//存在两个虚根，不返回值，无实际意义，为求根函数的完整性而做
//----------------第二种情况，一实根，两虚根---------------------------------------------------------//
		else if(L==0) {
			K=B/A;
			X1=-b/a+K;
			X2=X3=-K/2;
//			printf("X1=%.3f\n",X1);
//			printf("X2=X3=%.3f\n",X2);
			Xmax=max(X1,X2,X3);
//			printf("Xmax=%.3f",Xmax);

		}
//-----------------第三种情况，两个实根，其中一根为重根---------------------------------------------------------//
		else if(L<0) {
			T=(2*A*b-3*a*B)/(2*pow(A,1.5));
			rad=acos(T);
			X1=(-b-2*sqrt(A)*cos(rad/3))/(3*a);
			X2=(-b+sqrt(A)*(cos(rad/3)+sqrt(3)*sin(rad/3)))/(3*a);
			X3=(-b+sqrt(A)*(cos(rad/3)-sqrt(3)*sin(rad/3)))/(3*a);
//			printf("X1=%.3f\n",X1);
//			printf("X2=%.3f\n",X2);
//			printf("X3=%.3f",X3);
			Xmax=max(X1,X2,X3);
//			printf("\nXmax=%.3f",Xmax);
			Xmin=min(X1,X2,X3);
//			printf("\nXmin=%.3f",Xmin);
		}
	}
	return (Xmin);     //作为子函数时返回Xmin
}
//------------------------------------------------root_v函数-----------------------------------------//
double root_v(double a,double b,double c,double d) {
//double a,b,c,d; //作为主函数时需要显示，作为子函数是不用显示；
	double L,A,B,C; //L为判别式，A,B,C为判别符号
	double X1,X2,X3; //分别为三个解
	double Y1,Y2,Y3; //为在虚根情况下使用的符号，在运算中不体现，出于完整性考虑
	double P1,P2,P3; //在存在虚根的情况下使用的符号，分别表示虚根的实部和虚部，在运算中不体现，出于完整性考虑
	double GS1,GS2;  //表示三次根号下Y1，Y2,实际运算中不体现，出于完整性考虑
	double K;
	double T;
	double rad;
	double Xmax,Xmin;    //用到的变量声明
	double max(double a,double b,double c);
	double min(double a,double b,double c);   //两个子函数的声明
//printf("please input a b c d\n");         //作为子函数时隐藏
//scanf("%lf%lf%lf%lf",&a,&b,&c,&d);         //作为子函数时隐藏
//printf("a=%lf b=%lf c=%lf d=%lf\n",a,b,c,d);
	A=b*b-3*a*c;
	B=b*c-9*a*d;
	C=c*c-3*b*d;
	L=B*B-4*A*C;// 判别式
//	printf("A=%lf,B=%lf,C=%lf\n",A,B,C);
//	printf("L=%lf\n",L);
	if(A==0&&B==0) {
		X1=X2=X3=-b/(3*a);
//		printf("X1=X2=X3=%.3f\n",X1);
		Xmax=X1;
		Xmin=X1;


	}
//-------------------第一种情况，三重根 ------------------------------------------------------------//
	else    {
		if(L>0) {
			Y1=A*b+(3*a*(-B+sqrt(L)))/2;
			Y2=A*b+(3*a*(-B-sqrt(L)))/2;
//			printf("Y1=%.2f Y2=%.2f\n\n",Y1,Y2);
			if(Y1>0)
				GS1=pow(Y1,1.0/3);        //GS1为三次根号Y1
			else
				GS1=-pow(fabs(Y1),1.0/3);
			if(Y2>0)
				GS2=pow(Y2,(1.0/3));      //GS2为三次根号Y2
			else     GS2=-pow(fabs(Y2),(1.0/3));
			P1=(-2*b+GS1+GS2)/(6*a);
			P2=sqrt(3)*(GS1-GS2)/(6*a);
			X1=(-b-GS1-GS2)/(3*a);
//			printf("X1=%.3f\n",X1);
//			printf("X2=%.3f+%.3fi\n",P1,P2);
//			printf("X3=%.3f-%.3fi\n  X1and X2 are flase roots which can be overlooked",P1,P2);
			return 0;
		}
//----------------第二种情况，一实根，两虚根---------------------------------------------------------//
		else if(L==0) {
			K=B/A;
			X1=-b/a+K;
			X2=X3=-K/2;
			/*printf("X1=%.3f\n",X1);
			printf("X2=X3=%.3f\n",X2);*/
			Xmax=max(X1,X2,X3);
			//printf("Xmax=%.3f",Xmax);//     作为子函数时不用此条命令
			Xmin=min(X1,X2,X3);
			//printf("Xmin=%.3f",Xmin);  //   作为子函数时不用此条命令
		}
//-----------------第三种情况，两个实根，其中一根为重根---------------------------------------------------------//
		else if(L<0) {
			T=(2*A*b-3*a*B)/(2*pow(A,1.5));
			rad=acos(T);
			X1=(-b-2*sqrt(A)*cos(rad/3))/(3*a);
			X2=(-b+sqrt(A)*(cos(rad/3)+sqrt(3)*sin(rad/3)))/(3*a);
			X3=(-b+sqrt(A)*(cos(rad/3)-sqrt(3)*sin(rad/3)))/(3*a);
//			printf("X1=%.3f\n",X1);
//			printf("X2=%.3f\n",X2);
//			printf("X3=%.3f",X3);
			Xmax=max(X1,X2,X3);
			//printf("\nXmax=%.3f",Xmax);//作为子函数时隐藏此命令
			Xmin=min(X1,X2,X3);
//			printf("\nXmin=%.3f",Xmin);
		}
	}
	return (Xmax);     //作为子函数时返回Xmin
}
//--------------------------下为选出Xmax和Xmin的子函数-----------------------------------------//
double max(double a,double b,double c) {
	double t;
	t=a;
	if (t<b)
		t=b;
	if(t<c)
		t=c;
	return(t);
}
double min(double a,double b,double c) {
	double min;
	min=a;
	if (min>b)
		min=b;
	if(min>c)
		min=c;
	return(min);
}

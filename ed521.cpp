// Test program for ed521 scalar point multiplication
// Uses Bernstein et al.s point multiplication method from Curve41417 paper
// Fully Tested and debugged
// g++ -O3 ed521.cpp -o ed521
// M.Scott 27/8/2014

#include <iostream>
#include <ctime>
#include <inttypes.h>

using namespace std;

typedef __int128 type128;
typedef int64_t type64;

static const type64 bot58bits = 0x3ffffffffffffff;
static const type64 bot52bits = 0xfffffffffffff;

#include <stdint.h>

   __inline__ uint64_t rdtsc() {
   uint32_t lo, hi;
   __asm__ __volatile__ (      // serialize
     "xorl %%eax,%%eax \n        cpuid"
     ::: "%rax", "%rbx", "%rcx", "%rdx");
   
   __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
   return (uint64_t)hi << 32 | lo;
   }

// w=x+y
void gadd(type64 *x,type64 *y,type64 *w)
{
	w[0]=x[0]+y[0];
	w[1]=x[1]+y[1];
	w[2]=x[2]+y[2];
	w[3]=x[3]+y[3];
	w[4]=x[4]+y[4];
	w[5]=x[5]+y[5];
	w[6]=x[6]+y[6];
	w[7]=x[7]+y[7];
	w[8]=x[8]+y[8];
}

// w=x-y
void gsub(type64 *x,type64 *y,type64 *w)
{
	w[0]=x[0]-y[0];
	w[1]=x[1]-y[1];
	w[2]=x[2]-y[2];
	w[3]=x[3]-y[3];
	w[4]=x[4]-y[4];
	w[5]=x[5]-y[5];
	w[6]=x[6]-y[6];
	w[7]=x[7]-y[7];
	w[8]=x[8]-y[8];
}

// w-=x
void gdec(type64 *x,type64 *w)
{
	w[0]-=x[0];
	w[1]-=x[1];
	w[2]-=x[2];
	w[3]-=x[3];
	w[4]-=x[4];
	w[5]-=x[5];
	w[6]-=x[6];
	w[7]-=x[7];
	w[8]-=x[8];
}

// w=x
void gcopy(type64 *x,type64 *w)
{
	w[0]=x[0];
	w[1]=x[1];
	w[2]=x[2];
	w[3]=x[3];
	w[4]=x[4];
	w[5]=x[5];
	w[6]=x[6];
	w[7]=x[7];
	w[8]=x[8];
}

// w*=2
void gmul2(type64 *w)
{
	w[0]*=2;
	w[1]*=2;
	w[2]*=2;
	w[3]*=2;
	w[4]*=2;
	w[5]*=2;
	w[6]*=2;
	w[7]*=2;
	w[8]*=2;
}

// w-=2*x
void gsb2(type64 *x,type64 *w)
{
	w[0]-=2*x[0];
	w[1]-=2*x[1];
	w[2]-=2*x[2];
	w[3]-=2*x[3];
	w[4]-=2*x[4];
	w[5]-=2*x[5];
	w[6]-=2*x[6];
	w[7]-=2*x[7];
	w[8]-=2*x[8];
}

// reduce w - Short Coefficient Reduction
void scr(type64 *w)
{
	type64 t0,t1,t2;
	t0=w[0]&bot58bits;

	t1=w[1]+(w[0]>>58);
	w[1]=t1&bot58bits;

	t2=w[2]+(t1>>58);
	w[2]=t2&bot58bits;

	t1=w[3]+(t2>>58);
	w[3]=t1&bot58bits;

	t2=w[4]+(t1>>58);
	w[4]=t2&bot58bits;

	t1=w[5]+(t2>>58);
	w[5]=t1&bot58bits;

	t2=w[6]+(t1>>58);
	w[6]=t2&bot58bits;

	t1=w[7]+(t2>>58);
	w[7]=t1&bot58bits;

	t2=w[8]+(t1>>58);
	w[8]=t2&bot58bits;

	w[0]=t0+2*(t2>>58);	
}

// multiply w by a constant, w*=i

void gmuli(type64 *w,int i)
{
	type128 t;

	t=(type128)w[0]*i;
	w[0]=((type64)t)&bot58bits;

	t=(type128)w[1]*i+(t>>58);
	w[1]=((type64)t)&bot58bits;

	t=(type128)w[2]*i+(t>>58);
	w[2]=((type64)t)&bot58bits;

	t=(type128)w[3]*i+(t>>58);
	w[3]=((type64)t)&bot58bits;

	t=(type128)w[4]*i+(t>>58);
	w[4]=((type64)t)&bot58bits;

	t=(type128)w[5]*i+(t>>58);
	w[5]=((type64)t)&bot58bits;

	t=(type128)w[6]*i+(t>>58);
	w[6]=((type64)t)&bot58bits;

	t=(type128)w[7]*i+(t>>58);
	w[7]=((type64)t)&bot58bits;

	t=(type128)w[8]*i+(t>>58);
	w[8]=((type64)t)&bot58bits;

	w[0]+=2*(t>>58);
}

// z=x^2

void gsqr(type64 *x,type64 *z)
{
	type128 t0,t1,t2;
	type128 st0;
	t1=2*((type128)x[0]*x[8]+(type128)x[1]*x[7]+(type128)x[2]*x[6]+(type128)x[3]*x[5])   +(type128)x[4]*x[4];
	st0=((type64) t1)&bot58bits;
	t2=4*((type128)x[1]*x[8]+(type128)x[2]*x[7]+(type128)x[3]*x[6]+(type128)x[4]*x[5])   +(type128)x[0]*x[0]+2*(t1>>58);
	z[0]=((type64) t2)&bot58bits;
	t1=4*((type128)x[2]*x[8]+(type128)x[3]*x[7]+(type128)x[4]*x[6])     +2*((type128)x[0]*x[1]+(type128)x[5]*x[5])  +(t2>>58);
	z[1]=((type64) t1)&bot58bits;
	t2=4*((type128)x[3]*x[8]+(type128)x[4]*x[7]+(type128)x[5]*x[6])     +2*(type128)x[0]*x[2]+(type128)x[1]*x[1]+(t1>>58);
	z[2]=((type64) t2)&bot58bits;
	t1=4*((type128)x[4]*x[8]+(type128)x[5]*x[7])      +2*((type128)x[0]*x[3]+(type128)x[1]*x[2]+(type128)x[6]*x[6])  +(t2>>58);
	z[3]=((type64) t1)&bot58bits;
	t2=4*((type128)x[5]*x[8]+(type128)x[6]*x[7])+2*((type128)x[0]*x[4]+(type128)x[1]*x[3])+(type128)x[2]*x[2]+(t1>>58);
	z[4]=((type64) t2)&bot58bits;
	t1=4*(type128)x[6]*x[8]+2*((type128)x[0]*x[5]+(type128)x[1]*x[4]+(type128)x[2]*x[3]+(type128)x[7]*x[7])+(t2>>58);
	z[5]=((type64) t1)&bot58bits;
	t2=4*(type128)x[7]*x[8]+2*((type128)x[0]*x[6]+(type128)x[1]*x[5]+(type128)x[2]*x[4])+(type128)x[3]*x[3]+(t1>>58);
	z[6]=((type64) t2)&bot58bits;
	t1=2*((type128)x[0]*x[7]+(type128)x[1]*x[6]+(type128)x[2]*x[5]+(type128)x[3]*x[4]+(type128)x[8]*x[8])+(t2>>58);
	z[7]=((type64) t1)&bot58bits;
	st0+=(t1>>58);
	z[8]=((type64)st0)&bot58bits;
	z[0]+=2*(type64)(st0>>58);
}

// z=2x^2

void gsqr2(type64 *x,type64 *z)
{
	type128 t0,t1,t2;
	type128 st0;
	t1=4*((type128)x[0]*x[8]+(type128)x[1]*x[7]+(type128)x[2]*x[6]+(type128)x[3]*x[5])   +(type128)x[4]*(2*x[4]);
	st0=((type64) t1)&bot58bits;
	t2=8*((type128)x[1]*x[8]+(type128)x[2]*x[7]+(type128)x[3]*x[6]+(type128)x[4]*x[5])   +2*((type128)x[0]*x[0]+(t1>>58));
	z[0]=((type64) t2)&bot58bits;
	t1=8*((type128)x[2]*x[8]+(type128)x[3]*x[7]+(type128)x[4]*x[6])     +4*((type128)x[0]*x[1]+(type128)x[5]*x[5])  +(t2>>58);
	z[1]=((type64) t1)&bot58bits;
	t2=8*((type128)x[3]*x[8]+(type128)x[4]*x[7]+(type128)x[5]*x[6])     +4*(type128)x[0]*x[2]+(type128)x[1]*(2*x[1])+(t1>>58);
	z[2]=((type64) t2)&bot58bits;
	t1=8*((type128)x[4]*x[8]+(type128)x[5]*x[7])      +4*((type128)x[0]*x[3]+(type128)x[1]*x[2]+(type128)x[6]*x[6])  +(t2>>58);
	z[3]=((type64) t1)&bot58bits;
	t2=8*((type128)x[5]*x[8]+(type128)x[6]*x[7])  +4*((type128)x[0]*x[4]+(type128)x[1]*x[3])  +(type128)x[2]*(2*x[2])  +(t1>>58);
	z[4]=((type64) t2)&bot58bits;
	t1=8*(type128)x[6]*x[8]  +4*((type128)x[0]*x[5]+(type128)x[1]*x[4]+(type128)x[2]*x[3]+(type128)x[7]*x[7])  +(t2>>58);
	z[5]=((type64) t1)&bot58bits;
	t2=8*(type128)x[7]*x[8]  +4*((type128)x[0]*x[6]+(type128)x[1]*x[5]+(type128)x[2]*x[4])  +(type128)x[3]*(2*x[3])  +(t1>>58);
	z[6]=((type64) t2)&bot58bits;
	t1=4*((type128)x[0]*x[7]+(type128)x[1]*x[6]+(type128)x[2]*x[5]+(type128)x[3]*x[4]+(type128)x[8]*x[8])  +(t2>>58);
	z[7]=((type64) t1)&bot58bits;
	st0+=(t1>>58);
	z[8]=((type64)st0)&bot58bits;
	z[0]+=2*(type64)(st0>>58);
}


// z=x*y - Granger's method

void gmul(type64 *x,type64 *y,type64 *z)
{
	type128 t0=(type128)x[0]*y[0] + 
		(type128)x[1]*y[1] + 
		(type128)x[2]*y[2] +
		(type128)x[3]*y[3] +
		(type128)x[4]*y[4];
	type128 t2,t3;
	type128 t5=(type128)x[5]*y[5];
	type128 t6=(type128)x[6]*y[6];
	type128 t7=(type128)x[7]*y[7];
	type128 t8=(type128)x[8]*y[8];
	type128 t1=t5+t6+t7+t8;
	t2=t0+t1-(type128)(x[0]-x[8])*(y[0]-y[8])-(type128)(x[1]-x[7])*(y[1]-y[7])
		-(type128)(x[2]-x[6])*(y[2]-y[6])-(type128)(x[3]-x[5])*(y[3]-y[5]);
	t0+=4*t1;

	type128 st1=((type64) t2)&bot58bits;

	t3=t0-(type128)(x[4]-2*x[5])*(y[4]-2*y[5])-(type128)(x[3]-2*x[6])*(y[3]-2*y[6])
	-(type128)(x[2]-2*x[7])*(y[2]-2*y[7])-(type128)(x[1]-2*x[8])*(y[1]-2*y[8])+2*(t2>>58);
	z[0]=((type64) t3)&bot58bits;

	t0-=2*t5;
	t2=t0-(type128)(x[0]-x[1])*(y[0]-y[1])-(type128)(x[4]-2*x[6])*(y[4]-2*y[6])
		-(type128)(x[2]-2*x[8])*(y[2]-2*y[8])-(type128)(x[3]-2*x[7])*(y[3]-2*y[7])+(t3>>58);
	z[1]=((type64) t2)&bot58bits;

	t0-=t5;
	t3=t0-(type128)(x[0]-x[2])*(y[0]-y[2])-(type128)(x[5]-2*x[6])*(y[5]-2*y[6])
		-(type128)(x[3]-2*x[8])*(y[3]-2*y[8])-(type128)(x[4]-2*x[7])*(y[4]-2*y[7])+(t2>>58);
	z[2]=((type64) t3)&bot58bits;

	t0-=2*t6;
	t2=t0-(type128)(x[0]-x[3])*(y[0]-y[3])-(type128)(x[1]-x[2])*(y[1]-y[2])
		-(type128)(x[4]-2*x[8])*(y[4]-2*y[8])-(type128)(x[5]-2*x[7])*(y[5]-2*y[7])+(t3>>58);
	z[3]=((type64) t2)&bot58bits;

	t0-=t6;
	t3=t0-(type128)(x[0]-x[4])*(y[0]-y[4])-(type128)(x[1]-x[3])*(y[1]-y[3])
		-(type128)(x[5]-2*x[8])*(y[5]-2*y[8])-(type128)(x[6]-2*x[7])*(y[6]-2*y[7])+(t2>>58);
	z[4]=((type64) t3)&bot58bits;

	t0-=2*t7;
	t2=t0-(type128)(x[0]-x[5])*(y[0]-y[5])-(type128)(x[1]-x[4])*(y[1]-y[4])
		-(type128)(x[2]-x[3])*(y[2]-y[3])-(type128)(x[6]-2*x[8])*(y[6]-2*y[8])+(t3>>58);
	z[5]=((type64) t2)&bot58bits;

	t0-=t7;
	t3=t0-(type128)(x[0]-x[6])*(y[0]-y[6])-(type128)(x[1]-x[5])*(y[1]-y[5])
		-(type128)(x[2]-x[4])*(y[2]-y[4])-(type128)(x[7]-2*x[8])*(y[7]-2*y[8])+(t2>>58);
	z[6]=((type64) t3)&bot58bits;

	t0-=2*t8;
	t2=t0-(type128)(x[0]-x[7])*(y[0]-y[7])-(type128)(x[1]-x[6])*(y[1]-y[6])
		-(type128)(x[2]-x[5])*(y[2]-y[5])-(type128)(x[3]-x[4])*(y[3]-y[4])+(t3>>58);
	z[7]=((type64) t2)&bot58bits;

	st1+=(t2>>58);
	z[8]=((type64)st1)&bot58bits;
	z[0]+=2*(type64)(st1>>58);
}

// Inverse x = 1/x = x^(p-2) mod p

int ginv(type64 *x)
{
	type64 x127[9],w[9],t[9],z[9];
	gsqr(x,x127);       // x127=x^2
	gmul(x127,x,t);     // t=x^3
	gsqr(t,x127);       // x127=x^6
	gmul(x127,x,w);     // w=x^7
	gsqr(w,x127);       // 
	gsqr(x127,t);  
	gsqr(t,x127);       // x127=x^56
	gcopy(x127,t);		// t=x^56
	gmul(w,t,x127);     // x127=x^63    
	gsqr(x127,t);
	gmul(t,x,x127);     // x127=x^127

	gsqr(x127,t);
	gmul(t,x,z);        // z=x^255

	gcopy(z,w);
	for (int i=0;i<4;i++)
	{
		gsqr(z,t);
		gsqr(t,z);
	}         
	gmul(z,w,t);        // z=z16       
  
	gcopy(t,w);
	for (int i=0;i<8;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}    
	gmul(t,w,z);        // z=z32      

	gcopy(z,w);
	for (int i=0;i<16;i++)
	{
		gsqr(z,t);
		gsqr(t,z);
	}    
	gmul(z,w,t);        // z=z64      

	gcopy(t,w);
	for (int i=0;i<32;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}    
	gmul(t,w,z);        // z=z128     

	gcopy(z,w);
	for (int i=0;i<64;i++)
	{
		gsqr(z,t);
		gsqr(t,z);
	}    
	gmul(z,w,t);        // z=z256       

	gcopy(t,w);
	for (int i=0;i<128;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}    
	gmul(t,w,z);      // z=z512        

	gsqr(z,t);
	gsqr(t,z);
	gsqr(z,t);
	gsqr(t,z);
	gsqr(z,t);
	gsqr(t,z);
	gsqr(z,t);
	gmul(t,x127,z);
	gsqr(z,t);
	gsqr(t,z);
	gmul(z,x,t);
	gcopy(t,x);
}

// Point Structure

typedef struct {
type64 x[9];
type64 y[9];
type64 z[9];
type64 t[9];
} ECp;

// P+=P

void double_1(ECp *P)
{
	type64 a[9],b[9],e[9],f[9],g[9],h[9];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gcopy(P->t,e);
	gmul2(e);
	gadd(a,b,g);
	gcopy(g,f); f[0]-=2;
	gsub(a,b,h);
	gmul(e,f,P->x); 
	gmul(g,h,P->y); 
	gsqr(g,P->z);
	gsb2(g,P->z);
	gmul(e,h,P->t);  
	scr(P->z);
}

void double_2(ECp *P)
{
	type64 a[9],b[9],c[9],e[9],f[9],g[9],h[9];
	gsqr(P->x,a);
	gsqr(P->y,b); 
	gsqr2(P->z,c); //gsqr(P->z,c); gmul2(c); 
	gadd(P->x,P->y,g); gsqr(g,e); gdec(a,e); gdec(b,e); 
	gadd(a,b,g); 
	gsub(g,c,f); 
	gsub(a,b,h); 
	gmul(e,f,P->x); 
	gmul(g,h,P->y);
	gmul(f,g,P->z); 
}

void double_3(ECp *P)
{
	type64 a[9],b[9],c[9],e[9],f[9],g[9],h[9];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gsqr2(P->z,c); //gsqr(P->z,c); gmul2(c);
	gadd(P->x,P->y,g); gsqr(g,e); gdec(a,e); gdec(b,e); 
	gadd(a,b,g);
	gsub(g,c,f);
	gsub(a,b,h);
	gmul(e,f,P->x); 
	gmul(g,h,P->y); 
	gmul(f,g,P->z);
	gmul(e,h,P->t); 
}

//P+=Q;

void add_1(ECp *Q,ECp *P)
{
	type64 a[9],b[9],c[9],d[9],e[9],f[9],g[9],h[9];
	gmul(P->x,Q->x,a);
	gmul(P->y,Q->y,b);
	gmul(P->t,Q->t,c);
	gadd(P->z,c,f);  /* reversed sign as d is negative */
	gsub(P->z,c,g);
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); gmul(c,d,e); gdec(a,e); gdec(b,e);
	gmul(e,f,P->x); 
	gmul(g,h,P->y); 
	gmul(f,g,P->z); 
	gmul(e,h,P->t); 
}

void add_2(ECp *Q,ECp *P)
{
	type64 a[9],b[9],c[9],d[9],e[9],f[9],g[9],h[9];
	gmul(P->x,Q->x,a);
	gmul(P->y,Q->y,b);
	gmul(P->t,Q->t,c);
	gmul(P->z,Q->z,d);
	gadd(d,c,f);  /* reversed sign as d is negative */
	gsub(d,c,g);
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); gmul(c,d,e); gdec(a,e); gdec(b,e);
	gmul(e,f,P->x); 
	gmul(g,h,P->y); 
	gmul(f,g,P->z); 
}

//P=0

void inf(ECp *P)
{
	for (int i=0;i<=8;i++)
		P->x[i]=P->y[i]=P->z[i]=P->t[i]=0;
	P->y[0]=P->z[0]=1;
}

// Initialise P

void init(type64 *x,type64 *y,ECp *P)
{
	for (int i=0;i<=8;i++)
	{
		P->x[i]=x[i];
		P->y[i]=y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	gmul(x,y,P->t);
}

//P=Q

void copy(ECp *Q,ECp *P)
{
	for (int i=0;i<=8;i++)
	{
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		P->t[i]=Q->t[i];
	}
}

// P=-Q

void neg(ECp *Q,ECp *P)
{
	for (int i=0;i<=8;i++)
	{
		P->x[i]=-Q->x[i]; 
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		P->t[i]=-Q->t[i]; 
	}
}

/* Make Affine */

void norm(ECp *P)
{
	type64 w[9],t[9];
	gcopy(P->z,w);
	ginv(w);
	gmul(P->x,w,t); scr(t); gcopy(t,P->x);
	gmul(P->y,w,t); scr(t); gcopy(t,P->y);
	gmul(P->z,w,t); scr(t); gcopy(t,P->z);
	gmul(P->t,w,t); scr(t); gcopy(t,P->t);
}

/* Precomputation */

void precomp(ECp *P,ECp W[])
{
	inf(&W[0]);
	copy(P,&W[1]); gmuli(W[1].t,376014);
	copy(P,&W[2]); double_1(&W[2]);
	copy(&W[2],&W[3]); add_1(&W[1],&W[3]);

	copy(&W[2],&W[4]); double_3(&W[4]);
	copy(&W[4],&W[5]); add_1(&W[1],&W[5]);
	copy(&W[3],&W[6]); double_3(&W[6]);
	copy(&W[6],&W[7]); add_1(&W[1],&W[7]);
	copy(&W[4],&W[8]); double_3(&W[8]);
	copy(&W[8],&W[9]); add_1(&W[1],&W[9]);

	copy(&W[5],&W[10]); double_3(&W[10]);
	copy(&W[10],&W[11]); add_1(&W[1],&W[11]);

	copy(&W[6],&W[12]); double_3(&W[12]);
	copy(&W[12],&W[13]); add_1(&W[1],&W[13]);

	copy(&W[7],&W[14]); double_3(&W[14]);
	copy(&W[14],&W[15]); add_1(&W[1],&W[15]);

	copy(&W[8],&W[16]); double_3(&W[16]);
	copy(&W[16],&W[17]); add_1(&W[1],&W[17]);

	copy(&W[9],&W[18]); double_3(&W[18]);
	copy(&W[18],&W[19]); add_1(&W[1],&W[19]);

	copy(&W[10],&W[20]); double_3(&W[20]);
	copy(&W[20],&W[21]); add_1(&W[1],&W[21]);

	copy(&W[11],&W[22]); double_3(&W[22]);
	copy(&W[22],&W[23]); add_1(&W[1],&W[23]);

	copy(&W[12],&W[24]); double_3(&W[24]);
	copy(&W[24],&W[25]); add_1(&W[1],&W[25]);

	copy(&W[13],&W[26]); double_3(&W[26]);
	copy(&W[26],&W[27]); add_1(&W[1],&W[27]);

	copy(&W[14],&W[28]); double_3(&W[28]);
	copy(&W[28],&W[29]); add_1(&W[1],&W[29]);

	copy(&W[15],&W[30]); double_3(&W[30]);
	copy(&W[30],&W[31]); add_1(&W[1],&W[31]);

	copy(&W[16],&W[32]); double_3(&W[32]);

/* premultiply t parameter by curve constant */

	gmuli(W[2].t,376014);
	gmuli(W[3].t,376014);
	gmuli(W[4].t,376014);
	gmuli(W[5].t,376014);
	gmuli(W[6].t,376014);
	gmuli(W[7].t,376014);
	gmuli(W[8].t,376014);
	gmuli(W[9].t,376014);
	gmuli(W[10].t,376014);
	gmuli(W[11].t,376014);
	gmuli(W[12].t,376014);
	gmuli(W[13].t,376014);
	gmuli(W[14].t,376014);
	gmuli(W[15].t,376014);
	gmuli(W[16].t,376014);
	gmuli(W[17].t,376014);
	gmuli(W[18].t,376014);
	gmuli(W[19].t,376014);
	gmuli(W[20].t,376014);
	gmuli(W[21].t,376014);
	gmuli(W[22].t,376014);
	gmuli(W[23].t,376014);
	gmuli(W[24].t,376014);
	gmuli(W[25].t,376014);
	gmuli(W[26].t,376014);
	gmuli(W[27].t,376014);
	gmuli(W[28].t,376014);
	gmuli(W[29].t,376014);
	gmuli(W[30].t,376014);
	gmuli(W[31].t,376014);
	gmuli(W[32].t,376014);
}

/* Windows of width 6 */

void window(ECp *Q,ECp *P)
{
	double_2(P);
	double_2(P);
	double_2(P);
	double_2(P);
	double_2(P);
	double_3(P);
	add_2(Q,P);
}

void output(ECp *P)
{
	cout << "x[0]= " << hex << P->x[0] << endl;
	cout << "x[1]= " << hex << P->x[1] << endl;
	cout << "x[2]= " << hex << P->x[2] << endl;
	cout << "x[3]= " << hex << P->x[3] << endl;
	cout << "x[4]= " << hex << P->x[4] << endl;
	cout << "x[5]= " << hex << P->x[5] << endl;
	cout << "x[6]= " << hex << P->x[6] << endl;
	cout << "x[7]= " << hex << P->x[7] << endl;
	cout << "x[8]= " << hex << P->x[8] << endl;
	cout << endl;

	cout << "y[0]= " << hex << P->y[0] << endl;
	cout << "y[1]= " << hex << P->y[1] << endl;
	cout << "y[2]= " << hex << P->y[2] << endl;
	cout << "y[3]= " << hex << P->y[3] << endl;
	cout << "y[4]= " << hex << P->y[4] << endl;
	cout << "y[5]= " << hex << P->y[5] << endl;
	cout << "y[6]= " << hex << P->y[6] << endl;
	cout << "y[7]= " << hex << P->y[7] << endl;
	cout << "y[8]= " << hex << P->y[8] << endl;
	cout << endl;
}

/* Point Multiplication - exponent is 519 bits */

void mul(int *w,ECp *P)
{
	ECp W[33],Q;
	precomp(P,W);

	copy(&W[w[86]],P);  
	for (int i=85;i>=0;i--)
	{
		if (w[i]>=0) copy(&W[w[i]],&Q);
		else         neg(&W[-w[i]],&Q);
		window(&Q,P);
	}
	norm(P); 
}

//#define TEST  /* define to multiply by group order */

int main()
{
	uint64_t bef,aft;
	int i,w[87];
	int ii,lpz=10000;
	type64 xs[9],ys[9];
	ECp P;

/* Base point on Curve (from SafeCurves Website) */

	xs[0]=0x2A940A2F19BA6CLL;
	xs[1]=0x3EC4CD920E2A8CLL;
	xs[2]=0x1D568FC99C6059DLL;
	xs[3]=0x3331C90D2C6BA52LL;
	xs[4]=0xC6203913F6ECC5LL;
	xs[5]=0x1B2063B22FCF270LL;
	xs[6]=0x2878A3BFD9F42FCLL;
	xs[7]=0x6277E432C8A5ACLL;
	xs[8]=0x752CB45C48648BLL;

	ys[0]=0xcLL;
	ys[1]=0x0LL;
	ys[2]=0x0LL;
	ys[3]=0x0LL;
	ys[4]=0x0LL;
	ys[5]=0x0LL;
	ys[6]=0x0LL;
	ys[7]=0x0LL;
	ys[8]=0x0LL;


// random multiplier arranged in signed windows of size 6
#ifndef TEST
w[0]= 7; w[1]= -11; w[2]= -31; w[3]= 20; w[4]= -3; w[5]= -1;
w[6]= -7; w[7]= 17; w[8]= 26; w[9]= -30; w[10]= 29; w[11]= 5; 
w[12]= -20; w[13]= -8; w[14]= 7; w[15]= -22; w[16]= 1; w[17]= 8;
w[18]= -21; w[19]= -14; w[20]= -20; w[21]= -6; w[22]= 23; w[23]= 27;
w[24]= 4; w[25]= -22; w[26]= -25; w[27]= 24; w[28]= 13; w[29]= -21;
w[30]= 25; w[31]= -1; w[32]= -7; w[33]= -16; w[34]= -1; w[35]= 3;
w[36]= 17; w[37]= -10; w[38]= 5; w[39]= 8; w[40]= 17; w[41]= -7;
w[42]= 22; w[43]= -8; w[44]= 27; w[45]= 23; w[46]= 28; w[47]= -20;
w[48]= -3; w[49]= -30; w[50]= -14; w[51]= 18; w[52]= 27; w[53]= -15;
w[54]= 29; w[55]= 31; w[56]= 7; w[57]= -14; w[58]= 29; w[59]= -17; 
w[60]= -14; w[61]= -20; w[62]= -22; w[63]= 5; w[64]= 27; w[65]= -5;
w[66]= 2; w[67]= 16; w[68]= 24; w[69]= -1; w[70]= 5; w[71]= -30;
w[72]= -18; w[73]= 22; w[74]= 31; w[75]= 29; w[76]= 20; w[77]= -15;
w[78]= 21; w[79]= -3; w[80]= 16; w[81]= -1; w[82]= 12; w[83]= -2;
w[84]= -28; w[85]= 28; w[86]= 4;

// Group Order - for testing
#else
w[0]= -21; w[1]= -10; w[2]= 1; w[3]= 6; w[4]= -11; w[5]= 24;
w[6]= 3; w[7]= 9; w[8]= -22; w[9]= 4; w[10]= 20; w[11]= 17;
w[12]= 31; w[13]= -4; w[14]= -23; w[15]= -25; w[16]= 23; w[17]= 17;
w[18]= 12; w[19]= -10; w[20]= -4; w[21]= 20; w[22]= -16; w[23]= 16;
w[24]= 5; w[25]= -5; w[26]= -24; w[27]= 24; w[28]= -17; w[29]= -29;
w[30]= -20; w[31]= 14; w[32]= -9; w[33]= 24; w[34]= 8; w[35]= -1;
w[36]= 7; w[37]= 29; w[38]= -28; w[39]= -14; w[40]= -9; w[41]= 23;
w[42]= 17; w[43]= -1; w[44]= 0; w[45]= 0; w[46]= 0; w[47]= 0;
w[48]= 0; w[49]= 0; w[50]= 0; w[51]= 0; w[52]= 0; w[53]= 0;
w[54]= 0; w[55]= 0; w[56]= 0; w[57]= 0; w[58]= 0; w[59]= 0; 
w[60]= 0; w[61]= 0; w[62]= 0; w[63]= 0; w[64]= 0; w[65]= 0;
w[66]= 0; w[67]= 0; w[68]= 0; w[69]= 0; w[70]= 0; w[71]= 0;
w[72]= 0; w[73]= 0; w[74]= 0; w[75]= 0; w[76]= 0; w[77]= 0;
w[78]= 0; w[79]= 0; w[80]= 0; w[81]= 0; w[82]= 0; w[83]= 0;
w[84]= 0; w[85]= 0; w[86]= 8;
#endif

	bef=rdtsc();
	for (ii=0;ii<lpz;ii++)
	{
		init(xs,ys,&P);
		mul(w,&P);
	}
	aft=rdtsc();
	cout << "Clock cycles= " << (aft-bef)/(lpz) << endl;

	output(&P);

	return 0;
}



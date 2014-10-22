// Test program for ws521 scalar point multiplication
// This is NIST standard Weierstrass curve p-521 
// Fully Tested and debugged
// g++ -O3 ws521.cpp -o ws521
// M.Scott 18/9/2014

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

// fused operations

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

//w=w-x-y
void gtsb(type64 *x,type64 *y,type64 *w)
{
	w[0]-=x[0]+y[0];
	w[1]-=x[1]+y[1];
	w[2]-=x[2]+y[2];
	w[3]-=x[3]+y[3];
	w[4]-=x[4]+y[4];
	w[5]-=x[5]+y[5];
	w[6]-=x[6]+y[6];
	w[7]-=x[7]+y[7];
	w[8]-=x[8]+y[8];	
}

//w=w-2x-y
void gtsb2(type64 *x,type64 *y,type64 *w)
{
	w[0]-=2*x[0]+y[0];
	w[1]-=2*x[1]+y[1];
	w[2]-=2*x[2]+y[2];
	w[3]-=2*x[3]+y[3];
	w[4]-=2*x[4]+y[4];
	w[5]-=2*x[5]+y[5];
	w[6]-=2*x[6]+y[6];
	w[7]-=2*x[7]+y[7];
	w[8]-=2*x[8]+y[8];
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

// w=2(x-y)
void g2sb(type64 *x,type64 *y,type64 *w)
{
	w[0]=2*(x[0]-y[0]);
	w[1]=2*(x[1]-y[1]);
	w[2]=2*(x[2]-y[2]);
	w[3]=2*(x[3]-y[3]);
	w[4]=2*(x[4]-y[4]);
	w[5]=2*(x[5]-y[5]);
	w[6]=2*(x[6]-y[6]);
	w[7]=2*(x[7]-y[7]);
	w[8]=2*(x[8]-y[8]);
}


// w=3(x+y)
void g3ad(type64 *x,type64 *y,type64 *w)
{
	w[0]=3*(x[0]+y[0]);
	w[1]=3*(x[1]+y[1]);
	w[2]=3*(x[2]+y[2]);
	w[3]=3*(x[3]+y[3]);
	w[4]=3*(x[4]+y[4]);
	w[5]=3*(x[5]+y[5]);
	w[6]=3*(x[6]+y[6]);
	w[7]=3*(x[7]+y[7]);
	w[8]=3*(x[8]+y[8]);
}


// w=4*w-x
void g4sb(type64 *x,type64 *w)
{
	w[0]=4*w[0]-x[0];
	w[1]=4*w[1]-x[1];
	w[2]=4*w[2]-x[2];
	w[3]=4*w[3]-x[3];
	w[4]=4*w[4]-x[4];
	w[5]=4*w[5]-x[5];
	w[6]=4*w[6]-x[6];
	w[7]=4*w[7]-x[7];
	w[8]=4*w[8]-x[8];
}

// w*=4

void gmul4(type64 *w)
{
	w[0]*=4;
	w[1]*=4;
	w[2]*=4;
	w[3]*=4;
	w[4]*=4;
	w[5]*=4;
	w[6]*=4;
	w[7]*=4;
	w[8]*=4;
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

// w-=8*x
void gsb8(type64 *x,type64 *w)
{
	w[0]-=8*x[0];
	w[1]-=8*x[1];
	w[2]-=8*x[2];
	w[3]-=8*x[3];
	w[4]-=8*x[4];
	w[5]-=8*x[5];
	w[6]-=8*x[6];
	w[7]-=8*x[7];
	w[8]-=8*x[8];
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

// z=x^2
// Note t0=r8|r9, t1=r10|r11, t2=r12|r13, t3=r14|r15

void gsqr(type64 *x,type64 *z)
{
	type128 t0,t1,t2;
	type64 st0;
	t1=2*((type128)x[0]*x[8]+(type128)x[1]*x[7]+(type128)x[2]*x[6]+(type128)x[3]*x[5])+(type128)x[4]*x[4];
	st0=((type64) t1)&bot58bits;
	t2=4*((type128)x[1]*x[8]+(type128)x[2]*x[7]+(type128)x[3]*x[6]+(type128)x[4]*x[5])+(type128)x[0]*x[0]+2*(t1>>58);
	z[0]=((type64) t2)&bot58bits;
	t1=4*((type128)x[2]*x[8]+(type128)x[3]*x[7]+(type128)x[4]*x[6])+2*((type128)x[0]*x[1]+(type128)x[5]*x[5])+(t2>>58);
	z[1]=((type64) t1)&bot58bits;
	t2=4*((type128)x[3]*x[8]+(type128)x[4]*x[7]+(type128)x[5]*x[6])+2*(type128)x[0]*x[2]+(type128)x[1]*x[1]+(t1>>58);
	z[2]=((type64) t2)&bot58bits;
	t1=4*((type128)x[4]*x[8]+(type128)x[5]*x[7])+2*((type128)x[0]*x[3]+(type128)x[1]*x[2]+(type128)x[6]*x[6])+(t2>>58);
	z[3]=((type64) t1)&bot58bits;
	t2=4*((type128)x[5]*x[8]+(type128)x[6]*x[7])+2*((type128)x[0]*x[4]+(type128)x[1]*x[3])+(type128)x[2]*x[2]+(t1>>58);
	z[4]=((type64) t2)&bot58bits;
	t1=4*(type128)x[6]*x[8]+2*((type128)x[0]*x[5]+(type128)x[1]*x[4]+(type128)x[2]*x[3]+(type128)x[7]*x[7])+(t2>>58);
	z[5]=((type64) t1)&bot58bits;
	t2=4*(type128)x[7]*x[8]+2*((type128)x[0]*x[6]+(type128)x[1]*x[5]+(type128)x[2]*x[4])+(type128)x[3]*x[3]+(t1>>58);
	z[6]=((type64) t2)&bot58bits;
	t1=2*((type128)x[0]*x[7]+(type128)x[1]*x[6]+(type128)x[2]*x[5]+(type128)x[3]*x[4]+(type128)x[8]*x[8])+(t2>>58);
	z[7]=((type64) t1)&bot58bits;
	st0+=(type64)(t1>>58);
	z[8]=st0&bot58bits;
	z[0]+=2*(st0>>58);
}

// z=x*y
// Note t0=r8|r9, t1=r10|r11, t2=r12|r13, t3=r14|r15

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

	type64 st1=((type64) t2)&bot58bits;

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

	st1+=((type64)(t2>>58));
	z[8]=st1&bot58bits;
	z[0]+=2*(st1>>58);
}

//
// Inverse x = 1/x = x^(p-2) mod p
// 13 muls, 520 sqrs
//

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
	gmul(z,w,t);		// z=z256       

	gcopy(t,w);
	for (int i=0;i<128;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}    
	gmul(t,w,z);		// z=z512        

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
int inf;
} ECp;

// P=0

void inf(ECp *P)
{
	P->inf=1;
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
	P->inf=0;
}

// P=Q

void copy(ECp *Q,ECp *P)
{
	for (int i=0;i<=8;i++)
	{
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
	}
	P->inf=Q->inf;
}

// P=-Q

void neg(ECp *Q,ECp *P)
{
	for (int i=0;i<=8;i++)
	{
		P->x[i]=Q->x[i]; 
		P->y[i]=-Q->y[i];
		P->z[i]=Q->z[i]; 
	}
	P->inf=Q->inf;
}

/* Make Affine */

void norm(ECp *P)
{
	type64 iz2[9],iz3[9],w[9],t[9];
	if (P->inf) return;
	gcopy(P->z,w);
	ginv(w);
	gsqr(w,iz2);
	gmul(iz2,w,iz3);
	gmul(P->x,iz2,t); scr(t); gcopy(t,P->x);
	gmul(P->y,iz3,t); scr(t); gcopy(t,P->y);
	gmul(P->z,w,t); scr(t); gcopy(t,P->z);
}

void doubl(ECp *P)
{
	type64 r0[9],r1[9],r2[9],r3[9];

	gsqr(P->z,r0);
	gsqr(P->y,r1);
	gmul(P->x,r1,r2);
	gadd(P->y,P->z,r3);
	g3ad(P->x,r0,P->y);
	gdec(r0,P->x);
	gmul(P->x,P->y,P->z);
	gsqr(P->z,P->x);
	gsb8(r2,P->x);
	scr(P->x);
	g4sb(P->x,r2);
	gmul(P->z,r2,P->y);
	gsqr(r1,r2);
	gsb8(r2,P->y);
	scr(P->y);
	gsqr(r3,P->z);
	gtsb(r0,r1,P->z);
	scr(P->z);
}

// P+=Q, Q is affine

void add_a(ECp *Q,ECp *P)
{
	type64 r0[9],r1[9],r2[9],r3[9],r4[9];
// if zs=1 then r4=xt, r1=yt
	gsqr(P->z,r3);  // not if its one
	gmul(Q->x,r3,r4); // ditto
	gmul(Q->y,r3,r2); //ditto
	gmul(r2,P->z,r1); // not if its one

	gdec(P->x,r4);  //w1
	g2sb(r1,P->y,r2); //w8
	gsqr(r4,r0);

// if zs=1 then new zs=2.r4
	gadd(P->z,r4,r1);
	gsqr(r1,P->z);
	gtsb(r3,r0,P->z);

	scr(P->z);

	gmul4(r0);
	gmul(r4,r0,r1);
	gmul(r1,P->y,r3);
 	gmul(r0,P->x,P->y);
	gsqr(r2,P->x);
	gtsb2(P->y,r1,P->x);
	scr(P->x);

	gsub(P->y,P->x,r0);
	gmul(r0,r2,P->y);
	gsb2(r3,P->y);
	scr(P->y);
}

// P+=Q, Q is projective

void add_p(ECp *Q,ECp *P)
{
	type64 z1z1[9],z2z2[9],u1[9],u2[9],s1[9],s2[9],h[9],i[9],j[9];
	gsqr(P->z,z1z1);
	gsqr(Q->z,z2z2);
	gmul(P->x,z2z2,u1);
	gmul(Q->x,z1z1,u2);
	gmul(P->y,Q->z,h);
	gmul(h,z2z2,s1);
	gmul(Q->y,P->z,h);
	gmul(h,z1z1,s2);
	gsub(u2,u1,h);
	gcopy(h,j); gmul2(j); gsqr(j,i);
	gmul(h,i,j);

	gmul(u1,i,u2);
	g2sb(s2,s1,u1);
	gsqr(u1,i);
	gsub(i,j,P->x);
	gsb2(u2,P->x);
	scr(P->x);

	gmul(s1,j,i);
	gsub(u2,P->x,j);
	gmul(j,u1,P->y);
	gsb2(i,P->y);
	scr(P->y);

	gadd(P->z,Q->z,i);
	gsqr(i,j);

	gtsb(z1z1,z2z2,j);
	gmul(h,j,P->z);
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

/* Precomputation */

void precomp(ECp *P,ECp W[])
{
	ECp Q;
	copy(P,&Q);
	doubl(&Q);
	copy(P,&W[0]); //P
	
	for (int i=1;i<32;i++)
	{
		copy(&W[i-1],&W[i]);
		add_p(&Q,&W[i]);
	}
}

/* Windows of width 6 */

void window(ECp *Q,ECp *P)
{
	doubl(P);
	doubl(P);
	doubl(P);
	doubl(P);
	doubl(P);
	doubl(P);
	add_p(Q,P);
}

/* Point Multiplication - exponent is 521 bits */

void mul(int *w,ECp *P)
{
	int j;
	ECp W[32],Q;
	precomp(P,W);
	copy(&W[(w[87]-1)/2],P);  
	for (int i=86;i>=0;i--)
	{
		if (w[i]>=0) {j=(w[i]-1)/2; copy(&W[j],&Q); /*printf("j= %d\n",j);*/}
		else         {j=(-w[i]-1)/2; neg(&W[j],&Q); /*printf("j= -%d\n",j);*/}
		window(&Q,P);
	}

	norm(P); 
}

//#define TEST  /* define to multiply by group order */

int main()
{
	uint64_t bef,aft;
	int i,w[88];
	int ii,lpz=10000;
	ECp P;
	type64 xs[9],ys[9];

/* Base point on NIST Curve */

	xs[0]=0x17E7E31C2E5BD66LL;
	xs[1]=0x22CF0615A90A6FELL;
	xs[2]=0x127A2FFA8DE334LL;
	xs[3]=0x1DFBF9D64A3F877LL;
	xs[4]=0x6B4D3DBAA14B5ELL;
	xs[5]=0x14FED487E0A2BD8LL;
	xs[6]=0x15B4429C6481390LL;
	xs[7]=0x3A73678FB2D988ELL;
	xs[8]=0xC6858E06B70404LL;

	ys[0]=0xBE94769FD16650LL;
	ys[1]=0x31C21A89CB09022LL;
	ys[2]=0x39013FAD0761353LL;
	ys[3]=0x2657BD099031542LL;
	ys[4]=0x3273E662C97EE72LL;
	ys[5]=0x1E6D11A05EBEF45LL;
	ys[6]=0x3D1BD998F544495LL;
	ys[7]=0x3001172297ED0B1LL;
	ys[8]=0x11839296A789A3BLL;

#ifndef TEST
w[0]= 13; w[1]= 29; w[2]= -25; w[3]= -39; w[4]= -55; w[5]= 53; w[6]= -35; w[7]= 63; w[8]= -53; 
w[9]= -9; w[10]= 43; w[11]= -15; w[12]= 61; w[13]= -63; w[14]= -33; w[15]= -13; w[16]= 33; 
w[17]= -47; w[18]= -33; w[19]= -7; w[20]= -25; w[21]= 21; w[22]= -53; w[23]= -35; w[24]= -39; 
w[25]= -25; w[26]= -23; w[27]= -63; w[28]= -59; w[29]= -39; w[30]= 45; w[31]= -5; w[32]= 13; 
w[33]= -11; w[34]= 7; w[35]= 63; w[36]= 27; w[37]= -5; w[38]= -41; w[39]= 61; w[40]= -31; 
w[41]= -17; w[42]= 23; w[43]= -39; w[44]= 15; w[45]= 27; w[46]= -27; w[47]= 55; w[48]= 41; 
w[49]= -13; w[50]= 59; w[51]= -41; w[52]= 31; w[53]= 41; w[54]= 7; w[55]= 3; w[56]= 59; 
w[57]= -63; w[58]= 59; w[59]= 53; w[60]= -13; w[61]= -23; w[62]= 33; w[63]= 63; w[64]= 13; 
w[65]= -13; w[66]= -59; w[67]= 1; w[68]= -1; w[69]= 9; w[70]= -59; w[71]= 17; w[72]= -59; 
w[73]= 59; w[74]= 41; w[75]= 59; w[76]= 25; w[77]= -41; w[78]= 9; w[79]= 7; w[80]= -31; 
w[81]= -11; w[82]= 25; w[83]= 33; w[84]= 29; w[85]= 59; w[86]= -49; w[87]= 1;
#else
w[0]= -55; w[1]= -47; w[2]= -57; w[3]= 15; w[4]= -47; w[5]= 59; w[6]= 49; w[7]= 45; w[8]= 47; 
w[9]=45; w[10]= 43; w[11]= 43; w[12]= 7; w[13]= 49; w[14]= -39; w[15]= -29; w[16]= -7; 
w[17]= -25; w[18]= 29; w[19]= 45; w[20]= -5; w[21]= 1; w[22]= 29; w[23]= 41; w[24]= -55; 
w[25]= 29; w[26]= -49; w[27]= 19; w[28]= -63; w[29]= -15; w[30]= 61; w[31]= 31; w[32]= 43; 
w[33]= 25; w[34]= 57; w[35]= 11; w[36]= -1; w[37]= -49; w[38]= 57; w[39]= -31; w[40]= -57; 
w[41]= 7; w[42]= -27; w[43]= 63; w[44]= 63; w[45]= 63; w[46]= 63; w[47]= 63; w[48]= 63; 
w[49]= 63; w[50]= 63; w[51]= 63; w[52]= 63; w[53]= 63; w[54]= 63; w[55]= 63; w[56]= 63; 
w[57]= 63; w[58]= 63; w[59]= 63; w[60]= 63; w[61]= 63; w[62]= 63; w[63]= 63; w[64]= 63; 
w[65]= 63; w[66]= 63; w[67]= 63; w[68]= 63; w[69]= 63; w[70]= 63; w[71]= 63; w[72]= 63; 
w[73]= 63; w[74]= 63; w[75]= 63; w[76]= 63; w[77]= 63; w[78]= 63; w[79]= 63; w[80]= 63; 
w[81]=63; w[82]= 63; w[83]= 63; w[84]= 63; w[85]= 63; w[86]= -33; w[87]= 1;
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


